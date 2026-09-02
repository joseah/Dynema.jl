# ---------------------------------------------------------------------------- #
#                          matrix_market_expression.jl                         #
# ---------------------------------------------------------------------------- #
#
# Core single-gene extraction from a Matrix Market (.mtx) sparse count
# matrix -- genes as rows, cells as columns, as commonly exported from
# Seurat/scanpy pipelines -- usable directly from Julia code (not just via
# bin/dynema-map) -- see `extract_gene_expression` below, the public entry
# point.
#
# Matrix Market coordinate files store nonzero (row, col, value) triples,
# typically sorted by column rather than by row/gene, with no separate index:
# there is no way to jump straight to one gene's entries the way tabix does
# for a VCF. So this still means one sequential pass over every line -- but
# a *streaming* one. Unlike a generic Matrix Market reader (e.g.
# MatrixMarket.jl's mmread), it never builds the full sparse matrix in
# memory: it keeps only the (typically few, since single-cell data is
# sparse) entries belonging to the requested row and discards the rest as it
# goes. For a multi-GB matrix this is still a full read of the file (time
# dominated by decompression + line parsing -- for hundreds of millions of
# entries, expect low-single-digit minutes, not seconds), but with
# essentially constant memory use rather than needing to materialize the
# whole matrix.

"""
    open_maybe_gzip(path) -> IO

Opens `path` for reading, transparently decompressing it if it ends in
`.gz`. Internal helper for `extract_gene_expression`/`read_labels`; not
exported.
"""
function open_maybe_gzip(path::AbstractString)
    io = open(path)
    return endswith(path, ".gz") ? GzipDecompressorStream(io) : io
end

"""
    read_labels(path) -> Vector{String}

Reads a features/barcodes file: one label per line, gzip-transparent. If a
line has multiple tab-separated fields (e.g. a 10x-style features.tsv with
Ensembl id + gene symbol columns), the *last* field is used. Internal
helper for `extract_gene_expression`; not exported.
"""
function read_labels(path::AbstractString)
    io = open_maybe_gzip(path)
    labels = try
        [String(split(line, '\t')[end]) for line in eachline(io)]
    finally
        close(io)
    end
    return labels
end

"""
    resolve_mtx_triplet(prefix) -> (mtx, features, barcodes)

Given a filename prefix (e.g. `"expr_0.05"`), resolves the three companion
files a Matrix Market export usually comes as: `<prefix>.mtx`,
`<prefix>.features`, `<prefix>.barcodes` -- each tried gzipped
(`<prefix>.mtx.gz`, ...) first, then plain. Errors up front, naming exactly
which file(s) are missing and every path that was tried for each, rather
than failing deep inside the (potentially minutes-long) matrix scan.
"""
function resolve_mtx_triplet(prefix::AbstractString)

    function find_one(suffix::AbstractString)
        candidates = ["$(prefix).$(suffix).gz", "$(prefix).$(suffix)"]
        idx = findfirst(isfile, candidates)
        found = idx === nothing ? nothing : candidates[idx]
        return found, candidates
    end

    mtx, mtx_tried           = find_one("mtx")
    features, features_tried = find_one("features")
    barcodes, barcodes_tried = find_one("barcodes")

    missing_msgs = String[]
    mtx === nothing      && push!(missing_msgs, "matrix (.mtx): none of $(join(mtx_tried, ", ")) exist")
    features === nothing && push!(missing_msgs, "features: none of $(join(features_tried, ", ")) exist")
    barcodes === nothing && push!(missing_msgs, "barcodes: none of $(join(barcodes_tried, ", ")) exist")

    isempty(missing_msgs) ||
        error("Matrix Market prefix '$prefix' is missing expected file(s):\n  " * join(missing_msgs, "\n  "))

    return mtx, features, barcodes

end

"""
    extract_gene_expression(; mtx, features, barcodes, gene, id_col="cell_id", verbose=true)

Core Matrix Market single-gene extraction: streams through `mtx` once,
extracting the expression of `gene` (matched against `features`, one entry
per matrix row) across every cell (one entry per matrix column, named from
`barcodes`), without ever loading the whole sparse matrix into memory (see
module-level comment above). Matrix entries not present in the sparse file
are zero.

Returns a `NamedTuple` with:
- `expr::DataFrame`: two columns (`id_col`, `gene`), shaped exactly like a plain
  `dynema-map --expr` TSV would be.
- `target_row`, `n_genes`: which row `gene` was found at, and how many genes/rows
  `features` has in total.
- `n_cells`: number of cells/columns (from `barcodes`).
- `n_found`: number of nonzero entries read for `gene`.

Set `verbose=false` to suppress the built-in progress `println`s -- e.g. for
a caller (like the `dynema-map` CLI) that wants to render its own progress
messages from the returned counts instead.
"""
function extract_gene_expression(; mtx::AbstractString, features::AbstractString,
                                  barcodes::AbstractString, gene::AbstractString,
                                  id_col::AbstractString = "cell_id",
                                  verbose::Bool = true)

    printlnv("Reading gene list from $features..."; verbose)
    gene_names = read_labels(features)
    row_matches = findall(==(gene), gene_names)
    isempty(row_matches) && error("Gene '$gene' not found in $features")
    length(row_matches) > 1 &&
        error("Gene '$gene' matches $(length(row_matches)) rows in $features; expected exactly one")
    target_row = row_matches[1]

    printlnv("Reading cell barcodes from $barcodes..."; verbose)
    cell_ids = read_labels(barcodes)
    n_cells = length(cell_ids)

    printlnv("Scanning $mtx for gene '$gene' (row $target_row of $(length(gene_names)))..."; verbose)
    io = open_maybe_gzip(mtx)
    expr = zeros(Float64, n_cells)
    n_found = 0
    dims_read = false

    try
        for line in eachline(io)

            startswith(line, "%") && continue # %%MatrixMarket / comment lines

            if !dims_read
                dims = split(line)
                length(dims) == 3 || error("Unexpected Matrix Market dimensions line in $mtx: '$line'")
                n_rows, n_cols, _ = parse.(Int, dims)
                n_rows == length(gene_names) ||
                    error("$mtx has $n_rows rows but $features has $(length(gene_names)) genes -- files don't match")
                n_cols == n_cells ||
                    error("$mtx has $n_cols columns but $barcodes has $n_cells barcodes -- files don't match")
                dims_read = true
                continue
            end

            sp1 = findfirst(==(' '), line)
            sp1 === nothing && continue
            row = parse(Int, SubString(line, 1, sp1 - 1))

            row == target_row || continue

            rest = SubString(line, sp1 + 1)
            sp2 = findfirst(==(' '), rest)
            col = parse(Int, SubString(rest, 1, sp2 - 1))
            val = parse(Float64, SubString(rest, sp2 + 1))
            expr[col] = val
            n_found += 1

        end
    finally
        close(io)
    end

    dims_read || error("$mtx has no data after its comment header")
    printlnv("Found $n_found nonzero entries for '$gene' across $n_cells cells"; verbose)

    df = DataFrame()
    df[!, id_col] = cell_ids
    df[!, gene] = expr

    return (; expr = df, target_row, n_genes = length(gene_names), n_cells, n_found)

end
