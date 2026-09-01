# ---------------------------------------------------------------------------- #
#                                  mtx_expr.jl                                 #
# ---------------------------------------------------------------------------- #
#
# Extracts a single gene's expression vector from a (optionally gzipped)
# Matrix Market (.mtx) sparse count matrix -- genes as rows, cells as
# columns, as commonly exported from Seurat/scanpy pipelines -- plus its
# accompanying features (gene names, one per row) and barcodes (cell ids,
# one per column) files. `include()`d by dynema_map.jl.
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

using CodecZlib
using DataFrames

"""
    open_maybe_gzip(path) -> IO

Opens `path` for reading, transparently decompressing it if it ends in
`.gz`.
"""
function open_maybe_gzip(path::AbstractString)
    io = open(path)
    return endswith(path, ".gz") ? GzipDecompressorStream(io) : io
end

"""
    read_labels(path) -> Vector{String}

Reads a features/barcodes file: one label per line, gzip-transparent. If a
line has multiple tab-separated fields (e.g. a 10x-style features.tsv with
Ensembl id + gene symbol columns), the *last* field is used.
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
        error("--expr-prefix '$prefix' is missing expected file(s):\n  " * join(missing_msgs, "\n  "))

    return mtx, features, barcodes

end

"""
    extract_gene_expression(; mtx, features, barcodes, gene, id_col="cell_id") -> DataFrame

Streams through `mtx` once, extracting the expression of `gene` (matched
against `features`, one entry per matrix row) across every cell (one entry
per matrix column, named from `barcodes`). Matrix entries not present in the
sparse file are zero. Returns a two-column DataFrame (`id_col`, `gene`),
shaped exactly like a plain `--expr` TSV would be.
"""
function extract_gene_expression(; mtx::AbstractString, features::AbstractString,
                                  barcodes::AbstractString, gene::AbstractString,
                                  id_col::AbstractString = "cell_id")

    section("Reading gene expression (Matrix Market)")
    bullet("matrix:    $mtx")
    bullet("features:  $features")
    bullet("barcodes:  $barcodes")

    gene_names = read_labels(features)
    row_matches = findall(==(gene), gene_names)
    isempty(row_matches) && error("Gene '$gene' not found in $features")
    length(row_matches) > 1 &&
        error("Gene '$gene' matches $(length(row_matches)) rows in $features; expected exactly one")
    target_row = row_matches[1]
    bullet("gene '$gene' found at row $target_row of $(length(gene_names))")

    cell_ids = read_labels(barcodes)
    n_cells = length(cell_ids)

    bullet("scanning for '$gene' across $n_cells cell(s)...")
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
    bullet("found $n_found nonzero entries for '$gene' across $n_cells cell(s)")

    df = DataFrame()
    df[!, id_col] = cell_ids
    df[!, gene] = expr

    return df

end
