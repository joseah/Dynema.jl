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

Reads a barcodes file: one label per line, gzip-transparent. If a line has
multiple tab-separated fields, the *last* field is used. (Features files
are read with `read_feature_fields` instead, keeping every column so genes
can be matched by id or name.) Internal helper for
`extract_gene_expression`; not exported.
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
    read_feature_fields(path) -> Vector{Vector{String}}

Reads a features file as one vector of tab-separated fields per line,
gzip-transparent -- handling both single-column exports (gene names only)
and the 10x two/three-column layout (gene id, gene name[, modality]). Used
to match a gene specification against *every* identifier a row carries.
Internal helper for `extract_gene_expression`; not exported.
"""
function read_feature_fields(path::AbstractString)
    io = open_maybe_gzip(path)
    try
        return [String.(split(line, '\t')) for line in eachline(io)]
    finally
        close(io)
    end
end

"""
    match_feature_row(feats, gene; features_path = "the features file") -> Int

Finds `gene`'s row among parsed features rows (see
[`read_feature_fields`](@ref)). Two supported layouts:

  - single column (e.g. a Seurat export): gene names only -- `gene` must
    match that column;
  - 10x layout (gene_id, gene_name[, modality]): `gene` is searched against
    gene_name (column 2) first; if nothing matches, it is assumed to be a
    gene id and searched against column 1.

Ensembl id version suffixes are ignored; no match, or an ambiguous match, is
an error (`features_path` names the file in messages). Internal helper for
`extract_gene_expression` and the CLI's --check mode; not exported.
"""
function match_feature_row(feats::Vector{Vector{String}}, gene::AbstractString;
                           features_path::AbstractString = "the features file")

    g = stripver(gene)
    findmatches(col) = findall(row -> length(row) >= col && stripver(row[col]) == g, feats)
    tenx = any(row -> length(row) >= 2, feats)

    row_matches, searched = if tenx
        m = findmatches(2)
        isempty(m) ? (findmatches(1), "gene_id (column 1)") : (m, "gene_name (column 2)")
    else
        (findmatches(1), "gene name")
    end
    isempty(row_matches) &&
        error(tenx ?
            "Gene '$gene' not found in $features_path (searched gene_name -- column 2 -- then gene_id -- column 1)" :
            "Gene '$gene' not found in $features_path")
    length(row_matches) > 1 &&
        error("Gene '$gene' matches $(length(row_matches)) rows in $features_path ($searched); " *
              "provide its unique gene id instead")

    return row_matches[1]

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

# ---------------------------------------------------------------------------- #
#                  Fast byte-level Matrix Market scanning                       #
# ---------------------------------------------------------------------------- #
#
# eachline() allocates one String per line -- for hundreds of millions of
# entries that is pure GC churn and dominates the scan time. The scanner
# below reads the data section in ~1 MiB chunks and parses row/col/value
# straight from the byte buffer, allocating only a couple of buffers per
# chunk, which brings a full-file scan down to roughly decompression speed.

# Parse the unsigned integer starting at buf[i]; returns (value, next index).
@inline function _parse_uint(buf::Vector{UInt8}, i::Int, iend::Int)
    v = 0
    @inbounds while i <= iend
        b = buf[i]
        (b < UInt8('0') || b > UInt8('9')) && break
        v = 10v + Int(b - UInt8('0'))
        i += 1
    end
    return v, i
end

# Parse the floating-point token spanning buf[i:iend] (trailing spaces/tabs/\r
# allowed). Fast exact path for plain decimals ("3", "1.0", "1.25e-3");
# falls back to Base.parse for anything unusual (long mantissas, extreme
# exponents) so results always match a standard parse. Returns (val, ok).
function _parse_mm_float(buf::Vector{UInt8}, i::Int, iend::Int)
    istart = i
    # trim trailing whitespace
    @inbounds while iend >= i && (buf[iend] == UInt8(' ') || buf[iend] == UInt8('\t') || buf[iend] == UInt8('\r'))
        iend -= 1
    end
    iend < i && return 0.0, false
    neg = false
    @inbounds begin
        if buf[i] == UInt8('-') || buf[i] == UInt8('+')
            neg = buf[i] == UInt8('-')
            i += 1
        end
        mant = 0.0
        ndig = 0
        while i <= iend && UInt8('0') <= buf[i] <= UInt8('9')
            mant = 10.0 * mant + (buf[i] - UInt8('0'))
            i += 1; ndig += 1
        end
        frac = 0
        if i <= iend && buf[i] == UInt8('.')
            i += 1
            while i <= iend && UInt8('0') <= buf[i] <= UInt8('9')
                mant = 10.0 * mant + (buf[i] - UInt8('0'))
                i += 1; frac += 1; ndig += 1
            end
        end
        expo = 0
        if i <= iend && (buf[i] == UInt8('e') || buf[i] == UInt8('E'))
            i += 1
            eneg = false
            if i <= iend && (buf[i] == UInt8('-') || buf[i] == UInt8('+'))
                eneg = buf[i] == UInt8('-')
                i += 1
            end
            ev = 0
            edig = 0
            while i <= iend && UInt8('0') <= buf[i] <= UInt8('9')
                ev = 10ev + Int(buf[i] - UInt8('0'))
                i += 1; edig += 1
            end
            edig == 0 && return 0.0, false
            expo = eneg ? -ev : ev
        end
    end
    (ndig == 0 || i <= iend) && return 0.0, false  # no digits, or trailing junk
    e = expo - frac
    if ndig <= 15 && -22 <= e <= 22
        # mant is exact (< 2^53) and 10^|e| is exact: result correct to 1 ulp
        val = e >= 0 ? mant * (10.0^e) : mant / (10.0^(-e))
        return neg ? -val : val, true
    end
    # rare/extreme token: delegate to Base for full correctness
    return parse(Float64, String(buf[istart:iend])), true
end

"""
    read_mtx_header(io) -> (n_rows, n_cols, nnz)

Consumes the `%` comment lines and the dimensions line of an open Matrix
Market coordinate stream, leaving `io` positioned at the first data entry.
"""
function read_mtx_header(io::IO)
    line = readline(io)
    while startswith(line, "%")
        eof(io) && error("Matrix Market stream has no data after its comment header")
        line = readline(io)
    end
    dims = split(line)
    length(dims) == 3 || error("Unexpected Matrix Market dimensions line: '$line'")
    n_rows, n_cols, nnz = parse.(Int, dims)
    return n_rows, n_cols, nnz
end

"""
    scan_mtx_data(f, io; only_row = 0, rows_only = false)

Streams the data section of a Matrix Market coordinate stream (header
already consumed via [`read_mtx_header`](@ref)) in byte chunks, calling
`f(row, col, val)` per entry without per-line allocations. With
`rows_only = true`, only the row id is parsed and `f(row, 0, 0.0)` is
called for every entry; with `only_row > 0`, `f` is called only for that
row's entries (other lines are skipped after parsing their row id).
"""
function scan_mtx_data(f::F, io::IO; only_row::Int = 0, rows_only::Bool = false) where F

    chunksize = 1 << 20
    leftover = UInt8[]

    process = function (buf::Vector{UInt8}, stop::Int)
        i = 1
        @inbounds while i <= stop
            nl = i
            while nl <= stop && buf[nl] != 0x0a
                nl += 1
            end
            lend = nl - 1
            lend >= i && buf[lend] == 0x0d && (lend -= 1)
            if lend >= i
                row, j = _parse_uint(buf, i, lend)
                if rows_only
                    f(row, 0, 0.0)
                elseif only_row == 0 || row == only_row
                    while j <= lend && (buf[j] == UInt8(' ') || buf[j] == UInt8('\t'))
                        j += 1
                    end
                    col, j = _parse_uint(buf, j, lend)
                    while j <= lend && (buf[j] == UInt8(' ') || buf[j] == UInt8('\t'))
                        j += 1
                    end
                    val, ok = _parse_mm_float(buf, j, lend)
                    ok || error("Malformed Matrix Market entry: '$(String(buf[i:min(lend, i + 80)]))'")
                    f(row, col, val)
                end
            end
            i = nl + 1
        end
        nothing
    end

    while !eof(io)
        chunk = read(io, chunksize)
        buf = isempty(leftover) ? chunk : vcat(leftover, chunk)
        lastnl = findlast(==(0x0a), buf)
        if lastnl === nothing
            leftover = buf
            continue
        end
        process(buf, lastnl)
        leftover = buf[lastnl + 1:end]
    end
    isempty(leftover) || process(push!(leftover, 0x0a), length(leftover) + 1)

    return nothing
end

# ---------------------------------------------------------------------------- #
#              Gene-major indexed sidecar (.dgx) for instant loads             #
# ---------------------------------------------------------------------------- #
#
# Matrix Market exports are cell-major and (usually) gzipped, so extracting
# one gene inherently costs a full decompress-and-scan of the file. The
# .dgx sidecar produced by `prepare_gene_expression` (or the
# `dynema-prepare-expr` CLI) reorders the nonzero entries by gene once, with
# an offset table, so any single gene loads with two small seeks+reads --
# milliseconds instead of minutes. Layout (native-endian):
#
#   bytes 0-7    magic "DYNEMAGX"
#   Int64        format version (1)
#   Int64 × 3    n_genes, n_cells, nnz
#   Int64 × (n_genes + 1)  entry-offset table (0-based, cumulative)
#   Int32 × nnz  cell indices, grouped by gene
#   Float32 × nnz  values, grouped by gene (exact for counts < 2^24)

const DGX_MAGIC = b"DYNEMAGX"
const DGX_VERSION = 1

"""
    dgx_sidecar_path(mtx) -> String

Path of the gene-major indexed sidecar for a given Matrix Market file:
strips `.gz` and `.mtx` and appends `.dgx` (e.g. `expr.mtx.gz` -> `expr.dgx`).
"""
function dgx_sidecar_path(mtx::AbstractString)
    base = endswith(mtx, ".gz") ? String(mtx[1:end - 3]) : String(mtx)
    base = endswith(base, ".mtx") ? base[1:end - 4] : base
    return base * ".dgx"
end

"""
    prepare_gene_expression(; mtx, out = nothing, verbose = true) -> String

One-time conversion of a (possibly gzipped) cell-major Matrix Market count
matrix into a gene-major indexed binary sidecar (`.dgx`, written next to
`mtx` unless `out` is given), after which `extract_gene_expression` -- and
therefore `dynema-map --expr-prefix` -- loads any single gene in
milliseconds instead of rescanning the whole file.

Runs two streaming passes over `mtx` (count entries per gene, then place
them), so it needs roughly 8 bytes of RAM per nonzero entry (~8 GB per
billion entries) plus the file's decompression, and takes about twice one
full scan. Values are stored as `Float32` (exact for counts below 2^24).
"""
function prepare_gene_expression(; mtx::AbstractString,
                                 out::Union{Nothing, AbstractString} = nothing,
                                 verbose::Bool = true)

    outpath = out === nothing ? dgx_sidecar_path(mtx) : String(out)

    # ------------------------- Pass 1: count per gene ------------------------- #

    printlnv("Pass 1/2: counting entries per gene in $mtx..."; verbose)
    t0 = time()
    io = open_maybe_gzip(mtx)
    local n_rows, n_cols, nnz, counts
    try
        n_rows, n_cols, nnz = read_mtx_header(io)
        counts = zeros(Int64, n_rows)
        scan_mtx_data(io; rows_only = true) do row, _, _
            @inbounds counts[row] += 1
        end
    finally
        close(io)
    end
    total = sum(counts)
    total == nnz ||
        error("$mtx declares $nnz entries but contains $total")
    printlnv("  $n_rows genes × $n_cols cells, $nnz entries ($(round(time() - t0, digits = 1))s)"; verbose)

    offsets = Vector{Int64}(undef, n_rows + 1)
    offsets[1] = 0
    @inbounds for g in 1:n_rows
        offsets[g + 1] = offsets[g] + counts[g]
    end

    # -------------------------- Pass 2: place entries ------------------------- #

    printlnv("Pass 2/2: reordering entries by gene (~$(round((nnz * 8) / 2^30, digits = 1)) GiB in RAM)..."; verbose)
    t0 = time()
    cols = Vector{Int32}(undef, nnz)
    vals = Vector{Float32}(undef, nnz)
    cursor = offsets[1:end - 1]
    noninteger = Ref(0)
    io = open_maybe_gzip(mtx)
    try
        read_mtx_header(io)
        scan_mtx_data(io) do row, col, val
            @inbounds begin
                idx = (cursor[row] += 1)
                cols[idx] = col
                vals[idx] = val
            end
            val == round(val) || (noninteger[] += 1)
        end
    finally
        close(io)
    end
    cursor == offsets[2:end] || error("Entry counts changed between passes; was $mtx modified?")
    # Single-cell expression matrices for Dynema are raw counts; non-integer
    # values almost always mean a normalized/log-transformed matrix was
    # exported by mistake, which would silently violate the Poisson model.
    noninteger[] > 0 &&
        @warn "$(noninteger[]) of $nnz entries in $mtx are not integer counts. Dynema models " *
              "raw counts (Poisson); make sure this matrix is not normalized/log-transformed."
    printlnv("  done ($(round(time() - t0, digits = 1))s)"; verbose)

    # --------------------------------- Write --------------------------------- #

    open(outpath, "w") do out_io
        write(out_io, DGX_MAGIC)
        write(out_io, Int64(DGX_VERSION))
        write(out_io, Int64(n_rows), Int64(n_cols), Int64(nnz))
        write(out_io, offsets)
        write(out_io, cols)
        write(out_io, vals)
    end
    printlnv("Wrote $outpath ($(round(filesize(outpath) / 2^30, digits = 2)) GiB)"; verbose)

    return outpath

end

"""
    read_dgx_gene(path, row, n_genes, n_cells) -> (cols::Vector{Int32}, vals::Vector{Float32})

Reads one gene's nonzero entries from a `.dgx` sidecar (see
[`prepare_gene_expression`](@ref)), validating the header against the
expected gene/cell counts so a stale sidecar (built from different
features/barcodes) fails loudly instead of silently misaligning.
"""
function read_dgx_gene(path::AbstractString, row::Int, n_genes::Int, n_cells::Int)

    open(path, "r") do io
        magic = read(io, length(DGX_MAGIC))
        magic == DGX_MAGIC || error("$path is not a Dynema .dgx sidecar (bad magic)")
        ver = read(io, Int64)
        ver == DGX_VERSION || error("$path has .dgx format version $ver; this Dynema reads version $DGX_VERSION")
        nr = read(io, Int64)
        nc = read(io, Int64)
        nnz = read(io, Int64)
        nr == n_genes ||
            error("$path was built for $nr genes but the features file has $n_genes -- rebuild the sidecar (dynema-prepare-expr)")
        nc == n_cells ||
            error("$path was built for $nc cells but the barcodes file has $n_cells -- rebuild the sidecar (dynema-prepare-expr)")
        1 <= row <= nr || error("Gene row $row out of range 1:$nr in $path")

        hdr = length(DGX_MAGIC) + 4 * sizeof(Int64)
        seek(io, hdr + (row - 1) * sizeof(Int64))
        o1 = read(io, Int64)
        o2 = read(io, Int64)
        k = Int(o2 - o1)

        colsstart = hdr + (nr + 1) * sizeof(Int64)
        valsstart = colsstart + nnz * sizeof(Int32)

        cols = Vector{Int32}(undef, k)
        seek(io, colsstart + o1 * sizeof(Int32))
        read!(io, cols)
        vals = Vector{Float32}(undef, k)
        seek(io, valsstart + o1 * sizeof(Float32))
        read!(io, vals)

        return cols, vals
    end

end

"""
    extract_gene_expression(; mtx, features, barcodes, gene, id_col="cell_id", verbose=true)

Core Matrix Market single-gene extraction: streams through `mtx` once,
extracting the expression of `gene` (matched against `features`, one entry
per matrix row) across every cell (one entry per matrix column, named from
`barcodes`), without ever loading the whole sparse matrix into memory (see
module-level comment above). Matrix entries not present in the sparse file
are zero.

`gene` is a gene name/symbol or a gene id. The `features` file may be a
single column of gene names (e.g. a Seurat export), in which case `gene`
must match that column, or the two/three-column 10x layout
(gene_id, gene_name[, modality]), in which case `gene` is searched against
gene_name (column 2) first and, failing that, against gene_id (column 1).
Ensembl id version suffixes are ignored; no match is an error.

Returns a `NamedTuple` with:
- `expr::DataFrame`: two columns (`id_col`, then `gene` as provided), shaped
  exactly like a plain `dynema-map --expr` TSV would be.
- `target_row`, `n_genes`: which row `gene` was found at, and how many genes/rows
  `features` has in total.
- `n_cells`: number of cells/columns (from `barcodes`).
- `n_found`: number of nonzero entries read for `gene`.

If a gene-major indexed sidecar (`.dgx`, built once with
[`prepare_gene_expression`](@ref) or the `dynema-prepare-expr` CLI) exists
next to `mtx`, the gene is loaded from it in milliseconds; otherwise the
whole file is scanned once (streaming, roughly decompression-bound).

Set `verbose=false` to suppress the built-in progress `println`s -- e.g. for
a caller (like the `dynema-map` CLI) that wants to render its own progress
messages from the returned counts instead.
"""
function extract_gene_expression(; mtx::AbstractString, features::AbstractString,
                                  barcodes::AbstractString, gene::AbstractString,
                                  id_col::AbstractString = "cell_id",
                                  verbose::Bool = true)

    printlnv("Reading gene list from $features..."; verbose)
    feats = read_feature_fields(features)
    n_genes = length(feats)
    target_row = match_feature_row(feats, gene; features_path = features)

    printlnv("Reading cell barcodes from $barcodes..."; verbose)
    cell_ids = read_labels(barcodes)
    n_cells = length(cell_ids)

    expr = zeros(Float64, n_cells)
    n_found = 0

    dgx = dgx_sidecar_path(mtx)
    if isfile(dgx)

        # Indexed sidecar present: two seeks + two small reads.
        printlnv("Loading gene '$gene' (row $target_row) from indexed sidecar $dgx..."; verbose)
        gcols, gvals = read_dgx_gene(dgx, target_row, n_genes, n_cells)
        @inbounds for l in eachindex(gcols)
            expr[gcols[l]] = Float64(gvals[l])
        end
        n_found = length(gcols)

    else

        # No sidecar: full streaming scan via the byte-level scanner (see
        # scan_mtx_data) -- roughly decompression-bound. Building a sidecar
        # once with prepare_gene_expression/dynema-prepare-expr turns this
        # into a milliseconds-scale indexed read for every later gene.
        printlnv("Scanning $mtx for gene '$gene' (row $target_row of $n_genes)..."; verbose)
        printlnv("  (tip: run dynema-prepare-expr once to index this matrix for instant gene loads)"; verbose)
        io = open_maybe_gzip(mtx)
        try
            n_rows, n_cols, _ = read_mtx_header(io)
            n_rows == n_genes ||
                error("$mtx has $n_rows rows but $features has $n_genes genes -- files don't match")
            n_cols == n_cells ||
                error("$mtx has $n_cols columns but $barcodes has $n_cells barcodes -- files don't match")
            scan_mtx_data(io; only_row = target_row) do _, col, val
                expr[col] = val
                n_found += 1
            end
        finally
            close(io)
        end

    end

    printlnv("Found $n_found nonzero entries for '$gene' across $n_cells cells"; verbose)

    df = DataFrame()
    df[!, id_col] = cell_ids
    df[!, gene] = expr

    return (; expr = df, target_row, n_genes, n_cells, n_found)

end
