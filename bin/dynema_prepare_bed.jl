#!/usr/bin/env julia
# ---------------------------------------------------------------------------- #
#                             dynema_prepare_bed.jl                            #
# ---------------------------------------------------------------------------- #
#
# Builds the bed-like gene file(s) that drive `dynema-map` from a standard
# GTF annotation: one row per gene (chr, start, end, gene, strand), keeping
# only genes that are actually present in the expression data's features
# file (when --features is given), choosing for each gene an identifier that
# `dynema-map` will resolve unambiguously at map time, and optionally
# splitting the result into fixed-size chunk files -- one HPC job per chunk
# is the intended genome-wide pattern:
#
#   ./bin/dynema-prepare-bed --gtf gencode.gtf.gz --features expr.features.gz \
#       --chunk-size 100 --out beds/chunk
#   # then: one dynema-map job per beds/chunk_*.bed
#
# Identifier choice per gene record: the gene_name/symbol when it matches
# the features file uniquely (dynema-map searches features gene names
# first); the gene_id when the symbol is absent or ambiguous (duplicated in
# the features file, or annotated at multiple loci in the GTF). Genes not
# found in the features file at all are dropped (reported at the end) --
# they could never be mapped anyway. Gene order follows the GTF, so chunks
# are genomically contiguous.
#
# Run with --help for the full list of options.
#
# ---------------------------------------------------------------------------- #
#                        Self-bootstrap on first run                          #
# ---------------------------------------------------------------------------- #
# Shares bin/Project.toml with dynema_map.jl -- if that script has already
# been run once, this one needs no extra setup.

import Pkg

const CLI_DIR = @__DIR__
Pkg.activate(CLI_DIR)

first_run = !isfile(joinpath(CLI_DIR, "Manifest.toml"))
if first_run
    println("First run detected: setting up the shared bin/ environment " *
            "(installing dependencies and, if needed, updating Julia's package " *
            "registry -- this can take several minutes; please wait)...")
    flush(stdout)
end
Pkg.develop(Pkg.PackageSpec(path = joinpath(CLI_DIR, "..")); io = first_run ? stdout : devnull)
Pkg.instantiate()

using ArgParse
using Dynema # read_feature_fields, stripver, open_maybe_gzip (qualified)

# ---------------------------------------------------------------------------- #
#                              Argument parsing                                #
# ---------------------------------------------------------------------------- #

function parse_commandline()

    s = ArgParseSettings(
        prog = "dynema_prepare_bed.jl",
        description = "Build dynema-map's bed-like gene file(s) from a GTF annotation, matched against an expression features file and optionally split into chunks for HPC batching.",
    )

    @add_arg_table! s begin
        "--gtf"
            help = "GTF gene annotation (.gtf or .gtf.gz), e.g. GENCODE or the Cell Ranger reference's genes.gtf."
            arg_type = String
            required = true
        "--features"
            help = "The expression data's features file (single-column gene names, or the 10x gene_id/gene_name[/modality] layout; plain or gzipped). Only genes present here are written, and identifiers are chosen to resolve unambiguously against it at map time. Strongly recommended; without it, every GTF gene record is written using its gene_name."
            default = nothing
        "--feature-type"
            help = "GTF feature type (column 3) to take gene records from. The default, 'gene', skips transcript/exon/CDS records."
            arg_type = String
            default = "gene"
        "--chunk-size"
            help = "Split the output into files of at most this many genes each, named <out>_001.bed, <out>_002.bed, ... (one dynema-map job per chunk is the intended genome-wide pattern). 0 writes everything to a single <out>.bed."
            arg_type = Int
            default = 0
        "--out"
            help = "Output prefix (directories are created as needed)."
            arg_type = String
            default = "genes"
    end

    return parse_args(s)

end

# ---------------------------------------------------------------------------- #
#                                   Helpers                                    #
# ---------------------------------------------------------------------------- #

# Extract the value of `key "value"` from a GTF attributes field.
function gtf_attr(attrs::AbstractString, key::AbstractString)
    r = findfirst(key * " \"", attrs)
    r === nothing && return nothing
    vstart = last(r) + 1
    vend = findnext('"', attrs, vstart)
    vend === nothing && return nothing
    return String(attrs[vstart:vend - 1])
end

# ---------------------------------------------------------------------------- #
#                                     Main                                     #
# ---------------------------------------------------------------------------- #

function main()

    args = parse_commandline()

    # ------------------------ Features identifier sets ------------------------- #

    # name_counts: how often each (version-stripped) gene name appears in the
    # features file -- a name occurring more than once cannot be resolved at
    # map time, so those genes fall back to their gene_id. ids: the gene_id
    # column (10x layout only).
    name_counts = Dict{String, Int}()
    ids = Set{String}()
    tenx = false
    if args["features"] !== nothing
        feats = Dynema.read_feature_fields(args["features"])
        tenx = any(row -> length(row) >= 2, feats)
        for row in feats
            nm = Dynema.stripver(row[min(2, length(row))])  # col 2 (10x) or col 1 (single-column)
            name_counts[nm] = get(name_counts, nm, 0) + 1
            tenx && length(row) >= 1 && push!(ids, Dynema.stripver(row[1]))
        end
        println("Features file: $(length(feats)) gene(s) ($(tenx ? "10x gene_id/gene_name layout" : "single-column gene names"))")
    end
    have_features = args["features"] !== nothing

    # ------------------------------- Scan the GTF ------------------------------ #

    records = NamedTuple[]  # (chr, start, stop, name, id, strand)
    gtf_name_counts = Dict{String, Int}()  # duplicated symbols within the GTF itself

    io = Dynema.open_maybe_gzip(args["gtf"])
    n_records = 0
    try
        for line in eachline(io)
            (isempty(line) || startswith(line, "#")) && continue
            f = split(line, '\t')
            length(f) >= 9 || continue
            f[3] == args["feature-type"] || continue
            n_records += 1
            name = gtf_attr(f[9], "gene_name")
            id = gtf_attr(f[9], "gene_id")
            (name === nothing && id === nothing) && continue
            push!(records, (chr = String(f[1]), start = String(f[2]), stop = String(f[3]),
                            name = name === nothing ? nothing : Dynema.stripver(name),
                            id = id === nothing ? nothing : Dynema.stripver(id),
                            strand = String(f[7])))
            name !== nothing && (gtf_name_counts[Dynema.stripver(name)] = get(gtf_name_counts, Dynema.stripver(name), 0) + 1)
        end
    finally
        close(io)
    end
    println("GTF: $n_records '$(args["feature-type"])' record(s)")

    # ------------------------ Choose one identifier per gene -------------------- #

    rows = Vector{Tuple{String, String, String, String, String}}()  # chr start stop gene strand
    n_by_id = 0
    n_not_in_features = 0
    n_unresolvable = 0
    seen = Set{String}()

    for rec in records
        # A gene name is usable when it is unique within the GTF and (if a
        # features file was given) present in it exactly once -- the
        # conditions under which dynema-map's name-first search resolves it
        # unambiguously. Otherwise fall back to the gene id (resolvable only
        # against a 10x-layout features file, or when no features file was
        # given).
        name_in_features = rec.name !== nothing ? get(name_counts, rec.name, 0) : 0
        name_ok = rec.name !== nothing &&
                  get(gtf_name_counts, rec.name, 0) == 1 &&
                  (!have_features || name_in_features == 1)
        id_ok = rec.id !== nothing && (!have_features || (tenx && rec.id in ids))

        identifier = if name_ok
            rec.name
        elseif id_ok
            n_by_id += 1
            rec.id
        elseif have_features && name_in_features == 0
            n_not_in_features += 1
            continue
        else
            n_unresolvable += 1
            continue
        end

        # A repeated identifier would collide on output filenames and be
        # rejected by dynema-map's bed reader; keep the first occurrence.
        identifier in seen && (n_unresolvable += 1; continue)
        push!(seen, identifier)

        push!(rows, (rec.chr, rec.start, rec.stop, identifier, rec.strand))
    end

    isempty(rows) && error("No usable gene records found -- check that --gtf and --features belong to the same annotation/naming scheme")

    # ---------------------------------- Write ----------------------------------- #

    outdir = dirname(args["out"])
    isempty(outdir) || mkpath(outdir)

    writebed(path, rs) = open(path, "w") do out_io
        for r in rs
            println(out_io, join(r, '\t'))
        end
    end

    nchunk = args["chunk-size"]
    files = String[]
    if nchunk <= 0
        path = args["out"] * ".bed"
        writebed(path, rows)
        push!(files, path)
    else
        n_files = cld(length(rows), nchunk)
        pad = max(3, length(string(n_files)))
        for (ci, lo) in enumerate(1:nchunk:length(rows))
            hi = min(lo + nchunk - 1, length(rows))
            path = args["out"] * "_" * lpad(string(ci), pad, '0') * ".bed"
            writebed(path, rows[lo:hi])
            push!(files, path)
        end
    end

    println("Wrote $(length(rows)) gene(s) to $(length(files)) file(s): $(first(files))$(length(files) > 1 ? " ... $(last(files))" : "")")
    n_by_id > 0 && println("  $n_by_id gene(s) written by gene_id (symbol absent or ambiguous)")
    n_not_in_features > 0 && println("  $n_not_in_features GTF gene(s) skipped: not present in --features")
    n_unresolvable > 0 && println("  $n_unresolvable gene(s) skipped: no unambiguous identifier")

end

main()
