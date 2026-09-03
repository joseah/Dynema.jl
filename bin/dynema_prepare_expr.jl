#!/usr/bin/env julia
# ---------------------------------------------------------------------------- #
#                            dynema_prepare_expr.jl                            #
# ---------------------------------------------------------------------------- #
#
# One-time conversion of a (possibly gzipped) cell-major Matrix Market count
# matrix into a gene-major indexed binary sidecar (.dgx, written next to the
# .mtx). Afterwards, `dynema-map --expr-prefix` (and
# `Dynema.extract_gene_expression`) loads any single gene in milliseconds by
# seeking straight to its entries, instead of decompressing and scanning the
# whole multi-GB file for every gene.
#
# Thin CLI wrapper around `Dynema.prepare_gene_expression` (a core library
# function, callable from any Julia session without this wrapper). Runs two
# streaming passes over the matrix, needing ~8 bytes of RAM per nonzero
# entry (~8 GB per billion entries).
#
# Usage:
#   ./dynema-prepare-expr --expr-prefix expr_0.05
# or:
#   ./dynema-prepare-expr --mtx expr_0.05.mtx.gz
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

# Always re-develop/instantiate (not just on first run): see the matching
# comment in dynema_map.jl -- this environment is shared between scripts, and
# first-run progress (registry update + dependency resolution, which can
# silently take several minutes) is shown instead of silenced.
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
using Dynema # prepare_gene_expression, resolve_mtx_triplet, dgx_sidecar_path

# ---------------------------------------------------------------------------- #
#                              Argument parsing                                #
# ---------------------------------------------------------------------------- #

function parse_commandline()

    s = ArgParseSettings(
        prog = "dynema_prepare_expr.jl",
        description = "Index a Matrix Market gene expression export for instant per-gene loads by dynema-map. Run once per matrix.",
    )

    @add_arg_table! s begin
        "--expr-prefix"
            help = "Filename prefix of the Matrix Market export: resolves to <prefix>.mtx[.gz] (plus .features/.barcodes, which are checked to exist but not read). Exactly one of --expr-prefix or --mtx must be given."
            default = nothing
        "--mtx"
            help = "Path to the Matrix Market file itself (.mtx or .mtx.gz). Exactly one of --expr-prefix or --mtx must be given."
            default = nothing
        "--out"
            help = "Output path for the sidecar. Default: alongside the matrix, with a .dgx extension (which is where dynema-map auto-detects it)."
            default = nothing
    end

    return parse_args(s)

end

function main()

    args = parse_commandline()

    (args["expr-prefix"] === nothing) == (args["mtx"] === nothing) &&
        error("Provide exactly one of --expr-prefix or --mtx")

    mtx = args["mtx"] !== nothing ? args["mtx"] : first(resolve_mtx_triplet(args["expr-prefix"]))
    isfile(mtx) || error("Matrix file not found: $mtx")

    t0 = time()
    out = prepare_gene_expression(mtx = mtx, out = args["out"])
    println("Done in $(round(time() - t0, digits = 1))s. dynema-map --expr-prefix will now use $out automatically.")

end

main()
