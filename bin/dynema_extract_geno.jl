#!/usr/bin/env julia
# ---------------------------------------------------------------------------- #
#                             dynema_extract_geno.jl                           #
# ---------------------------------------------------------------------------- #
#
# Extracts a donor x variant dosage matrix for the cis-window around one
# gene's TSS from a bgzipped, tabix-indexed VCF (vcf.gz + .tbi/.csi), in
# exactly the format expected by `dynema_map.jl --geno`. (`dynema_map.jl` can
# also do this extraction on the fly via its own `--vcf` option, using the
# same underlying code in vcf_geno.jl -- use this script instead when you
# want the extracted matrix saved to a file, e.g. to reuse across runs or
# inspect it directly.)
#
# Assumes the VCF already carries per-sample genotype dosages (FORMAT field
# `DS`) and/or genotype probabilities (FORMAT field `GP`) -- as produced by
# standard imputation pipelines (Minimac, IMPUTE2/5, Beagle, etc.). Hard-call
# genotypes (`GT`) are not read; a variant with neither `DS` nor `GP` is
# skipped.
#
# The VCF must be indexed (`tabix -p vcf file.vcf.gz`), but `tabix` itself
# does not need to be installed system-wide: this script depends on
# `htslib_jll`, which downloads a self-contained htslib build (including
# `tabix`/`bgzip`) as a regular Julia package artifact, for every common
# platform. No conda/brew/apt step required.
#
# Usage:
#   ./dynema-extract-geno --vcf genotypes.vcf.gz --chr chr17 --tss 43044295 \
#       --window 250000 --out BRCA1_geno.tsv
#
# or, looking up the TSS from an annotation file:
#   ./dynema-extract-geno --vcf genotypes.vcf.gz --tss-file genes.tsv \
#       --gene BRCA1 --window 250000 --out BRCA1_geno.tsv
#
# The console transcript, plus the exact command run, is also saved to
# --log (default: --out with its extension swapped for .log).
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
# comment in dynema_map.jl -- this environment is shared between scripts.
!isfile(joinpath(CLI_DIR, "Manifest.toml")) &&
    println("First run detected: setting up the shared bin/ environment (this can take a few minutes)...")
Pkg.develop(Pkg.PackageSpec(path = joinpath(CLI_DIR, "..")); io = devnull)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames

include(joinpath(CLI_DIR, "cli_output.jl")) # section, bullet -- used by vcf_geno.jl's extract_geno_dataframe
include(joinpath(CLI_DIR, "vcf_geno.jl")) # resolve_region, run_tabix, variant_dosage, extract_geno_dataframe

# ---------------------------------------------------------------------------- #
#                              Argument parsing                                #
# ---------------------------------------------------------------------------- #

function parse_commandline()

    s = ArgParseSettings(
        prog = "dynema_extract_geno.jl",
        description = "Extract a donor x variant dosage matrix around a gene's TSS from a tabix-indexed VCF, ready for dynema_map.jl --geno.",
    )

    @add_arg_table! s begin
        "--vcf"
            help = "Path to a bgzipped, tabix-indexed VCF (.vcf.gz with a .vcf.gz.tbi or .vcf.gz.csi index)."
            arg_type = String
            required = true
        "--chr"
            help = "Chromosome of the gene's TSS (must match the VCF's chromosome naming, e.g. 'chr1' vs '1'). Ignored if --gene/--tss-file are given."
            default = nothing
        "--tss"
            help = "Transcription start site position (1-based). Ignored if --gene/--tss-file are given."
            arg_type = Int
        "--tss-file"
            help = "TSV/CSV with columns 'gene_id', 'chr', 'tss' to look up --gene's TSS from, instead of passing --chr/--tss directly."
            default = nothing
        "--gene"
            help = "Gene id to look up in --tss-file."
            default = nothing
        "--window"
            help = "Cis-window half-width in bp around the TSS (region queried is TSS +/- window)."
            arg_type = Int
            default = 1_000_000
        "--field"
            help = "Which FORMAT field to convert to dosage: 'auto' (prefer GP, fall back to DS per-variant), 'GP', or 'DS'. With GP/DS forced, a variant lacking that field is skipped."
            arg_type = String
            default = "auto"
            range_tester = x -> x in ("auto", "GP", "DS")
        "--samples"
            help = "Optional TSV/CSV with columns 'vcf_id','donor_id' mapping VCF sample names to single-cell donor ids. Default: VCF sample names are used as donor ids as-is."
            default = nothing
        "--donor-col"
            help = "Column name to use for the donor id column in the output (should match dynema_map.jl's --donor-col)."
            arg_type = String
            default = "donor_id"
        "--maf"
            help = "Minimum minor allele frequency (computed from dosage, among the output samples) to retain a variant."
            arg_type = Float64
            default = 0.0
        "--max-missing"
            help = "Maximum fraction of samples allowed to have a missing dosage/GP value before a variant is dropped. Retained missing values are mean-imputed."
            arg_type = Float64
            default = 0.1
        "--out"
            help = "Output path for the donor x variant dosage matrix (TSV)."
            arg_type = String
            required = true
        "--log"
            help = "Path to save a full transcript of this run (the exact command invoked, plus everything printed to the console). Default: --out with its extension replaced by .log."
            default = nothing
    end

    return parse_args(s)

end

# ---------------------------------------------------------------------------- #
#                                     Main                                     #
# ---------------------------------------------------------------------------- #

function run_extract(args)

    out_df = extract_geno_dataframe(
        vcf = args["vcf"],
        chr = args["chr"],
        tss = args["tss"],
        tss_file = args["tss-file"],
        gene = args["gene"],
        window = args["window"],
        field = args["field"],
        samples_file = args["samples"],
        donor_col = args["donor-col"],
        maf = args["maf"],
        max_missing = args["max-missing"],
    )

    CSV.write(args["out"], out_df; delim = '\t')
    section("Done")
    bullet("wrote $(ncol(out_df) - 1) variant(s) x $(nrow(out_df)) donor(s) to $(args["out"])")

end

function main()

    args = parse_commandline()

    log_path = args["log"] !== nothing ? args["log"] : splitext(args["out"])[1] * ".log"
    command  = "dynema_extract_geno.jl " * join(shell_quote.(ARGS), " ")

    section("Command")
    bullet(command)

    with_tee_log(log_path; command = command) do
        run_extract(args)
    end

end

main()
