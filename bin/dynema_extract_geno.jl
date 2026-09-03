#!/usr/bin/env julia
# ---------------------------------------------------------------------------- #
#                             dynema_extract_geno.jl                           #
# ---------------------------------------------------------------------------- #
#
# Extracts a donor x variant dosage matrix for the cis-window around one
# gene's TSS from a bgzipped, tabix-indexed VCF (vcf.gz + .tbi/.csi), in
# exactly the format expected by `dynema_map.jl --geno`. Thin CLI wrapper
# around `Dynema.extract_geno_dataframe` -- a core library function, also
# used on the fly by `dynema_map.jl --vcf`, and callable directly from any
# Julia session/script (`using Dynema`) without either CLI wrapper. Use
# this script instead of `dynema_map.jl --vcf` when you want the extracted
# matrix saved to a file, e.g. to reuse across runs or inspect it directly.
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
#   ./dynema-extract-geno --vcf genotypes.vcf.gz --bed BRCA1.bed \
#       --window 250000 --out BRCA1_geno.tsv
#
# The gene's TSS is derived from the bed-like annotation the way FastQTL
# does it: the gene's start position on the + strand, its end position on
# the - strand.
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
using CSV
using DataFrames
using Dynema # extract_geno_dataframe

include(joinpath(CLI_DIR, "cli_output.jl")) # section, bullet, with_tee_log, shell_quote

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
        "--bed"
            help = "Single-gene bed-like file that specifies which gene's cis-window to extract: a plain-text (optionally gzipped) table with exactly one data row and columns chr, start, end, gene, strand (standard 6-column BED with a score column also works; header/# lines are skipped). The TSS is derived FastQTL-style: start on the + strand, end on the - strand. Chromosome naming must match the VCF's (e.g. 'chr1' vs '1')."
            arg_type = String
            required = true
        "--window"
            help = "Cis-window half-width in bp around the TSS (region queried is TSS +/- window)."
            arg_type = Int
            default = 500_000
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

    section("Extracting genotypes from VCF")
    bullet("file: $(args["vcf"])")

    r = extract_geno_dataframe(
        vcf = args["vcf"],
        bed = args["bed"],
        window = args["window"],
        field = args["field"],
        samples_file = args["samples"],
        donor_col = args["donor-col"],
        maf = args["maf"],
        max_missing = args["max-missing"],
        verbose = false,
    )
    bullet("cis-window: $(r.chr):$(r.start_pos)-$(r.end_pos) (TSS $(r.tss) +/- $(args["window"]) bp)")
    bullet("samples: $(r.n_samples_vcf) in VCF, $(r.n_samples_matched) retained after sample matching")
    bullet("variants in region: $(r.n_variants_total)")
    bullet("skipped (multiallelic): $(r.n_multiallelic)", indent = 2)
    bullet("skipped (missing GP/DS field): $(r.n_no_field)", indent = 2)
    bullet("skipped (missingness > $(args["max-missing"])): $(r.n_high_missing)", indent = 2)
    bullet("retained after MAF >= $(args["maf"]) filter: $(r.n_retained)", indent = 2)

    CSV.write(args["out"], r.geno; delim = '\t')
    section("Done")
    bullet("wrote $(r.n_retained) variant(s) x $(nrow(r.geno)) donor(s) to $(args["out"])")

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
