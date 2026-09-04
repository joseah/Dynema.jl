#!/usr/bin/env julia
# ---------------------------------------------------------------------------- #
#                                 dynema_map.jl                                #
# ---------------------------------------------------------------------------- #
#
# Command-line wrapper around `Dynema.map_locus` for users who don't want to
# write any Julia code. Reads single-cell gene expression either from a plain
# --expr TSV/CSV or extracted on the fly for one gene from an --expr-prefix
# Matrix Market triplet (via Dynema.extract_gene_expression), reads
# single-cell metadata from a plain TSV/CSV, gets donor-level genotype
# dosages either from a pre-extracted --geno matrix or extracted on the fly
# for the tested gene's cis-window from a --vcf (via
# Dynema.extract_geno_dataframe), builds the eQTL model formula from a
# handful of flags, runs Dynema for one gene against one or more genetic
# variants, and writes a summary statistics table. The full console
# transcript, plus the exact command run, is also saved to --log (default:
# --out with its extension swapped for .log) -- see cli_output.jl.
#
# The VCF and Matrix Market extraction this script uses are core Dynema
# library functions, not CLI-only code -- `Dynema.extract_geno_dataframe`
# and `Dynema.extract_gene_expression`/`Dynema.resolve_mtx_triplet` can be
# called directly from any Julia session/script (`using Dynema`), without
# this CLI wrapper at all.
#
# Usage (pre-extracted genotype matrix):
#   julia --project=<this-dir> dynema_map.jl --geno geno.tsv [options]
#
# or extracting genotypes on the fly from a tabix-indexed VCF:
#   julia --project=<this-dir> dynema_map.jl --vcf genotypes.vcf.gz \
#       --bed BRCA1.bed --window 250000 [options]
#
# or, via the bundled launcher (auto-activates this environment):
#   ./dynema-map [options]
#
# Run with --help for the full list of options and see README.md /
# docs/src/tutorials for a description of the expected input file layout.
#
# ---------------------------------------------------------------------------- #
#                        Self-bootstrap on first run                          #
# ---------------------------------------------------------------------------- #

import Pkg

const CLI_DIR = @__DIR__
Pkg.activate(CLI_DIR)

# Always re-develop/instantiate (not just on first run): this environment is
# shared with dynema_extract_geno.jl, so a dependency added for one script
# must not be left unresolved when only the other has been run before. Both
# calls are fast no-ops once the manifest already satisfies Project.toml --
# but not on a genuine first run: Pkg.develop's own dependency resolution
# (and, on a brand-new Julia/juliaup install, updating the General registry
# before it can resolve anything) can silently take several minutes. That
# work is normally silenced below (`io = devnull`) so a repeat, already-set-up
# run doesn't print noisy "Resolving package versions..." on every
# invocation -- but on first run we deliberately let it print to stdout
# instead, so the user sees real progress rather than nothing.
first_run = !isfile(joinpath(CLI_DIR, "Manifest.toml"))
if first_run
    println("First run detected: setting up the shared bin/ environment " *
            "(installing dependencies and, if needed, updating Julia's package " *
            "registry -- this can take several minutes; please wait)...")
    flush(stdout)
end
Pkg.develop(Pkg.PackageSpec(path = joinpath(CLI_DIR, "..")); io = first_run ? stdout : devnull)
Pkg.instantiate()

# --workers N must take effect *before* `using Dynema` below, so that the
# package is automatically loaded on the new worker processes too. It is
# therefore scanned from the raw ARGS here; ArgParse still declares it for
# --help/validation. Equivalent to launching with `julia -p N`.
using Distributed
using LinearAlgebra: BLAS
let i = findfirst(==("--workers"), ARGS)
    if i !== nothing && i < length(ARGS)
        nw = tryparse(Int, ARGS[i + 1])
        if nw !== nothing && nw > 0 && nworkers() < nw
            println("Starting $nw worker process(es)...")
            flush(stdout)
            addprocs(nw - (nworkers() > 1 ? nworkers() : 0); exeflags = "--project=$CLI_DIR")
            # Divide the machine's BLAS threads across the workers: left at
            # its default, every worker's BLAS assumes it owns all cores, and
            # nw competing thread pools oversubscribe the CPU during the
            # GEMM/IRLS steps.
            blas_threads = max(1, Sys.CPU_THREADS ÷ nworkers())
            @everywhere workers() begin
                using LinearAlgebra: BLAS
                BLAS.set_num_threads($blas_threads)
            end
        end
    end
end

# --boot needs the optional WildBootTests package (Dynema's own direct score
# test covers everything else). On first --boot use, install it into this
# shared bin/ environment, then load it so Dynema's bootstrap extension
# activates. Scanned from the raw ARGS (like --workers) since ArgParse runs
# later; placed after the --workers block so a `using` here also reaches any
# freshly started worker processes.
if "--boot" in ARGS
    if !haskey(Pkg.project().dependencies, "WildBootTests")
        println("--boot uses the optional WildBootTests package; installing it into the " *
                "bin/ environment (one-time, may take a minute)...")
        flush(stdout)
        Pkg.add("WildBootTests"; io = stdout)
    end
    @eval using WildBootTests
end

using ArgParse
using CSV
using DataFrames
using StatsModels
using Dynema # also brings in extract_geno_dataframe, extract_gene_expression, resolve_mtx_triplet

include(joinpath(CLI_DIR, "cli_output.jl")) # section, bullet, elapsed_str, with_tee_log, shell_quote

# ---------------------------------------------------------------------------- #
#                              Argument parsing                                #
# ---------------------------------------------------------------------------- #

function parse_commandline()

    s = ArgParseSettings(
        prog = "dynema_map.jl",
        description = "Map single-cell eQTL effects for one gene with Dynema, from plain text input files. No Julia code required.",
    )

    @add_arg_table! s begin
        "--expr"
            help = "TSV/CSV file with single-cell gene expression counts. Must contain a cell-id column and one column per gene (or a single expression column if no --bed is given). Exactly one of --expr or --expr-prefix must be given."
            arg_type = String
        "--expr-prefix"
            help = "Filename prefix of a Matrix Market gene expression export (genes as rows, cells as columns): resolves to <prefix>.mtx, <prefix>.features, <prefix>.barcodes, each tried gzipped first then plain (e.g. 'expr_0.05' -> expr_0.05.mtx.gz/.features.gz/.barcodes.gz). If a <prefix>.dgx index built by dynema-prepare-expr exists, the gene loads from it in milliseconds; otherwise the matrix is streamed once for the tested gene, never loaded in full. Requires --bed (whose gene column names the gene to extract). Exactly one of --expr or --expr-prefix must be given."
            default = nothing
        "--geno"
            help = "TSV/CSV file with donor-level genotype dosages: one donor-id column plus one column per variant. Exactly one of --geno or --vcf must be given."
            arg_type = String
        "--vcf"
            help = "Bgzipped, tabix-indexed VCF (.vcf.gz) to extract the gene's cis-window genotypes from on the fly, instead of a pre-extracted --geno file. Requires --bed to locate the cis-window, plus --window/--field/--samples/--maf/--max-missing below. Exactly one of --geno or --vcf must be given."
            default = nothing
        "--bed"
            help = "Bed-like file that specifies WHAT to map and WHERE: a plain-text (optionally gzipped) table with columns chr, start, end, gene, strand (standard 6-column BED with a score column also works; header/# lines are skipped). Each data row is one gene to map, so a multi-row file defines a batch mapped sequentially in this run -- split a genome-wide bed into chunks and submit one dynema-map job per chunk to parallelize on a cluster. The gene column -- a name/symbol or a gene id -- names each gene: with a single-column features file (e.g. a Seurat export) it must match that column; with a 10x features file (gene_id, gene_name, ...) it is searched against gene_name first, then gene_id (Ensembl id version suffixes ignored). The TSS is derived FastQTL-style: start on the + strand, end on the - strand. Chromosome naming must match the VCF's."
            default = nothing
        "--window"
            help = "(--vcf only) Cis-window half-width in bp around the TSS."
            arg_type = Int
            default = 500_000
        "--field"
            help = "(--vcf only) Which FORMAT field to convert to dosage: 'auto' (prefer GP, fall back to DS per-variant), 'GP', or 'DS'."
            arg_type = String
            default = "auto"
            range_tester = x -> x in ("auto", "GP", "DS")
        "--samples"
            help = "(--vcf only) Optional TSV/CSV with columns 'vcf_id','donor_id' mapping VCF sample names to single-cell donor ids."
            default = nothing
        "--maf"
            help = "(--vcf only) Minimum minor allele frequency (among the matched donors) to retain a variant."
            arg_type = Float64
            default = 0.0
        "--max-missing"
            help = "(--vcf only) Maximum fraction of donors allowed to have a missing dosage/GP value before a variant is dropped; retained missing values are mean-imputed."
            arg_type = Float64
            default = 0.1
        "--meta"
            help = "TSV/CSV file with single-cell metadata: cell id, donor id, cell-state contexts, and any covariates -- one row per cell."
            arg_type = String
            required = true
        "--cell-id-col"
            help = "Column name for the cell identifier, shared across --expr and --meta."
            arg_type = String
            default = "cell_id"
        "--donor-col"
            help = "Column name for the donor identifier in --meta and --geno."
            arg_type = String
            default = "donor_id"
        "--variants"
            help = "Comma-separated list of variant/column names to test, from --geno or the --vcf extraction. Default: test every variant found."
            default = nothing
        "--covariates"
            help = "Comma-separated list of column names in --meta to include as additive covariates (e.g. age,sex,nUMI,percent_mito,gPC1,gPC2)."
            arg_type = String
            default = ""
        "--effect"
            help = "Effect(s) to test, comma-separated: any of 'main' (context-independent), 'interaction' (single- or multi-context), and 'total' (main + interaction jointly), e.g. 'main,interaction,total'. Always required -- no default or shortcut -- so what is being tested is always explicit. Multiple effects share each gene's data extraction and write separate output files, <out>_<effect>_<gene>.tsv. In a multi-effect run, 'main' uses the classic model without G x context terms; a standalone --effect main can still adjust for (untested) interactions via --interaction-with."
            arg_type = String
            required = true
            range_tester = x -> !isempty(x) && all(t -> strip(t) in ("main", "interaction", "total"), split(x, ","))
        "--interaction-with"
            help = "Comma-separated cell-state/context column names in --meta (e.g. C1,C2,C3) to include as G-by-context interaction terms; their main effects are added to the model automatically. Required (and must be non-empty) when --effect includes interaction/total -- never defaulted or guessed. Optional with a standalone --effect main: if given, the interaction terms are added to the formula (adjusted for) but not tested. Contexts you only want to adjust for (no interaction) belong in --covariates."
            default = nothing
        "--boot"
            help = "Also compute bootstrap p-values via adaptive score bootstrapping (recommended for small or imbalanced cohorts). Adds two columns: p_boot, the empirical bootstrap p-value (floored at ~2/B), and p_boot_approx, a FastQTL-style beta approximation fitted to the bootstrap distribution that extrapolates smoothly below that floor. Uses the optional WildBootTests package, installed automatically into the bin/ environment on first use."
            action = :store_true
        "--B"
            help = "Comma-separated adaptive bootstrap iteration schedule. Only used with --boot."
            arg_type = String
            default = "200,200,1600,2000,16000,20000"
        "--ptype"
            help = "Bootstrap p-value type: equaltail or symmetric. Only used with --boot."
            arg_type = String
            default = "equaltail"
            range_tester = x -> x in ("equaltail", "symmetric")
        "--betas"
            help = "Which variants to fit the unrestricted model for, attaching per-variant effect-size estimates as extra output columns (the score test itself never needs them): 'lead' (default) fits only the lead variant(s) -- smallest analytical p-value, including all exactly tied variants -- leaving the estimate columns empty for the rest; 'all' fits every variant (substantially slower); 'none' skips estimates entirely (fastest)."
            arg_type = String
            default = "lead"
            range_tester = x -> x in ("all", "none", "lead")
        "--parallel"
            help = "Distribute variants across worker processes with Distributed.jl. Use with --workers N, or start workers yourself (e.g. `julia -p 4 --project=bin bin/dynema_map.jl ...`)."
            action = :store_true
        "--workers"
            help = "Number of local worker processes to start for parallel mapping. Implies --parallel."
            arg_type = Int
            default = 0
        "--positions"
            help = "Optional TSV/CSV with two columns (variant, position) giving genomic positions to attach to the output. Only needed with --geno: with --vcf, positions are taken from the VCF automatically (this file overrides them if given)."
            default = nothing
        "--skip-existing"
            help = "Skip any gene x effect whose output file already exists -- lets a killed or partially completed batch job be resubmitted without recomputing finished genes. Skipped entries appear in the batch summary with status 'skipped' (lead statistics re-read from the existing file when possible)."
            action = :store_true
        "--check"
            help = "Validate all inputs and exit without mapping: files exist and parse, metadata columns are present, every --bed gene is found in the expression data, the annotation's chromosomes exist in the VCF index, and metadata donors are covered by the VCF samples. Run this before submitting large batches."
            action = :store_true
        "--out"
            help = "Output prefix. Each gene writes its own summary statistics table, named <out>_<gene>.tsv (or <out><gene>.tsv if the prefix ends with '/'; a trailing .tsv on the prefix is stripped; directories are created). Useful to distinguish analyses of the same genes, e.g. --out main vs --out interaction. Default: <gene>.tsv in the current directory."
            arg_type = String
            default = ""
        "--log"
            help = "Path to save a full transcript of this run (the exact command invoked, plus everything printed to the console -- progress, warnings, and the results tables, for every gene in the batch). Default: <out prefix>.log, or dynema_map.log if no --out is given."
            default = nothing
    end

    return parse_args(s)

end

# ---------------------------------------------------------------------------- #
#                                   Helpers                                    #
# ---------------------------------------------------------------------------- #

splitcsv(x::Nothing) = String[]
splitcsv(x::AbstractString) = isempty(strip(x)) ? String[] : String.(strip.(split(x, ",")))

readtable(path::AbstractString) = CSV.read(path, DataFrame)

"""
    resolve_interactions(effect, interaction_with_arg) -> Vector{String}

Resolves the context columns that should appear as `G & context` terms in
the model formula (passed to `build_formula`, and used to validate --meta's
columns before that):

- If `--effect` is `interaction`/`total` (which always test at least one
  `G & context` term), `--interaction-with` must be given explicitly and
  non-empty -- errors immediately otherwise, rather than guessing.
- If `--effect` is `main`, `--interaction-with` is optional: given
  explicitly, its interaction terms are still added to the formula (e.g. so
  `--effect main --interaction-with CV1` can adjust for a G x CV1
  interaction without testing it); omitted, no interaction terms are added.

Every interaction context automatically gets a main-effect term in the
model (see `build_formula`); contexts to adjust for *without* an
interaction belong in `--covariates`.
"""
function resolve_interactions(effect::String,
                               interaction_with_arg::Union{Nothing,Vector{String}})

    if effect in ("interaction", "total")
        (interaction_with_arg === nothing || isempty(interaction_with_arg)) &&
            error("--effect $effect requires --interaction-with (a non-empty, comma-separated " *
                  "list of context columns to test G-by-context interactions for)")
    end

    # Duplicates within --interaction-with would create identical G & context
    # terms (and a rank-deficient joint test). Error rather than silently
    # deduplicate: the user must be explicit about which interactions are
    # being tested.
    if interaction_with_arg !== nothing
        dupes = unique([c for c in interaction_with_arg if count(==(c), interaction_with_arg) > 1])
        isempty(dupes) ||
            error("--interaction-with lists context(s) more than once ($(join(dupes, ", "))); " *
                  "each context to test must appear exactly once")
    end

    return interaction_with_arg === nothing ? String[] : interaction_with_arg

end

"""
    build_formula(covariates, effect, interaction_with)

Builds a `FormulaTerm` and the `termtest` argument for `map_locus`
programmatically (equivalent to writing e.g. `@formula(0 ~ 1 + G + C1 + G & C1)`
by hand), from plain lists of column names supplied on the command line.
Each context in `interaction_with` (see `resolve_interactions`) contributes
both its main-effect term and a `G & context` term (an interaction without
its main effect would be an unusual, generally statistically improper
model); the interaction terms are only added to `termtest` (i.e. actually
tested, alongside/instead of `G`) when `effect` is `interaction`/`total`.
"""
function build_formula(covariates::Vector{String},
                        effect::String, interaction_with::Vector{String})

    rhs_terms = Any[term(1), term(:G)]
    for c in unique(covariates)
        push!(rhs_terms, term(Symbol(c)))
    end
    # A context listed in both --covariates and --interaction-with (easy to
    # do by mistake) is harmless: its main effect enters the model once.
    for c in interaction_with
        c in covariates || push!(rhs_terms, term(Symbol(c)))
    end
    for c in interaction_with
        push!(rhs_terms, term(:G) & term(Symbol(c)))
    end

    termtest = String[]

    if effect in ("main", "total")
        push!(termtest, "G")
    end

    if effect in ("interaction", "total")
        isempty(interaction_with) &&
            error("--effect $effect requires a non-empty --interaction-with (this should have been caught earlier by resolve_interactions)")
        for c in interaction_with
            push!(termtest, "G & $(c)")
        end
    end

    rhs = reduce(+, rhs_terms)
    f = FormulaTerm(term(0), rhs)

    return f, (length(termtest) == 1 ? termtest[1] : termtest)

end

# ---------------------------------------------------------------------------- #
#                                     Main                                     #
# ---------------------------------------------------------------------------- #

"""
    run_map(args; term_size = displaysize(stdout))

The actual body of `dynema_map.jl`, factored out from `main()` so it can be
run inside `with_tee_log` (see cli_output.jl) -- everything it prints via
`section`/`bullet`/`@warn`/`show` etc. is thereby captured to `--log` as
well as the console.

`term_size` should be the real terminal's `displaysize(stdout)`, captured in
`main()` *before* `with_tee_log` redirects stdout to a pipe -- a pipe isn't
a TTY and can't report a terminal size, so without this the results table
would fall back to wrapping at Julia's narrow non-interactive default width
instead of what actually fits the user's terminal.
"""
function run_map(args; term_size = displaysize(stdout))

    (args["geno"] === nothing) == (args["vcf"] === nothing) &&
        error("Provide exactly one of --geno (pre-extracted matrix) or --vcf (extract on the fly)")
    (args["expr"] === nothing) == (args["expr-prefix"] === nothing) &&
        error("Provide exactly one of --expr (TSV/CSV) or --expr-prefix (Matrix Market)")

    covariates           = splitcsv(args["covariates"])
    interaction_with_arg = args["interaction-with"] === nothing ? nothing : splitcsv(args["interaction-with"])

    # Effects to test: a comma-separated subset of main/interaction/total.
    # Each effect gets its own model; in a multi-effect run, 'main' is the
    # classic main-effect model *without* G x context terms (matching a
    # standalone --effect main run), while interaction/total include them.
    effects = unique(String.(strip.(split(args["effect"], ","))))

    # Fails immediately (before any file loading below) if interaction/total
    # was requested without a non-empty --interaction-with -- see
    # resolve_interactions's docstring.
    effect_iw = Dict(e => resolve_interactions(e,
                        (e == "main" && length(effects) > 1) ? nothing : interaction_with_arg)
                     for e in effects)
    all_iw = unique(reduce(vcat, values(effect_iw); init = String[]))
    variant_filter       = splitcsv(args["variants"])
    B                = parse.(Int, splitcsv(args["B"]))
    ptype            = Symbol(args["ptype"])
    id_col           = args["cell-id-col"]
    donor_col        = args["donor-col"]

    # ------------------------------- Gene batch ------------------------------- #

    # The --bed file specifies both WHAT to map (its gene column drives the
    # expression lookup) and WHERE (positions/strand give each TSS for --vcf
    # extraction). Each of its rows is one gene: a multi-row file is a batch,
    # mapped sequentially within this run.
    args["vcf"] !== nothing && args["bed"] === nothing &&
        error("--vcf requires --bed (bed-like gene file) to locate each gene's cis-window")
    args["expr-prefix"] !== nothing && args["bed"] === nothing &&
        error("--expr-prefix requires --bed (whose gene column names the gene(s) to extract; there's no header row to infer a single gene from)")
    genes = args["bed"] === nothing ? nothing : Dynema.read_gene_bed(args["bed"])
    genes !== nothing && length(genes) > 1 && args["geno"] !== nothing &&
        error("--geno provides a single pre-extracted cis-window, so it only combines with a " *
              "single-gene --bed ($(length(genes)) genes given); use --vcf for multi-gene batches")

    if genes !== nothing
        section("Gene batch")
        bullet("file: $(args["bed"])")
        shown = join([g.gene for g in first(genes, 8)], ", ") * (length(genes) > 8 ? ", ..." : "")
        bullet("$(length(genes)) gene(s): $shown")
    end

    # Validate --meta's columns (cell id, donor id, contexts, covariates) up
    # front, before the potentially slow genotype/expression reads below --
    # so a typo'd column name fails fast instead of after a multi-minute
    # VCF/Matrix Market extraction.
    section("Reading metadata")
    bullet("file: $(args["meta"])")
    meta = readtable(args["meta"])
    id_col in names(meta) || error("--meta is missing cell-id column '$id_col'")
    donor_col in names(meta) || error("--meta is missing donor-id column '$donor_col'")
    for c in vcat(covariates, all_iw)
        c in names(meta) || error("Column '$c' not found in --meta")
    end
    bullet("cell-id column: $id_col, donor-id column: $donor_col")
    isempty(covariates) || bullet("covariates: $(join(covariates, ", "))")
    isempty(all_iw) || bullet("interaction context(s): $(join(all_iw, ", "))")
    bullet("found $(nrow(meta)) cell(s)")
    bullet("found $(length(unique(meta[:, donor_col]))) donor(s)")

    # --------------------- Batch-shared inputs, read once --------------------- #

    expr_table = nothing
    mtx = features = barcodes = nothing
    if args["expr"] !== nothing
        section("Reading gene expression")
        bullet("file: $(args["expr"])")
        expr_t0 = time()
        expr_table = readtable(args["expr"])
        id_col in names(expr_table) || error("Expression table is missing cell-id column '$id_col'")
        bullet("done in $(elapsed_str(expr_t0))")
    else
        mtx, features, barcodes = resolve_mtx_triplet(args["expr-prefix"])
        section("Gene expression (Matrix Market)")
        bullet("matrix:    $mtx")
        bullet("features:  $features")
        bullet("barcodes:  $barcodes")
    end

    geno_table = nothing
    if args["geno"] !== nothing
        section("Reading genotypes")
        bullet("file: $(args["geno"])")
        geno_table = readtable(args["geno"])
        donor_col in names(geno_table) || error("--geno is missing donor-id column '$donor_col'")
    end

    pos_map = nothing
    if args["positions"] !== nothing
        pos_df = readtable(args["positions"])
        ncol(pos_df) >= 2 || error("--positions file must have at least two columns: variant, position")
        pos_map = Dict(zip(pos_df[:, 1], pos_df[:, 2]))
    end

    # The model formulas are gene-independent: build and report them once.
    models = map(effects) do e
        f_e, tt_e = build_formula(covariates, e, effect_iw[e])
        (effect = e, f = f_e, termtest = tt_e)
    end
    section(length(models) > 1 ? "Models" : "Model")
    for mdl in models
        bullet("[$(mdl.effect)] formula: $(mdl.f)")
        bullet("testing term(s): $(mdl.termtest isa AbstractVector ? join(mdl.termtest, ", ") : mdl.termtest)", indent = 2)
    end

    # ------------------------------ Output naming ----------------------------- #

    out_prefix = args["out"]
    endswith(out_prefix, ".tsv") && (out_prefix = out_prefix[1:end - 4])
    # With a single effect, filenames stay <out>_<gene>.tsv; with several,
    # the effect is included: <out>_<effect>_<gene>.tsv.
    outpath_for(glabel, effect) = begin
        tag = length(effects) > 1 ? "$(effect)_$(glabel)" : glabel
        isempty(out_prefix) ? "$(tag).tsv" :
            endswith(out_prefix, "/") ? out_prefix * tag * ".tsv" :
            out_prefix * "_" * tag * ".tsv"
    end
    outdir = dirname(out_prefix)
    isempty(outdir) || mkpath(outdir)
    endswith(out_prefix, "/") && mkpath(out_prefix)

    # ------------------------------- Check mode -------------------------------- #

    # --check: validate everything a batch job will need -- in seconds, without
    # mapping -- so a typo fails one interactive command instead of a hundred
    # submitted jobs. Hard failures above (missing files, bad meta columns,
    # malformed bed) have already surfaced by this point; here the per-gene
    # and cross-file consistency checks run, accumulating every problem found.
    if args["check"]
        section("Check mode (--check): validating inputs without mapping")
        issues = String[]

        # Expression: every bed gene must be resolvable
        if expr_table === nothing
            feats = Dynema.read_feature_fields(features)
            isfile(Dynema.dgx_sidecar_path(mtx)) ||
                bullet("note: no .dgx index next to $mtx -- each gene will trigger a full matrix scan; run dynema-prepare-expr once")
            for g in (genes === nothing ? [] : genes)
                try
                    Dynema.match_feature_row(feats, g.gene; features_path = String(features))
                catch err
                    push!(issues, "expression: $(sprint(showerror, err))")
                end
            end
        elseif genes !== nothing
            expr_gene_cols = setdiff(names(expr_table), [id_col])
            for g in genes
                nhits = count(c -> Dynema.stripver(c) == Dynema.stripver(g.gene), expr_gene_cols)
                nhits == 1 ||
                    push!(issues, "expression: gene '$(g.gene)' matches $nhits column(s) of --expr")
            end
        end

        # Genotypes: index, chromosomes, and donor coverage
        if args["vcf"] !== nothing
            try
                isfile(args["vcf"]) || error("VCF not found: $(args["vcf"])")
                (isfile(args["vcf"] * ".tbi") || isfile(args["vcf"] * ".csi")) ||
                    error("no tabix index (.tbi/.csi) next to $(args["vcf"])")
                for c in unique([g.chr for g in genes])
                    try
                        Dynema.verify_chr(args["vcf"], c)
                    catch err
                        push!(issues, "vcf: $(sprint(showerror, err))")
                    end
                end
                vsamples = Dynema.vcf_samples(args["vcf"])
                donor_set = if args["samples"] !== nothing
                    smap = readtable(args["samples"])
                    for c in ("vcf_id", "donor_id")
                        c in names(smap) || error("--samples is missing required column '$c'")
                    end
                    Set(smap.donor_id[map(v -> v in Set(vsamples), smap.vcf_id)])
                else
                    Set(vsamples)
                end
                missing_donors = setdiff(unique(meta[:, donor_col]), donor_set)
                isempty(missing_donors) ||
                    push!(issues, "donors: $(length(missing_donors)) --meta donor(s) have no genotype " *
                                  "(e.g. $(join(first(collect(missing_donors), 3), ", "))); check --vcf/--samples")
            catch err
                push!(issues, "vcf: $(sprint(showerror, err))")
            end
        else
            missing_donors = setdiff(unique(meta[:, donor_col]), geno_table[:, donor_col])
            isempty(missing_donors) ||
                push!(issues, "donors: $(length(missing_donors)) --meta donor(s) missing from --geno " *
                              "(e.g. $(join(first(collect(missing_donors), 3), ", ")))")
        end

        if isempty(issues)
            bullet("all checks passed -- ready to map $(genes === nothing ? 1 : length(genes)) gene(s) x $(length(effects)) effect(s)")
            return
        end
        for i in issues
            bullet("PROBLEM: $i")
        end
        error("--check found $(length(issues)) problem(s); nothing was mapped")
    end

    # --------------------------- Resume support helpers ------------------------- #

    # Summary row for a gene x effect skipped via --skip-existing: lead
    # statistics are re-read from the existing output so a resumed batch still
    # produces a complete summary table.
    function skipped_row(glabel, g, effect, path)
        base = (gene = glabel, effect = effect, status = "skipped",
                chr = g === nothing ? missing : g.chr,
                tss = g === nothing ? missing : g.tss)
        try
            df = readtable(path)
            li = argmin(df.p)
            statcol_i = findfirst(c -> c in ("z", "χ²"), names(df))
            statcol = statcol_i === nothing ? nothing : names(df)[statcol_i]
            return merge(base, (n_variants = nrow(df),
                lead_variant = df.variant[li],
                lead_pos = "pos" in names(df) ? df.pos[li] : missing,
                stat_type = statcol === nothing ? missing : statcol,
                lead_stat = statcol === nothing ? missing : df[li, statcol],
                lead_p = df.p[li],
                lead_p_boot = "p_boot" in names(df) ? df.p_boot[li] : missing,
                lead_p_boot_approx = "p_boot_approx" in names(df) ? df.p_boot_approx[li] : missing,
                out_file = path))
        catch
            return merge(base, (n_variants = missing, lead_variant = missing, lead_pos = missing,
                stat_type = missing, lead_stat = missing, lead_p = missing,
                lead_p_boot = missing, lead_p_boot_approx = missing, out_file = path))
        end
    end

    # -------------------------------- Map one gene ---------------------------- #

    function map_gene(g)

        # --skip-existing: when every requested effect's output already
        # exists, skip the gene before paying for any extraction.
        if args["skip-existing"] && g !== nothing &&
           all(isfile(outpath_for(g.gene, mdl.effect)) for mdl in models)
            bullet("all output file(s) already exist; skipping (--skip-existing)")
            return [skipped_row(g.gene, g, mdl.effect, outpath_for(g.gene, mdl.effect)) for mdl in models]
        end

        gene = g === nothing ? nothing : g.gene

        # Genotypes first: --vcf extraction is normally much faster than an
        # unindexed Matrix Market scan, so genotype failures surface early.
        geno_t0 = time()
        vcf_pos_map = nothing
        geno_df = if args["vcf"] !== nothing
            bullet("extracting genotypes from $(args["vcf"])")
            r = extract_geno_dataframe(
                vcf = args["vcf"],
                chr = g.chr,
                tss = g.tss,
                window = args["window"],
                field = args["field"],
                samples_file = args["samples"],
                donor_col = donor_col,
                maf = args["maf"],
                max_missing = args["max-missing"],
                verbose = false,
            )
            bullet("cis-window: $(r.chr):$(r.start_pos)-$(r.end_pos) (TSS $(r.tss) +/- $(args["window"]) bp, from --bed)", indent = 2)
            bullet("samples: $(r.n_samples_vcf) in VCF, $(r.n_samples_matched) retained after sample matching", indent = 2)
            bullet("variants in region: $(r.n_variants_total); retained: $(r.n_retained) " *
                   "(skipped: $(r.n_multiallelic) multiallelic, $(r.n_no_field) missing GP/DS, " *
                   "$(r.n_high_missing) high-missingness)", indent = 2)
            vcf_pos_map = Dict(zip(names(r.geno)[2:end], r.positions))
            r.geno
        else
            geno_table
        end
        donor_col in names(geno_df) || error("Genotype table is missing donor-id column '$donor_col'")
        ncol(geno_df) == 1 &&
            error("No variants in the genotype table (0 columns besides '$donor_col'); nothing to test. " *
                  "Check --bed/--window (and that the annotation's chromosome naming matches the VCF's), " *
                  "and --maf/--max-missing.")
        bullet("genotypes ready in $(elapsed_str(geno_t0))")

        expr_t0 = time()
        expr_df = if expr_table === nothing
            r = extract_gene_expression(
                mtx = mtx,
                features = features,
                barcodes = barcodes,
                gene = gene,
                id_col = id_col,
                verbose = false,
            )
            bullet("gene '$gene' found at expression row $(r.target_row) of $(r.n_genes); " *
                   "$(r.n_found) nonzero entries across $(r.n_cells) cell(s)")
            r.expr
        else
            expr_table
        end
        id_col in names(expr_df) || error("Expression table is missing cell-id column '$id_col'")

        gene_cols = setdiff(names(expr_df), [id_col])
        isempty(gene_cols) && error("Expression table has no gene columns besides '$id_col'")
        gene_col = if gene === nothing
            length(gene_cols) == 1 ? gene_cols[1] :
                error("Expression table has multiple gene columns ($(join(first(gene_cols, 10), ", "))...); provide --bed naming the gene(s)")
        else
            hits = filter(c -> Dynema.stripver(c) == Dynema.stripver(gene), gene_cols)
            isempty(hits) && error("Gene '$gene' not found in the expression table")
            length(hits) > 1 &&
                error("Gene '$gene' matches multiple expression columns ($(join(hits, ", ")))")
            hits[1]
        end
        gene_label = gene === nothing ? gene_col : gene
        bullet("expression ready in $(elapsed_str(expr_t0))")

        # ----------------------- Align expression to metadata ------------------ #

        expr_lookup = Dict(zip(expr_df[:, id_col], expr_df[:, gene_col]))
        missing_expr = setdiff(meta[:, id_col], expr_df[:, id_col])
        isempty(missing_expr) ||
            error("$(length(missing_expr)) cell(s) in --meta have no matching row in --expr (e.g. $(first(missing_expr)))")
        pheno = Float64.([expr_lookup[cid] for cid in meta[:, id_col]])

        # ------------------------ Expand genotypes to cells -------------------- #

        donor_ids = geno_df[:, donor_col]
        snp_cols  = isempty(variant_filter) ? setdiff(names(geno_df), [donor_col]) : variant_filter
        for v in snp_cols
            v in names(geno_df) || error("Variant '$v' not found in the genotype table (--geno or --vcf extraction)")
        end

        missing_donors = setdiff(unique(meta[:, donor_col]), donor_ids)
        isempty(missing_donors) ||
            error("$(length(missing_donors)) donor(s) in --meta have no matching genotype (e.g. $(first(missing_donors))); check --geno/--vcf and --samples")

        geno_mat = Matrix(geno_df[:, snp_cols])
        ex_geno = expand_genotypes(geno_mat, donor_ids, meta[:, donor_col], snp_cols)

        # Positions: with --vcf they come from the VCF automatically; an
        # explicit --positions file overrides them (and is the only source
        # for --geno). Resolved once per gene, shared across effects.
        effective_pos_map = pos_map !== nothing ? pos_map : vcf_pos_map
        gene_pos = nothing
        if effective_pos_map !== nothing
            pos = get.(Ref(effective_pos_map), snp_cols, missing)
            if any(ismissing, pos)
                missing_pos = snp_cols[ismissing.(pos)]
                @warn "No position found for $(length(missing_pos)) variant(s) (e.g. $(first(missing_pos))); positions not attached"
            else
                gene_pos = Int.(pos)
            end
        end

        # ------------------------- Run Dynema (per effect) ---------------------- #

        rows = NamedTuple[]
        for mdl in models
            tag = length(models) > 1 ? "[$(mdl.effect)] " : ""
            outpath = outpath_for(gene_label, mdl.effect)
            if args["skip-existing"] && isfile(outpath)
                bullet("$(tag)$outpath exists; skipping (--skip-existing)")
                push!(rows, skipped_row(gene_label, g, mdl.effect, outpath))
                continue
            end
            try
                bullet("$(tag)mapping $(length(snp_cols)) variant(s)...")
                res = map_locus(mdl.f;
                    pheno = pheno,
                    geno = ex_geno,
                    meta = meta,
                    groups = meta[:, donor_col],
                    termtest = mdl.termtest,
                    parallel = args["parallel"] || args["workers"] > 0,
                    betas = Symbol(args["betas"]),
                    boot = args["boot"],
                    B = B,
                    ptype = ptype,
                    gene = gene_label,
                    chr = g === nothing ? nothing : g.chr,
                )
                gene_pos === nothing || set_pos!(res, gene_pos)

                # Dynema only defines the MIME"text/plain" show method (used
                # for REPL auto-display), so invoke it explicitly to get the
                # pretty-printed summary instead of Julia's generic fallback.
                show(IOContext(stdout, :displaysize => term_size), MIME("text/plain"), res)
                println()

                summ = copy(get_summary(res))
                if get_pos(res) !== nothing
                    insertcols!(summ, 2, :pos => get_pos(res))
                    g === nothing || insertcols!(summ, 2, :chr => fill(g.chr, nrow(summ)))
                end
                CSV.write(outpath, summ; delim = '\t')
                bullet("$(tag)wrote summary statistics for $(nrow(summ)) variant(s) to $outpath")

                # One summary row per gene x effect, for the batch summary table.
                li = argmin(summ.p)
                statcol = get_stattype(res)
                push!(rows, (gene = gene_label,
                    effect = mdl.effect,
                    status = "ok",
                    chr = g === nothing ? missing : g.chr,
                    tss = g === nothing ? missing : g.tss,
                    n_variants = nrow(summ),
                    lead_variant = summ.variant[li],
                    lead_pos = "pos" in names(summ) ? summ.pos[li] : missing,
                    stat_type = statcol,
                    lead_stat = summ[li, statcol],
                    lead_p = summ.p[li],
                    lead_p_boot = "p_boot" in names(summ) ? summ.p_boot[li] : missing,
                    lead_p_boot_approx = "p_boot_approx" in names(summ) ? summ.p_boot_approx[li] : missing,
                    out_file = outpath))
            catch err
                # Single gene + single effect: fail loudly (rethrown again by
                # the batch loop); otherwise one failed effect must not sink
                # the others.
                (length(batch) == 1 && length(models) == 1) && rethrow()
                bt = catch_backtrace()
                @warn "Mapping failed for gene $gene_label, effect $(mdl.effect); continuing" exception = (err, bt)
                push!(rows, (gene = gene_label, effect = mdl.effect, status = "failed",
                    chr = g === nothing ? missing : g.chr,
                    tss = g === nothing ? missing : g.tss,
                    n_variants = missing, lead_variant = missing, lead_pos = missing,
                    stat_type = missing, lead_stat = missing, lead_p = missing,
                    lead_p_boot = missing, lead_p_boot_approx = missing, out_file = missing))
            end
        end

        return rows

    end

    # ---------------------------- Map the whole batch -------------------------- #

    batch = genes === nothing ? [nothing] : genes
    n_ok = 0
    failed = String[]
    gene_rows = NamedTuple[]
    batch_t0 = time()

    for (gi, g) in enumerate(batch)
        glabel = g === nothing ? "single --expr gene" : g.gene
        section(length(batch) > 1 ? "Gene $glabel [$gi/$(length(batch))]" : "Gene $glabel")
        try
            rows = map_gene(g)
            append!(gene_rows, rows)
            all(r -> r.status in ("ok", "skipped"), rows) ? (n_ok += 1) : push!(failed, glabel)
        catch err
            # For a single-gene run, fail loudly (as before); in a batch, one
            # bad gene (e.g. absent from the features file) must not sink the
            # rest of the job.
            length(batch) == 1 && rethrow()
            bt = catch_backtrace()
            @warn "Mapping failed for gene $glabel; continuing with the remaining genes" exception = (err, bt)
            push!(failed, glabel)
            for mdl in models
                push!(gene_rows, (gene = glabel, effect = mdl.effect, status = "failed",
                    chr = g === nothing ? missing : g.chr,
                    tss = g === nothing ? missing : g.tss,
                    n_variants = missing, lead_variant = missing, lead_pos = missing,
                    stat_type = missing, lead_stat = missing, lead_p = missing,
                    lead_p_boot = missing, lead_p_boot_approx = missing, out_file = missing))
            end
        end
    end

    # ------------------------- Batch-level gene summary ------------------------- #

    # One row per gene (lead variant and its statistics), so a batch's -- or,
    # concatenated across chunks, a whole study's -- top associations can be
    # inspected without parsing every per-gene table.
    summary_path = isempty(out_prefix) ? "dynema_summary.tsv" :
        endswith(out_prefix, "/") ? out_prefix * "summary.tsv" : out_prefix * "_summary.tsv"
    CSV.write(summary_path, DataFrame(gene_rows); delim = '\t')

    section("Done")
    bullet("$n_ok/$(length(batch)) gene(s) fully mapped ($(join(effects, ", "))) in $(elapsed_str(batch_t0))")
    n_skipped = count(r -> r.status == "skipped", gene_rows)
    n_skipped > 0 && bullet("$n_skipped gene x effect output(s) skipped (--skip-existing)")
    isempty(failed) || bullet("gene(s) with failures: $(join(failed, ", "))")
    bullet("per-gene summary (lead variants): $summary_path")

end

"""
    main()

Parses CLI arguments, prints the invoking command, then runs `run_map(args)`
wrapped in `with_tee_log` so the whole run's console output -- and that same
command -- is saved to `--log` (default: the --out prefix with `.log`
appended, or `dynema_map.log`) in addition to being shown on the terminal as
usual.
"""
function main()

    args = parse_commandline()

    out_stem = args["out"]
    endswith(out_stem, ".tsv") && (out_stem = out_stem[1:end - 4])
    log_path = args["log"] !== nothing ? args["log"] :
        isempty(out_stem)          ? "dynema_map.log" :
        endswith(out_stem, "/")    ? out_stem * "dynema_map.log" :
                                     out_stem * ".log"
    command  = "dynema_map.jl " * join(shell_quote.(ARGS), " ")

    section("Command")
    bullet(command)

    # Captured before with_tee_log redirects stdout to a (non-TTY) pipe, so
    # the results table below still wraps at the real terminal's width.
    term_size = displaysize(stdout)

    with_tee_log(log_path; command = command) do
        run_map(args; term_size = term_size)
    end

end

main()
