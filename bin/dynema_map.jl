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
        "--contexts"
            help = "Comma-separated list of cell-state/context column names in --meta (e.g. C1,C2,C3). Always included as covariates. Required (and must be non-empty) when --effect is interaction/total, alongside --interaction-with."
            arg_type = String
            default = ""
        "--effect"
            help = "Effect to test: 'main' (context-independent), 'interaction' (single- or multi-context), or 'total' (main + interaction jointly). Always required -- no default -- to avoid silently testing the wrong effect."
            arg_type = String
            required = true
            range_tester = x -> x in ("main", "interaction", "total")
        "--interaction-with"
            help = "Comma-separated subset of --contexts to include as G-by-context interaction terms in the model. Required (and must be non-empty), alongside a non-empty --contexts, when --effect is interaction/total -- neither is ever defaulted or guessed. Optional with --effect main: if given there, the interaction terms are still added to the formula (e.g. to adjust for/include an interaction without testing it), just not tested."
            default = nothing
        "--boot"
            help = "Also compute bootstrap p-values via adaptive score bootstrapping (recommended for small or imbalanced cohorts). Adds two columns: p_boot, the empirical bootstrap p-value (floored at ~2/B), and p_boot_approx, a FastQTL-style beta approximation fitted to the bootstrap distribution that extrapolates smoothly below that floor."
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
        "--no-betas"
            help = "Skip fitting the unrestricted model per variant and omit its coefficient estimates from the output. The score test does not need those estimates -- they are only reported for convenience -- so this substantially reduces per-variant fitting work."
            action = :store_true
        "--parallel"
            help = "Distribute variants across worker processes with Distributed.jl. Use with --workers N, or start workers yourself (e.g. `julia -p 4 --project=bin bin/dynema_map.jl ...`)."
            action = :store_true
        "--workers"
            help = "Number of local worker processes to start for parallel mapping. Implies --parallel."
            arg_type = Int
            default = 0
        "--positions"
            help = "Optional TSV/CSV with two columns (variant, position) giving genomic positions to attach to the output."
            default = nothing
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
    resolve_interactions(effect, contexts, interaction_with_arg) -> Vector{String}

Resolves the contexts that should actually appear as `G & context` terms in
the model formula (passed to `build_formula`, and used to validate --meta's
columns before that):

- If `--effect` is `interaction`/`total` (which always test at least one
  `G & context` term), both `--contexts` and `--interaction-with` must be
  given explicitly and non-empty -- errors immediately otherwise, rather
  than guessing a set of interactions or silently modeling a G-by-context
  interaction with no declared context main effects.
- If `--effect` is `main`, `--interaction-with` is optional: given
  explicitly, its contexts are still added to the formula (e.g. so
  `--effect main --interaction-with CV1` can adjust for/include a G x CV1
  interaction in the model without testing it); omitted, no interaction
  terms at all are added.
"""
function resolve_interactions(effect::String, contexts::Vector{String},
                               interaction_with_arg::Union{Nothing,Vector{String}})

    needs_interaction = effect in ("interaction", "total")

    if needs_interaction
        isempty(contexts) &&
            error("--effect $effect requires a non-empty --contexts (the context(s) to test G-by-context interactions for)")
        (interaction_with_arg === nothing || isempty(interaction_with_arg)) &&
            error("--effect $effect requires --interaction-with (a non-empty, comma-separated subset of --contexts)")
    end

    return interaction_with_arg === nothing ? String[] : interaction_with_arg

end

"""
    build_formula(covariates, contexts, effect, interaction_with)

Builds a `FormulaTerm` and the `termtest` argument for `map_locus`
programmatically (equivalent to writing e.g. `@formula(0 ~ 1 + G + C1 + G & C1)`
by hand), from plain lists of column names supplied on the command line.
`interaction_with` (see `resolve_interactions`) is always added to the
formula as `G & context` terms when non-empty; only added to `termtest`
(i.e. actually tested, alongside/instead of `G`) when `effect` is
`interaction`/`total`.

Any context named in `interaction_with` that isn't also in `contexts` is
still given a main-effect term in the formula (a `G & context` interaction
term without that context's own main effect would be an unusual, generally
statistically improper model), and a warning is printed -- this most often
means `--contexts`/`--interaction-with` have a typo relative to each other.
"""
function build_formula(covariates::Vector{String}, contexts::Vector{String},
                        effect::String, interaction_with::Vector{String})

    extra_contexts = setdiff(interaction_with, contexts)
    if !isempty(extra_contexts)
        @warn "--interaction-with includes context(s) not in --contexts ($(join(extra_contexts, ", "))); " *
              "adding them as main-effect covariates too, since an interaction term needs its main effect " *
              "in the model. Double check --contexts/--interaction-with for a typo if this is unexpected."
    end
    all_contexts = union(contexts, interaction_with)

    rhs_terms = Any[term(1), term(:G)]
    for c in covariates
        push!(rhs_terms, term(Symbol(c)))
    end
    for c in all_contexts
        push!(rhs_terms, term(Symbol(c)))
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
    contexts             = splitcsv(args["contexts"])
    interaction_with_arg = args["interaction-with"] === nothing ? nothing : splitcsv(args["interaction-with"])
    # Fails immediately (before any file loading below) if --effect
    # interaction/total was given without a non-empty --interaction-with --
    # see resolve_interactions's docstring.
    interaction_with     = resolve_interactions(args["effect"], contexts, interaction_with_arg)
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
    for c in vcat(covariates, union(contexts, interaction_with))
        c in names(meta) || error("Column '$c' not found in --meta")
    end
    bullet("cell-id column: $id_col, donor-id column: $donor_col")
    isempty(covariates) || bullet("covariates: $(join(covariates, ", "))")
    isempty(contexts) || bullet("contexts: $(join(contexts, ", "))")
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

    # The model formula is gene-independent: build and report it once.
    f, termtest = build_formula(covariates, contexts, args["effect"], interaction_with)
    section("Model")
    bullet("formula: $f")
    bullet("testing term(s): $(termtest isa AbstractVector ? join(termtest, ", ") : termtest)")

    # ------------------------------ Output naming ----------------------------- #

    out_prefix = args["out"]
    endswith(out_prefix, ".tsv") && (out_prefix = out_prefix[1:end - 4])
    outpath_for(glabel) = isempty(out_prefix) ? "$(glabel).tsv" :
        endswith(out_prefix, "/") ? out_prefix * glabel * ".tsv" :
        out_prefix * "_" * glabel * ".tsv"
    outdir = dirname(out_prefix)
    isempty(outdir) || mkpath(outdir)
    endswith(out_prefix, "/") && mkpath(out_prefix)

    # -------------------------------- Map one gene ---------------------------- #

    function map_gene(g)

        gene = g === nothing ? nothing : g.gene

        # Genotypes first: --vcf extraction is normally much faster than an
        # unindexed Matrix Market scan, so genotype failures surface early.
        geno_t0 = time()
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

        # -------------------------------- Run Dynema --------------------------- #

        bullet("mapping $(length(snp_cols)) variant(s)...")
        res = map_locus(f;
            pheno = pheno,
            geno = ex_geno,
            meta = meta,
            groups = meta[:, donor_col],
            termtest = termtest,
            parallel = args["parallel"] || args["workers"] > 0,
            betas = !args["no-betas"],
            boot = args["boot"],
            B = B,
            ptype = ptype,
            gene = gene_label,
            chr = g === nothing ? nothing : g.chr,
        )

        # --------------------- Optional genomic positions ---------------------- #

        if pos_map !== nothing
            pos = get.(Ref(pos_map), snp_cols, missing)
            if any(ismissing, pos)
                missing_pos = snp_cols[ismissing.(pos)]
                @warn "No position found for $(length(missing_pos)) variant(s) (e.g. $(first(missing_pos))); positions not attached"
            else
                set_pos!(res, Int.(pos))
            end
        end

        # Dynema only defines the MIME"text/plain" show method (used for REPL
        # auto-display), so invoke it explicitly to get the pretty-printed
        # summary instead of Julia's generic struct fallback.
        show(IOContext(stdout, :displaysize => term_size), MIME("text/plain"), res)
        println()

        # ------------------------------ Write output --------------------------- #

        summ = get_summary(res)
        outpath = outpath_for(gene_label)
        CSV.write(outpath, summ; delim = '\t')
        bullet("wrote summary statistics for $(nrow(summ)) variant(s) to $outpath")

        return nothing

    end

    # ---------------------------- Map the whole batch -------------------------- #

    batch = genes === nothing ? [nothing] : genes
    n_ok = 0
    failed = String[]
    batch_t0 = time()

    for (gi, g) in enumerate(batch)
        glabel = g === nothing ? "single --expr gene" : g.gene
        section(length(batch) > 1 ? "Gene $glabel [$gi/$(length(batch))]" : "Gene $glabel")
        try
            map_gene(g)
            n_ok += 1
        catch err
            # For a single-gene run, fail loudly (as before); in a batch, one
            # bad gene (e.g. absent from the features file) must not sink the
            # rest of the job.
            length(batch) == 1 && rethrow()
            bt = catch_backtrace()
            @warn "Mapping failed for gene $glabel; continuing with the remaining genes" exception = (err, bt)
            push!(failed, glabel)
        end
    end

    section("Done")
    bullet("$n_ok/$(length(batch)) gene(s) mapped successfully in $(elapsed_str(batch_t0))")
    isempty(failed) || bullet("failed: $(join(failed, ", "))")

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
