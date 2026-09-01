#!/usr/bin/env julia
# ---------------------------------------------------------------------------- #
#                                 dynema_map.jl                                #
# ---------------------------------------------------------------------------- #
#
# Command-line wrapper around `Dynema.map_locus` for users who don't want to
# write any Julia code. Reads single-cell gene expression either from a plain
# --expr TSV/CSV or extracted on the fly for one gene from an --expr-prefix
# Matrix Market triplet (see mtx_expr.jl), reads single-cell metadata from a
# plain TSV/CSV, gets donor-level genotype dosages either from a
# pre-extracted --geno matrix or extracted on the fly for the tested gene's
# cis-window from a --vcf (see vcf_geno.jl), builds the eQTL model formula
# from a handful of flags, runs Dynema for one gene against one or more
# genetic variants, and writes a summary statistics table. The full console
# transcript, plus the exact command run, is also saved to --log (default:
# --out with its extension swapped for .log) -- see cli_output.jl.
#
# Usage (pre-extracted genotype matrix):
#   julia --project=<this-dir> dynema_map.jl --geno geno.tsv [options]
#
# or extracting genotypes on the fly from a tabix-indexed VCF:
#   julia --project=<this-dir> dynema_map.jl --vcf genotypes.vcf.gz \
#       --gene BRCA1 --tss-file genes.tsv --window 250000 [options]
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
# calls are fast no-ops once the manifest already satisfies Project.toml.
!isfile(joinpath(CLI_DIR, "Manifest.toml")) &&
    println("First run detected: setting up the shared bin/ environment (this can take a few minutes)...")
Pkg.develop(Pkg.PackageSpec(path = joinpath(CLI_DIR, "..")); io = devnull)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using StatsModels
using CategoricalArrays
using Dynema

include(joinpath(CLI_DIR, "cli_output.jl")) # section, bullet, elapsed_str -- used here and by the includes below
include(joinpath(CLI_DIR, "vcf_geno.jl")) # extract_geno_dataframe, for --vcf on-the-fly extraction
include(joinpath(CLI_DIR, "mtx_expr.jl"))  # resolve_mtx_triplet, extract_gene_expression, for --expr-prefix on-the-fly extraction

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
            help = "TSV/CSV file with single-cell gene expression counts. Must contain a cell-id column and one column per gene (or a single expression column if --gene is omitted). Exactly one of --expr or --expr-prefix must be given."
            arg_type = String
        "--expr-prefix"
            help = "Filename prefix of a Matrix Market gene expression export (genes as rows, cells as columns): resolves to <prefix>.mtx, <prefix>.features, <prefix>.barcodes, each tried gzipped first then plain (e.g. 'expr_0.05' -> expr_0.05.mtx.gz/.features.gz/.barcodes.gz). Streamed for just the tested gene, never loaded in full. Requires --gene. Exactly one of --expr or --expr-prefix must be given."
            default = nothing
        "--gene"
            help = "Name of the gene to test: the --expr column to use (required if --expr has more than one non-id column), the row to look up in --expr-prefix's features file (required with --expr-prefix), and/or the id to look up in --tss-file for --vcf (only used, and only required, if --tss-file is given instead of --chr/--tss directly)."
            default = nothing
        "--geno"
            help = "TSV/CSV file with donor-level genotype dosages: one donor-id column plus one column per variant. Exactly one of --geno or --vcf must be given."
            arg_type = String
        "--vcf"
            help = "Bgzipped, tabix-indexed VCF (.vcf.gz) to extract the gene's cis-window genotypes from on the fly, instead of a pre-extracted --geno file. Uses --gene (or --chr/--tss) for the TSS, plus --tss-file/--window/--field/--samples/--maf/--max-missing below. Exactly one of --geno or --vcf must be given."
            default = nothing
        "--tss-file"
            help = "(--vcf only) TSV/CSV with columns 'gene_id','chr','tss' to look up --gene's TSS from. Alternative to passing --chr/--tss directly."
            default = nothing
        "--tss"
            help = "(--vcf only) Transcription start site position (1-based), if not using --tss-file."
            arg_type = Int
        "--window"
            help = "(--vcf only) Cis-window half-width in bp around the TSS."
            arg_type = Int
            default = 1_000_000
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
            help = "Comma-separated list of cell-state/context column names in --meta (e.g. C1,C2,C3). Always included as covariates and available for interaction testing."
            arg_type = String
            default = ""
        "--test"
            help = "Effect to test: 'main' (context-independent), 'interaction' (single- or multi-context), or 'total' (main + interaction jointly)."
            arg_type = String
            default = "main"
            range_tester = x -> x in ("main", "interaction", "total")
        "--interaction-with"
            help = "Comma-separated subset of --contexts to test a G-by-context interaction for. Required for --test interaction/total unless --contexts is set, in which case it defaults to all of --contexts."
            default = nothing
        "--boot"
            help = "Also compute empirical p-values via adaptive score bootstrapping (recommended for small or imbalanced cohorts)."
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
        "--wald"
            help = "Use a Wald test at the unrestricted MLE instead of the default Lagrange-multiplier/score test."
            action = :store_true
        "--parallel"
            help = "Distribute variants across worker processes with Distributed.jl (workers must already be started, e.g. via `julia -p 4`)."
            action = :store_true
        "--positions"
            help = "Optional TSV/CSV with two columns (variant, position) giving genomic positions to attach to the output."
            default = nothing
        "--chr"
            help = "Chromosome label attached to the output. Also doubles as the VCF query chromosome for --vcf (with --tss), if not using --tss-file."
            default = nothing
        "--out"
            help = "Output path for the summary statistics table (TSV)."
            arg_type = String
            required = true
        "--log"
            help = "Path to save a full transcript of this run (the exact command invoked, plus everything printed to the console -- progress, warnings, and the results table). Default: --out with its extension replaced by .log."
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
    build_formula(covariates, contexts, test, interaction_with)

Builds a `FormulaTerm` and the `termtest` argument for `map_locus`
programmatically (equivalent to writing e.g. `@formula(0 ~ 1 + G + C1 + G & C1)`
by hand), from plain lists of column names supplied on the command line.

Any context named in `interaction_with` that isn't also in `contexts` is
still given a main-effect term in the formula (a `G & context` interaction
term without that context's own main effect would be an unusual, generally
statistically improper model), and a warning is printed -- this most often
means `--contexts`/`--interaction-with` have a typo relative to each other.
"""
function build_formula(covariates::Vector{String}, contexts::Vector{String},
                        test::String, interaction_with::Vector{String})

    extra_contexts = setdiff(interaction_with, contexts)
    if !isempty(extra_contexts) && test in ("interaction", "total")
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

    termtest = String[]

    if test in ("main", "total")
        push!(termtest, "G")
    end

    if test in ("interaction", "total")
        isempty(interaction_with) &&
            error("--test $test requires --interaction-with (or a non-empty --contexts to default to)")
        for c in interaction_with
            push!(rhs_terms, term(:G) & term(Symbol(c)))
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

    covariates       = splitcsv(args["covariates"])
    contexts         = splitcsv(args["contexts"])
    interaction_with = args["interaction-with"] === nothing ? contexts : splitcsv(args["interaction-with"])
    variant_filter   = splitcsv(args["variants"])
    B                = parse.(Int, splitcsv(args["B"]))
    ptype            = Symbol(args["ptype"])
    id_col           = args["cell-id-col"]
    donor_col        = args["donor-col"]

    # -------------------------------- Load data ------------------------------ #

    gene = args["gene"]

    # Validate --meta's columns (cell id, donor id, contexts, covariates) up
    # front, before the potentially slow --expr/--expr-prefix read below --
    # so a typo'd column name fails fast instead of after a multi-minute
    # Matrix Market scan.
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
    bullet("$(nrow(meta)) cell(s)")

    expr_t0 = time()
    expr_df = if args["expr-prefix"] !== nothing
        gene === nothing &&
            error("--expr-prefix requires --gene (there's no header row to infer a single gene from)")
        mtx, features, barcodes = resolve_mtx_triplet(args["expr-prefix"])
        extract_gene_expression(
            mtx = mtx,
            features = features,
            barcodes = barcodes,
            gene = gene,
            id_col = id_col,
        )
    else
        section("Reading gene expression")
        bullet("file: $(args["expr"])")
        readtable(args["expr"])
    end
    bullet("done in $(elapsed_str(expr_t0))")
    id_col in names(expr_df) || error("Expression table is missing cell-id column '$id_col'")

    gene_cols = setdiff(names(expr_df), [id_col])
    isempty(gene_cols) && error("Expression table has no gene columns besides '$id_col'")
    if gene === nothing
        gene = length(gene_cols) == 1 ? gene_cols[1] :
            error("Expression table has multiple gene columns ($(join(gene_cols, ", "))); specify --gene")
    end
    gene in names(expr_df) || error("Gene '$gene' not found in the expression table")

    geno_t0 = time()
    geno_df = if args["vcf"] !== nothing
        extract_geno_dataframe(
            vcf = args["vcf"],
            chr = args["chr"],
            tss = args["tss"],
            tss_file = args["tss-file"],
            gene = args["tss-file"] === nothing ? nothing : gene,
            window = args["window"],
            field = args["field"],
            samples_file = args["samples"],
            donor_col = donor_col,
            maf = args["maf"],
            max_missing = args["max-missing"],
        )
    else
        section("Reading genotypes")
        bullet("file: $(args["geno"])")
        readtable(args["geno"])
    end
    bullet("done in $(elapsed_str(geno_t0))")
    donor_col in names(geno_df) || error("Genotype table is missing donor-id column '$donor_col'")

    # ----------------------- Align expression to metadata --------------------- #

    expr_lookup = Dict(zip(expr_df[:, id_col], expr_df[:, gene]))
    missing_expr = setdiff(meta[:, id_col], expr_df[:, id_col])
    isempty(missing_expr) ||
        error("$(length(missing_expr)) cell(s) in --meta have no matching row in --expr (e.g. $(first(missing_expr)))")
    pheno = Float64.([expr_lookup[cid] for cid in meta[:, id_col]])

    # ------------------------ Expand genotypes to cells ------------------------ #

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

    # --------------------------- Build model formula --------------------------- #

    f, termtest = build_formula(covariates, contexts, args["test"], interaction_with)
    section("Model")
    bullet("formula: $f")
    bullet("testing term(s): $(termtest isa AbstractVector ? join(termtest, ", ") : termtest)")
    bullet("N cells = $(nrow(meta)), N donors = $(length(unique(meta[:, donor_col]))), N variants = $(length(snp_cols))")

    # -------------------------------- Run Dynema -------------------------------- #

    section("Running Dynema")
    res = map_locus(f;
        pheno = pheno,
        geno = ex_geno,
        meta = meta,
        groups = meta[:, donor_col],
        termtest = termtest,
        parallel = args["parallel"],
        imposenull = !args["wald"],
        boot = args["boot"],
        B = B,
        ptype = ptype,
        gene = gene,
        chr = args["chr"],
    )

    # ------------------------- Optional genomic positions ----------------------- #

    if args["positions"] !== nothing
        pos_df = readtable(args["positions"])
        ncol(pos_df) >= 2 || error("--positions file must have at least two columns: variant, position")
        pos_map = Dict(zip(pos_df[:, 1], pos_df[:, 2]))
        pos = get.(Ref(pos_map), snp_cols, missing)
        if any(ismissing, pos)
            missing_pos = snp_cols[ismissing.(pos)]
            @warn "No position found for $(length(missing_pos)) variant(s) (e.g. $(first(missing_pos))); positions not attached"
        else
            set_pos!(res, Int.(pos))
        end
    end

    # Dynema only defines the MIME"text/plain" show method (used for REPL
    # auto-display), so invoke it explicitly to get the pretty-printed summary
    # instead of Julia's generic struct fallback.
    section("Results")
    show(IOContext(stdout, :displaysize => term_size), MIME("text/plain"), res)
    println()

    # ------------------------------- Write output ------------------------------- #

    summ = get_summary(res)
    CSV.write(args["out"], summ; delim = '\t')
    section("Done")
    bullet("wrote summary statistics for $(nrow(summ)) variant(s) to $(args["out"])")

end

"""
    main()

Parses CLI arguments, prints the invoking command, then runs `run_map(args)`
wrapped in `with_tee_log` so the whole run's console output -- and that same
command -- is saved to `--log` (default: `--out` with its extension
replaced by `.log`) in addition to being shown on the terminal as usual.
"""
function main()

    args = parse_commandline()

    log_path = args["log"] !== nothing ? args["log"] : splitext(args["out"])[1] * ".log"
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
