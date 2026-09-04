"""

    map_locus(f::FormulaTerm; pheno::AbstractVector, geno::Union{AbstractDataFrame, AbstractVector},
              meta::AbstractDataFrame, groups::Union{AbstractDataFrame, AbstractVector},
              termtest::Union{String, Vector{String}}) -> Dynema_struct

Performs a single-cell eQTL mappping for a gene. All data arguments are defined at the single-cell level observation.

# Mandatory arguments

- `f`: Modelling formula focused on predictors. Should include all the terms to be tested and covariates. As genotype information 
can change when iterating across variants, `G` is reserved to represent genotype information for each genetic variant. For example,
a minimal formula can be defined as `0 ~ 1 + G`, where `1` is the intercept term and `G` represents the genetic variant being tested.
Gene expression can be represented with no name simply as `0`.
- `pheno`: Gene expression counts for a particular gene. A vector of length n, where _n_ is the number of cells
- `geno`: Genotype information in allele dosage format or genotype probabilities (e.g. [0, 2, 1, 0, 1]) 
expanded at the single-cell level (length _n_). If a DataFrame is provided, each column should correspond to an individual
genetic variant and all variants will be tested. The dimensions of `geno` should be (_n_ x _gv_), where _gv_ is the total number 
of variants
-  `meta`: Single-cell level metadata including single-cell cell states of interest. It also should include
single-cell and donor-level covariates expanded at the single-cell level. Dimensions _n_ x _p_, where _p_ is the
number of variables
- `groups`: A Vector or DataFrame indicating how cells should be grouped. At least one grouping variable should be 
specified (e.g. donor structure). For example ["donor1", "donor2", "donor1", "donor3"].
- `termtest`: Term included in the formula that is tested. In the case of dynamic effects, this should be the interaction
`G` and a single-cell covariate (e.g. "G & CV1")

# Optional arguments

- `parallel`: Distributes contiguous blocks of variants across worker processes
(`Distributed.jl`; start workers first, e.g. `julia -p 4`). Each worker builds its
own locus workspace, so the fast path applies in parallel too. Warm-start chains
restart at block boundaries, so interaction-test results are reproducible for a
fixed worker count but can differ at IRLS convergence tolerance (~1e-7) across
worker counts or vs. a serial run
- `H0`: Null hypothesis value. By default `0`

Testing always uses a cluster-robust (CRVE) score/Lagrange-multiplier test at
the restricted (null) fit. If the null model contains no genotype-dependent
terms (as in main/total effect tests, where every tested term involves `G`),
it is fit only once and reused across all variants.

- `betas`: Which variants get unrestricted-model coefficient estimates attached as extra summary columns (they are not
needed for the score test itself): `:lead` (default) fits only the lead variant(s) -- smallest analytical p-value,
including exact ties -- leaving the estimate columns `missing` for all other variants; `:all` fits every variant
(slowest); `:none` skips estimates entirely (fastest)
- `boot`: Apply score bootstrapping? If false, analytical p-values using a CRVE are returned. Bootstrapping
requires the optional WildBootTests package (`Pkg.add("WildBootTests"); using WildBootTests` alongside `using Dynema`)
- `B``: Number of bootstrap iterations to apply. By default 39999 iterations at apply to achieve a p-value of 5 x 10^-5 for a two-tail test. 
For adaptive bootstrapping, the number of iterations might be specified as a vector indicating the number of iterations
to perform in each step
- `ptype`: Type of bootstrap p-value to return (:equaltail, :symmetric). Only applicable when one term is being tested.
- `rboot`: Whether to return the bootstrap distributions for each variant. Useful for direct assessment of bootstrap statistics and
internal debugging
- `pos`: A numeric value specifying a genomic location for each genetic variant. Stored in final output for convenience
- `gene`: Name of the gene being tested. Stored in final output for convenience
- `chr`: Chromosome position of gene being tested. Stored in final output for convenience.
"""
function map_locus(f::FormulaTerm; pheno::AbstractVector, geno::Union{AbstractMatrix, AbstractVector}, 
                    meta::AbstractDataFrame, groups::Union{AbstractDataFrame, AbstractVector}, termtest::Union{String, Vector{String}}, 
                    parallel = false,
                    H0::Float64 = Float64(0),
                    betas::Union{Symbol, AbstractString} = :lead,
                    boot::Bool = false,
                    B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    ptype::Symbol = :equaltail, rboot = false,
                    pos::Union{Nothing, Vector{Int64}, Vector{Float64}} = nothing,
                    gene::Union{Nothing, String} = nothing,
                    chr::Union{Nothing, String, Int} = nothing)

    # ------------------------------- Betas mode -------------------------------- #

    betamode = Symbol(betas)
    betamode in (:all, :none, :lead) ||
        error("betas must be :all, :none, or :lead (got $betas)")

    # --------------------- Vectorize arguments if necessary --------------------- #

    # Vectorize genotype
    geno = geno isa AbstractVector ? reshape(geno, :, 1) : geno

    # Retrieve variant names
    variant_names = geno isa Matrix ? ["variant_$(i)" for i in 1:size(geno, 2)] : names(geno, 2)

    # Vectorize term tested if only one is provided
    termtest = termtest isa Vector{String} ? termtest : [termtest]

     # Vectorize genotype
    groups = groups isa AbstractVector ? DataFrame(cluster = groups) : groups

    # Vectorize B
    B = B isa AbstractVector ? B : [B]

    # ------------- Validate dimensionality of input data structures ------------- #

    verify_nobs_map_locus(geno, pheno, meta, groups)

    # ---------------------------- Add phenotype data ---------------------------- #

    design = deepcopy(meta)

    # -------------------- Set R and r for hypothesis testing -------------------- #

    # Define R matrix
    terms = termnames(f)[2]
    R = falses(length(termtest), length(terms))
    i_terms = [first(findall(bt .== terms)) for bt in termtest]

    if length(i_terms) != length(termtest)
        throw("Some termtest variables are not found in the formula. Verify all bootstrap terms tested are included")
    else
        for i in 1:length(termtest)
            R[i, i_terms[i]] = true
        end
    end

    r = fill(H0, length(i_terms))

    # ----------------------- Add response term to formula ----------------------- #

    f_re = f.rhs
    f = FormulaTerm(term(:E), f_re)

    # --------------------------- Add groups to design --------------------------- #

    for col_name in names(groups)
        design[!, col_name] = CategoricalArray(groups[!, col_name])
    end
    

    # Recode grouping variables
    for i in 1:ncol(groups)
        groups[!, i]  = levelcode.(CategoricalArray(groups[:, i]))
    end

    groups = Matrix(groups)

    # ------------------------------- Add phenotype ------------------------------ #

    design[!, :E] = pheno

    # ------------------ Shared null fit across variants (score test) ------------ #

    # For main/total effect tests every tested term involves the genotype G,
    # so the restricted (null) model contains no genotype-dependent terms and
    # is *identical for every cis variant* -- fit it once here and reuse it
    # across the whole locus, instead of refitting it per variant. An
    # interaction test keeps G in its null, so it falls back to the
    # per-variant behaviour automatically.
    v0 = vec(any(R, dims = 1))
    f0 = FormulaTerm(f.lhs, f.rhs[.!v0])
    null_is_shared = !(:G in StatsModels.termvars(f0.rhs))
    m0shared = null_is_shared ? glm(f0, design, Poisson(), LogLink()) : nothing
    # Only the fitted values and coefficients of the shared null are needed
    # downstream -- extracted here so the (large) fitted model object never
    # has to be serialized to parallel workers.
    μ̂0    = m0shared === nothing ? nothing : fitted(m0shared)
    beta0 = m0shared === nothing ? nothing : coef(m0shared)

    # A per-locus workspace avoids all per-variant DataFrame/model-matrix
    # machinery: the model matrix X differs between variants only in its
    # genotype-dependent column(s), so it is built once and only those
    # columns are updated in place per variant (see build_locus_workspace).
    # With a shared null (main/total), scores and X'diag(μ̂)X are also only
    # column-patched; with G in the null (interaction), the null is refit per
    # variant -- but directly on the patched matrix with warm starts, instead
    # of through the formula interface. Under `parallel`, each worker builds
    # its own private workspace for a contiguous block of variants (see
    # map_chunk); serially, one workspace is built here. Automatically
    # skipped (nothing) when the formula's genotype terms are not plain
    # linear-in-G terms or the workspace fails its self-check.
    use_parallel = parallel && nworkers() > 1
    ws = (!use_parallel && size(geno, 2) > 0) ?
        build_locus_workspace(f, design, μ̂0, beta0, v0, geno[:, 1], groups) : nothing

    # ---------- Run association for each variant in input genotype data --------- #

    t0 = time()

    results = if use_parallel

        # Contiguous blocks, one per worker: the workspace (and its warm-start
        # chain) lives on the worker, and the closure below -- including the
        # large `design`/`geno` captures -- is serialized once per worker via
        # the CachingPool instead of once per variant.
        nv = size(geno, 2)
        nchunks = min(nworkers(), nv)
        bounds = round.(Int, range(0, nv; length = nchunks + 1))
        chunks = [bounds[c]+1:bounds[c+1] for c in 1:nchunks if bounds[c+1] > bounds[c]]

        chunk_results = @showprogress pmap(CachingPool(workers()), chunks) do idxs
            map_chunk(idxs; f = f, design = design, μ̂0 = μ̂0, beta0 = beta0, v0 = v0,
                      groups = groups, geno = geno, R = R, r = r, boot = boot, B = B,
                      ptype = ptype, rboot = rboot, betas = betamode == :all)
        end
        reduce(vcat, chunk_results)

    elseif ws !== nothing
        @showprogress map(1:size(geno, 2)) do i

            safe_map_variant_shared(ws, geno[:, i]; groups = groups,
                    R = R, r = r, boot = boot, B = B, ptype = ptype, rboot = rboot,
                    rng = StableRNG(1322), betas = betamode == :all)
        end

    else
        @showprogress map(1:size(geno, 2)) do i

            safe_map_variant(geno[:, i]; f = f,  d = design, groups = groups, R = R, r = r,
                    boot = boot, B = B, ptype = ptype, rboot = rboot, rng = StableRNG(1322),
                    μ̂0 = μ̂0, beta0 = beta0, betas = betamode == :all)
        end

    end
    t1 = time()
    timewait = t1 - t0

    # ------------------------ Collect summary statistics ------------------------ #

    failed_variants = isnothing.(results)
    kept_idx = collect(1:length(results))  # geno column index per retained result row

    if any(failed_variants)
        
        failed_variants_names = variant_names[failed_variants]
        println(Crayon(foreground = :yellow), "The following variants failed:\n$(join(failed_variants_names, '\n'))")
        println(Crayon(foreground = :red), "Removing failed variants from output...")
        print(Crayon(reset = true))  # don't leak the red styling into subsequent output
        results = results[Not(failed_variants)]
        variant_names = variant_names[Not(failed_variants)]
        kept_idx = kept_idx[Not(failed_variants)]

    end


    # reduce(vcat, ...) hits DataFrames' specialized fast path; splatting
    # (vcat(results...)) recompiles for every result count and is slow for
    # loci with many variants.
    summ_stats = rboot ? reduce(vcat, [first(x) for x in results]) : reduce(vcat, results)
    insertcols!(summ_stats, 1, :variant => variant_names)

    # --------------------------- Lead-variant betas ----------------------------- #

    # betamode == :lead: effect sizes are usually only inspected for the top
    # association(s), so instead of fitting the unrestricted model for every
    # variant (:all), fit it only for the lead variant(s) -- smallest
    # analytical p-value, including all exactly tied variants -- and attach
    # their coefficient estimates, leaving every other row `missing`.
    if betamode == :lead && nrow(summ_stats) > 0
        pmin = minimum(summ_stats.p)
        lead_rows = findall(==(pmin), summ_stats.p)
        lead_fits = map(lead_rows) do i
            design[!, :G] = Float64.(geno[:, kept_idx[i]])
            m = glm(f, design, Poisson(), LogLink())
            (coefnames(m), coef(m))
        end
        cnames = lead_fits[1][1]
        for (ci, cn) in enumerate(cnames)
            col = Vector{Union{Missing, Float64}}(missing, nrow(summ_stats))
            for (l, i) in enumerate(lead_rows)
                col[i] = lead_fits[l][2][ci]
            end
            summ_stats[!, cn] = col
        end
    end
    
    # ---------------- Collect bootstrap distribution if necessary --------------- #

    boot_dist = rboot ? [last(x) for x in results] : []


    # --------------------------- Create Dynema object --------------------------- #
    res = DynemaModel(f, termtest, nrow(design), length(unique(groups[:, 1])), summ_stats,
                        B, boot_dist, timewait, true, boot, names(summ_stats)[2], pos, gene, chr)


    return(res)


end


function map_variant(variant::AbstractVector; f::FormulaTerm, d::AbstractDataFrame,
                    groups::Matrix, R::BitMatrix, r::Vector{Float64},
                    boot::Bool = true,
                    B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000],
                    ptype::Symbol = :equaltail, rboot = true, rng::AbstractRNG = StableRNG(66),
                    μ̂0 = nothing, beta0 = nothing, betas::Bool = true)

        # ------------- Add expression and genotype data to model matrix ------------- #

        design = deepcopy(d)
        design[!, :G] = Float64.(variant)

        # --------------------- Fit unrestricted model if needed ---------------------- #

        # The score test never needs the unrestricted fit; it is fit purely
        # to report per-variant coefficient estimates, so `betas = false`
        # skips it entirely.
        m = betas ? glm(f, design, Poisson(), LogLink()) : nothing

        # ---------------------- Extract predictors and response --------------------- #

        if m === nothing
            mf = ModelFrame(f, design)
            X = ModelMatrix(mf).m
            y = response(mf)
        else
            X = modelmatrix(m)
            y = response(m)
        end

        # ------------------------------ Fit null model ------------------------------- #

        # Impose null on model: reuse the locus-shared restricted fit's
        # fitted values/coefficients when provided (main/total effect tests
        # -- the null has no genotype-dependent terms, so it is identical for
        # every variant); otherwise fit it for this variant.
        v = vec(any(R, dims=1))
        μ̂, betas0 = if μ̂0 === nothing
            f0 = FormulaTerm(f.lhs, f.rhs[.!vec(v)])
            m0 = glm(f0, design, Poisson(), LogLink())
            fitted(m0), coef(m0)
        else
            μ̂0, beta0
        end

        # Calculate betas and variance-covariance matrix of the unrestricted model
        # at the solution of the restricted model.
        # Scores will be calculated the same way too
        betas_vec = insert_zeros(vec(v), betas0)

        A = X' * (Diagonal(μ̂) * X) \ I

        # ------------------------------- Build scores ------------------------------- #

        scores = (y .- μ̂) .* X

        # -------------------------- Calculate CRVE p-value -------------------------- #

        p_analytical_res = crvetest(R, r; resp = y, scores = scores, betas = betas_vec,
                    A = A, clustid = groups)

        res = DataFrame(p_analytical_res.stattype => p_analytical_res.stat)
        res[!, :p] = [p_analytical_res.p]

        # --------------------------------- Bootstrap -------------------------------- #

        if boot
            pboot, pbeta, bootdist = scorebootstrap(R, r; resp = y, scores = scores, betas = betas_vec,
                        A = A, clustid = groups,
                        rng = rng, ptype = ptype, B = B)
            res[!, :p_boot] = [pboot]
            res[!, :p_boot_approx] = [pbeta]
        end


        # ------------------- Extract betas from unrestricted model ------------------ #

        if betas
            res = hcat(res, DataFrame(transpose(coef(m)), coefnames(m)))
        end

        # ---------------- Return bootstrap distributions if required ---------------- #

        if rboot & boot
           res = (res, boot, bootdist)
        end

        return res

end


function safe_map_variant(variant; kwargs...)
    try
        map_variant(variant; kwargs...)
    catch err
        bt = catch_backtrace()
        @warn "map_variant failed for a variant" exception = (err, bt)
        nothing
    end
end


# ---------------------------------------------------------------------------- #
#                  Locus fast path (per-locus workspace)                       #
# ---------------------------------------------------------------------------- #
#
# The model matrix X differs between variants only in its genotype-dependent
# column(s), and since each such column is linear in G (`g .* multiplier`),
# the multipliers can be captured once by building X with G = 1. All
# per-variant DataFrame copies and model-matrix rebuilds are thereby
# eliminated; what remains depends on where G sits:
#
#   - Shared null (main/total effect tests): μ̂ comes from the single
#     locus-wide restricted fit, so residuals, the fixed columns of the
#     scores, and the fixed block of H = X'diag(μ̂)X are constant -- per
#     variant only q columns are patched, one O(n·k·q) product and a k×k
#     solve run, with zero large allocations.
#   - Per-variant null (interaction tests -- G stays in the null): the
#     restricted model is refit for each variant, but directly on its
#     patched model matrix with warm-started IRLS (previous variant's
#     solution), and μ̂/scores/H are recomputed into preallocated buffers.

"""
    LocusWorkspace

Preallocated per-locus buffers for the fast path. Built once by
[`build_locus_workspace`](@ref); the genotype-dependent columns of `X` are
overwritten in place for each variant by `map_variant_shared`.

Two modes, distinguished by `shared`:

  - `shared = true` (main/total effect tests): the null model is fit once for
    the locus, so `μ̂`, `resid`, `betas0`, the fixed columns of `scores`/`Sg`,
    and the fixed block of `H` are all constant -- only genotype
    columns/rows are patched per variant.
  - `shared = false` (interaction tests -- `G` remains in the null): the null
    model is refit per variant, directly on the patched null model matrix
    `X0` with IRLS warm-started from the previous variant's solution; `μ̂`,
    `resid`, `scores`, `betas0`, `H`, and `Sg` are then per-variant scratch
    buffers, recomputed in place.
"""
struct LocusWorkspace
    X::Matrix{Float64}       # n × k model matrix; genotype columns patched per variant
    scores::Matrix{Float64}  # n × k score matrix (resid .* X)
    M::Matrix{Float64}       # n × q multipliers of G for each genotype column
    gcols::Vector{Int}       # indices of the genotype-dependent columns of X
    H::Matrix{Float64}       # k × k X'diag(μ̂)X
    tmp::Matrix{Float64}     # n × q buffer for μ̂ .* X[:, gcols] (shared mode)
    Hq::Matrix{Float64}      # q × k buffer for the updated rows of H (shared mode)
    y::Vector{Float64}       # response
    μ̂::Vector{Float64}       # fitted values of the null model
    resid::Vector{Float64}   # y .- μ̂
    betas0::Vector{Float64}  # null coefficients, zero-padded to full length
    coefnames::Vector{String} # column names of X, for reporting unrestricted betas
    warm::Vector{Float64}    # warm-start coefficients for the unrestricted fit
    warmready::Base.RefValue{Bool} # has `warm` been populated by a fit yet?
    # -- per-variant null fit (shared = false only) --
    shared::Bool               # is the null model shared across variants?
    nullcols::Vector{Int}      # X columns belonging to the restricted (null) model
    X0::Matrix{Float64}        # n × k0 null model matrix; fixed columns constant
    g0full::Vector{Int}        # genotype columns of X that are also in the null...
    g0pos::Vector{Int}         # ...and their positions in X0/nullcols
    warm0::Vector{Float64}     # warm-start coefficients for the null fit
    warm0ready::Base.RefValue{Bool}
    Wbuf::Matrix{Float64}      # n × k buffer for μ̂ .* X (test H and IRLS Hessians)
    # -- custom warm-started IRLS (see irls_poisson!) --
    fitη::Vector{Float64}      # n: linear predictor scratch
    fitμ::Vector{Float64}      # n: mean scratch
    fitr::Vector{Float64}      # n: y - μ scratch
    Hfit::Matrix{Float64}      # k × k Hessian scratch for the unrestricted fit
    H0fit::Matrix{Float64}     # k0 × k0 Hessian scratch for the null fit
    irls_betas::Base.RefValue{Bool}         # use custom IRLS for the unrestricted fit?
    irls_betas_checked::Base.RefValue{Bool} # one-time self-check done?
    irls_null::Base.RefValue{Bool}          # use custom IRLS for the null fit?
    irls_null_checked::Base.RefValue{Bool}  # one-time self-check done?
    # -- direct CRVE score test (one-way clustering only) --
    clustid::Vector{Int}       # cluster code (1..G) per observation; empty if multiway
    clustshare::Vector{Float64} # each cluster's share of observations
    Sg::Matrix{Float64}        # G × k per-cluster score sums
    direct::Base.RefValue{Bool}  # use crvetest_direct instead of the library call?
    checked::Base.RefValue{Bool} # has the direct test been self-checked yet?
end

_contains_G(t) = :G in StatsModels.termvars(t)

# A term is "linear in G" when it is the bare `G` term or an interaction in
# which `G` appears exactly once as a bare component -- i.e. its
# model-matrix column(s) equal `g .* M` for a genotype-independent
# multiplier M. Anything else involving G (e.g. a FunctionTerm like log(G))
# disqualifies the fast path.
_linear_in_G(t::Term) = t.sym == :G
_linear_in_G(t::InteractionTerm) =
    count(_contains_G, t.terms) == 1 &&
    all(c -> !_contains_G(c) || (c isa Term && c.sym == :G), t.terms)
_linear_in_G(t) = false

"""
    build_locus_workspace(f, design, μ̂0, beta0, v, g1, groups) -> Union{Nothing, LocusWorkspace}

Builds the [`LocusWorkspace`](@ref) for the fast path -- in shared-null mode
when `μ̂0`/`beta0` are the locus-shared restricted fit's fitted values and
coefficients, or in per-variant-null mode when they are `nothing`
(interaction tests). Returns `nothing` (falling back to the generic
per-variant path) when any genotype term in `f` is not linear in `G`, when
construction fails, or when the patched model matrix does not reproduce a
from-scratch build for the first variant `g1` (numeric self-check).
Adds/overwrites a `:G` column on `design`.
"""
function build_locus_workspace(f::FormulaTerm, design::AbstractDataFrame, μ̂0, beta0,
                               v::AbstractVector{Bool}, g1::AbstractVector, groups::Matrix)

    rhs_terms = f.rhs isa Tuple ? collect(f.rhs) : [f.rhs]
    all(t -> !_contains_G(t) || _linear_in_G(t), rhs_terms) || return nothing

    ws = try

        # Build X once with G = 1: the genotype-dependent columns then hold
        # exactly their G-multipliers (1 for `G` itself, the context values
        # for a `G & C` column), so thereafter column(g) = g .* multiplier.
        design[!, :G] = ones(Float64, nrow(design))
        mf = ModelFrame(f, design)
        mm = ModelMatrix(mf)
        X = mm.m
        y = Float64.(response(mf))

        gterm_idx = findall(_contains_G, rhs_terms)
        gcols = findall(a -> a in gterm_idx, mm.assign)
        M = X[:, gcols]

        n, k, q = size(X, 1), size(X, 2), length(gcols)
        shared = μ̂0 !== nothing

        if shared
            μ̂ = μ̂0
            resid = y .- μ̂
            scores = resid .* X
            # Same expression as the generic path, so the fixed (covariate)
            # block of H is numerically identical to a per-variant computation.
            H = X' * (Diagonal(μ̂) * X)
            betas0 = insert_zeros(v, beta0)
            warm = copy(betas0)
            nullcols = Int[]
            X0 = Matrix{Float64}(undef, 0, 0)
            g0full = Int[]; g0pos = Int[]
            warm0 = Float64[]
        else
            # Per-variant null: allocate scratch buffers; the null model
            # matrix X0 is the subset of X's columns whose terms are not
            # tested (columns constant across variants, except G's).
            μ̂ = zeros(Float64, n)
            resid = zeros(Float64, n)
            scores = zeros(Float64, n, k)
            H = zeros(Float64, k, k)
            betas0 = zeros(Float64, k)
            warm = zeros(Float64, k)
            tested_terms = findall(v)
            nullcols = findall(a -> !(a in tested_terms), mm.assign)
            X0 = X[:, nullcols]
            g0full = [j for j in gcols if j in nullcols]
            g0pos = Int[findfirst(==(j), nullcols) for j in g0full]
            warm0 = zeros(Float64, length(nullcols))
        end
        Wbuf = Matrix{Float64}(undef, n, k)

        # Cluster structures for the direct CRVE score test. Only one-way
        # clustering is implemented directly; with multiway clustering the
        # fast path keeps using the WildBootTests-based crvetest instead.
        onecluster = size(groups, 2) == 1
        clustid = onecluster ? Int.(vec(groups[:, 1])) : Int[]
        G = onecluster ? maximum(clustid) : 0
        Sg = zeros(Float64, G, k)
        clustshare = zeros(Float64, G)
        if onecluster
            @inbounds for i in eachindex(clustid)
                clustshare[clustid[i]] += 1.0
            end
            clustshare ./= n
            if shared
                @inbounds for j in 1:k, i in eachindex(clustid)
                    Sg[clustid[i], j] += scores[i, j]
                end
            end
        end

        LocusWorkspace(X, scores, M, gcols, H,
                       Matrix{Float64}(undef, n, q), Matrix{Float64}(undef, q, k),
                       y, μ̂, resid, betas0, String.(coefnames(mf)), warm, Ref(shared),
                       shared, nullcols, X0, g0full, g0pos, warm0, Ref(false), Wbuf,
                       Vector{Float64}(undef, n), Vector{Float64}(undef, n),
                       Vector{Float64}(undef, n), Matrix{Float64}(undef, k, k),
                       Matrix{Float64}(undef, length(nullcols), length(nullcols)),
                       Ref(true), Ref(false), Ref(true), Ref(false),
                       clustid, clustshare, Sg, Ref(onecluster), Ref(false))

    catch err
        @warn "Could not set up the locus fast path; falling back to per-variant model matrices" exception = err
        nothing
    end
    ws === nothing && return nothing

    # Numeric self-check: patching the genotype columns with the first
    # variant must reproduce a from-scratch model-matrix build. Guards the
    # term-ordering/linearity assumptions above against StatsModels changes.
    for (l, j) in enumerate(ws.gcols)
        @views ws.X[:, j] .= g1 .* ws.M[:, l]
    end
    design.G .= g1
    X_ref = ModelMatrix(ModelFrame(f, design)).m
    if !(size(X_ref) == size(ws.X) && isapprox(X_ref, ws.X; rtol = 1e-10))
        @warn "Locus fast path failed its model-matrix self-check; falling back to per-variant model matrices"
        return nothing
    end

    return ws

end

# Fused mean update and Poisson/log-link deviance: μ .= exp.(η), returning
# 2 Σᵢ [yᵢ log(yᵢ/μᵢ) - (yᵢ - μᵢ)] (the yᵢ = 0 terms reduce to μᵢ).
# Overflow/underflow in exp propagates to a non-finite deviance, which the
# caller treats as a rejected step.
function _poisson_dev!(μ::Vector{Float64}, η::Vector{Float64}, y::Vector{Float64})
    dev = 0.0
    @inbounds for i in eachindex(η)
        m = exp(η[i])
        μ[i] = m
        yi = y[i]
        dev += yi > 0 ? yi * log(yi / m) - (yi - m) : m
    end
    2 * dev
end

"""
    irls_poisson!(β, X, y, η, μ, r, W, Hf; maxiter = 30, tol = 1e-9) -> Bool

In-place Newton/IRLS for a Poisson GLM with log link (canonical, so Fisher
scoring = Newton), warm-started from the incoming `β` and using only the
caller's preallocated buffers: `η`, `μ`, `r` (length n) and `W` (n × k),
`Hf` (k × k) scratch. Each iteration solves `X'diag(μ)X δ = X'(y - μ)` by
Cholesky and step-halves on the deviance. Returns `true` on convergence
(relative deviance change below `tol`), leaving the estimates in `β` and the
matching means in `μ`; returns `false` (caller falls back to `GLM.fit`) on
non-convergence, a failed Cholesky, or a collapsed line search.
"""
function irls_poisson!(β::Vector{Float64}, X::AbstractMatrix{Float64}, y::Vector{Float64},
                       η::Vector{Float64}, μ::Vector{Float64}, r::Vector{Float64},
                       W::AbstractMatrix{Float64}, Hf::AbstractMatrix{Float64};
                       maxiter::Int = 30, tol::Float64 = 1e-9, minstepfac::Float64 = 1e-3)

    k = size(X, 2)
    g = Vector{Float64}(undef, k)
    δ = Vector{Float64}(undef, k)
    βtry = Vector{Float64}(undef, k)

    mul!(η, X, β)
    dev = _poisson_dev!(μ, η, y)
    isfinite(dev) || return false

    for _ in 1:maxiter

        r .= y .- μ
        mul!(g, transpose(X), r)
        W .= μ .* X
        mul!(Hf, transpose(X), W)
        ch = cholesky!(Symmetric(Hf); check = false)
        issuccess(ch) || return false
        ldiv!(δ, ch, g)

        # Step-halving on the deviance (accept non-increase up to roundoff)
        devold = dev
        f = 1.0
        while true
            @. βtry = β + f * δ
            mul!(η, X, βtry)
            dev = _poisson_dev!(μ, η, y)
            (isfinite(dev) && dev <= devold * (1 + 1e-12) + 1e-12) && break
            f /= 2
            f < minstepfac && return false
        end
        β .= βtry

        abs(devold - dev) < tol * (abs(dev) + 0.1) && return true

    end

    return false

end

"""
    map_variant_shared(ws, variant; ...)

Fast-path equivalent of [`map_variant`](@ref): updates the genotype-dependent
columns of `ws` in place and runs the same CRVE score test (and optional
bootstrap) as the generic path. With a shared null only those columns are
touched; with a per-variant null (interaction tests) the restricted model is
refit on the patched null matrix `ws.X0` with warm-started IRLS. When
`betas = true`, the unrestricted model is also fit per variant to report
coefficient estimates -- directly on `ws.X`/`ws.y`, skipping the per-variant
DataFrame -> schema -> model-matrix rebuild that the formula interface would
repeat for every variant.
"""
function map_variant_shared(ws::LocusWorkspace, variant::AbstractVector;
                    groups::Matrix, R::BitMatrix, r::Vector{Float64},
                    boot::Bool, B::Vector{Int64}, ptype::Symbol, rboot::Bool,
                    rng::AbstractRNG, betas::Bool)

        # ---------------- Patch genotype columns of X in place ----------------------- #

        for (l, j) in enumerate(ws.gcols)
            @views ws.X[:, j] .= variant .* ws.M[:, l]
        end

        if ws.shared

            # Shared null: only the genotype columns of scores/H change.
            for (l, j) in enumerate(ws.gcols)
                @views ws.scores[:, j] .= ws.resid .* ws.X[:, j]
                @views ws.tmp[:, l] .= ws.μ̂ .* ws.X[:, j]
            end

            # Update the genotype rows/columns of H = X'diag(μ̂)X; the fixed
            # (covariate) block was computed once for the locus.
            mul!(ws.Hq, transpose(ws.tmp), ws.X)
            for (l, j) in enumerate(ws.gcols)
                @views ws.H[j, :] .= ws.Hq[l, :]
                @views ws.H[:, j] .= ws.Hq[l, :]
            end

        else

            # Per-variant null (interaction tests): refit the restricted
            # model directly on its patched model matrix via the custom
            # warm-started IRLS (previous variant's solution). The first
            # variant is fit cold with GLM.fit -- which also seeds the warm
            # start -- and the first warm IRLS solution is self-checked
            # against a cold GLM.fit once; any IRLS failure falls back to
            # GLM.fit for that variant.
            for (jf, j0) in zip(ws.g0full, ws.g0pos)
                @views ws.X0[:, j0] .= ws.X[:, jf]
            end
            k0 = length(ws.nullcols)
            W0 = view(ws.Wbuf, :, 1:k0)
            β0hat = copy(ws.warm0)
            ok = ws.warm0ready[] && ws.irls_null[] &&
                 irls_poisson!(β0hat, ws.X0, ws.y, ws.fitη, ws.fitμ, ws.fitr, W0, ws.H0fit)

            if ok && !ws.irls_null_checked[]
                ws.irls_null_checked[] = true
                β0ref = coef(GLM.fit(GeneralizedLinearModel, ws.X0, ws.y, Poisson(), LogLink()))
                if !isapprox(β0hat, β0ref; rtol = 1e-4, atol = 1e-8)
                    @warn "Custom IRLS failed its self-check against GLM.fit for the null fit; using GLM.fit"
                    ws.irls_null[] = false
                    ok = false
                end
            end

            if ok
                ws.μ̂ .= ws.fitμ
            else
                m0 = GLM.fit(GeneralizedLinearModel, ws.X0, ws.y, Poisson(), LogLink())
                β0hat = coef(m0)
                copyto!(ws.μ̂, fitted(m0))
            end
            copyto!(ws.warm0, β0hat)
            ws.warm0ready[] = true

            ws.resid .= ws.y .- ws.μ̂
            ws.scores .= ws.resid .* ws.X
            fill!(ws.betas0, 0.0)
            ws.betas0[ws.nullcols] .= β0hat

            # Full recompute of H = X'diag(μ̂)X: μ̂ changed, so every block did.
            ws.Wbuf .= ws.μ̂ .* ws.X
            mul!(ws.H, transpose(ws.X), ws.Wbuf)

        end

        A = ws.H \ I

        # -------------------------- Calculate CRVE p-value -------------------------- #

        p_analytical_res = if ws.direct[]
            # Per-cluster score sums: with a shared null only the genotype
            # columns change (fixed columns were precomputed at build time);
            # with a per-variant null the residuals changed, so all columns
            # are re-accumulated.
            updatecols = ws.shared ? ws.gcols : collect(1:size(ws.Sg, 2))
            @inbounds for j in updatecols
                Sgj = view(ws.Sg, :, j)
                fill!(Sgj, 0.0)
                scj = view(ws.scores, :, j)
                for i in eachindex(ws.clustid)
                    Sgj[ws.clustid[i]] += scj[i]
                end
            end
            direct = crvetest_direct(R, A, ws.Sg, ws.clustshare)
            if !ws.checked[]
                # One-time self-check of the direct statistic/p-value against
                # the WildBootTests implementation, on the first variant --
                # skipped when the optional WildBootTests package isn't
                # installed (the direct implementation is then trusted
                # as-is).
                ws.checked[] = true
                if HAS_WILDBOOTTESTS[]
                    ref = crvetest(R, r; resp = ws.y, scores = ws.scores, betas = ws.betas0,
                                A = A, clustid = groups)
                    if !(isapprox(ref.stat, direct.stat; rtol = 1e-8) &&
                         isapprox(ref.p, direct.p; rtol = 1e-6))
                        @warn "Direct CRVE score test failed its self-check against WildBootTests " *
                              "(stat $(direct.stat) vs $(ref.stat)); falling back to the library implementation"
                        ws.direct[] = false
                        direct = (stat = ref.stat, p = ref.p, stattype = ref.stattype)
                    end
                end
            end
            direct
        else
            crvetest(R, r; resp = ws.y, scores = ws.scores, betas = ws.betas0,
                        A = A, clustid = groups)
        end

        res = DataFrame(p_analytical_res.stattype => p_analytical_res.stat)
        res[!, :p] = [p_analytical_res.p]

        # --------------------------------- Bootstrap -------------------------------- #

        if boot
            pboot, pbeta, bootdist = scorebootstrap(R, r; resp = ws.y, scores = ws.scores, betas = ws.betas0,
                        A = A, clustid = groups,
                        rng = rng, ptype = ptype, B = B)
            res[!, :p_boot] = [pboot]
            res[!, :p_boot_approx] = [pbeta]
        end

        # ------------------- Extract betas from unrestricted model ------------------ #

        if betas
            # Unrestricted fit on the already-patched model matrix via the
            # custom IRLS, warm-started from the previous variant's solution
            # (or the current null solution for the first variant) --
            # neighbouring cis variants are in LD, so the optimum barely
            # moves and 2-3 iterations typically suffice. The first custom
            # solution is self-checked once against GLM.fit; any IRLS failure
            # falls back to GLM.fit for that variant.
            βhat = copy(ws.warmready[] ? ws.warm : ws.betas0)
            ok = ws.irls_betas[] &&
                 irls_poisson!(βhat, ws.X, ws.y, ws.fitη, ws.fitμ, ws.fitr, ws.Wbuf, ws.Hfit)

            if ok && !ws.irls_betas_checked[]
                ws.irls_betas_checked[] = true
                βref = coef(GLM.fit(GeneralizedLinearModel, ws.X, ws.y, Poisson(), LogLink()))
                if !isapprox(βhat, βref; rtol = 1e-4, atol = 1e-8)
                    @warn "Custom IRLS failed its self-check against GLM.fit for the unrestricted fit; using GLM.fit"
                    ws.irls_betas[] = false
                    ok = false
                end
            end

            if !ok
                βhat = coef(GLM.fit(GeneralizedLinearModel, ws.X, ws.y, Poisson(), LogLink()))
            end
            copyto!(ws.warm, βhat)
            ws.warmready[] = true
            res = hcat(res, DataFrame(transpose(βhat), ws.coefnames))
        end

        # ---------------- Return bootstrap distributions if required ---------------- #

        if rboot & boot
           res = (res, boot, bootdist)
        end

        return res

end


function safe_map_variant_shared(ws, variant; kwargs...)
    try
        map_variant_shared(ws, variant; kwargs...)
    catch err
        bt = catch_backtrace()
        @warn "map_variant failed for a variant" exception = (err, bt)
        nothing
    end
end


"""
    map_chunk(idxs; ...)

Processes a contiguous block of variants on one worker process (or serially,
if called directly): builds a private [`LocusWorkspace`](@ref) for the block
and runs [`map_variant_shared`](@ref) per variant, falling back to the
generic [`map_variant`](@ref) path if the workspace cannot be built. A
worker's in-place workspace mutation is safe because each worker holds its
own deserialized copy and `pmap` never runs two tasks concurrently on one
worker; keeping blocks contiguous also keeps each worker's warm-start chain
deterministic for a fixed worker count.
"""
function map_chunk(idxs::AbstractUnitRange; f::FormulaTerm, design::AbstractDataFrame,
                   μ̂0, beta0, v0::AbstractVector{Bool}, groups::Matrix, geno,
                   R::BitMatrix, r::Vector{Float64}, boot::Bool, B::Vector{Int64},
                   ptype::Symbol, rboot::Bool, betas::Bool)

    ws = isempty(idxs) ? nothing :
        build_locus_workspace(f, design, μ̂0, beta0, v0, geno[:, first(idxs)], groups)

    map(idxs) do i
        if ws === nothing
            safe_map_variant(geno[:, i]; f = f, d = design, groups = groups, R = R, r = r,
                    boot = boot, B = B, ptype = ptype, rboot = rboot, rng = StableRNG(1322),
                    μ̂0 = μ̂0, beta0 = beta0, betas = betas)
        else
            safe_map_variant_shared(ws, geno[:, i]; groups = groups, R = R, r = r,
                    boot = boot, B = B, ptype = ptype, rboot = rboot,
                    rng = StableRNG(1322), betas = betas)
        end
    end

end
