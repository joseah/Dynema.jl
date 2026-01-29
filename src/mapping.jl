"""

`map_locus([f::FormulaTerm; pheno::AbstractVector, geno::Union{AbstractDataFrame, AbstractVector}, 
            meta::AbstractDataFrame, groups::Union{AbstractDataFrame, AbstractVector}, termtest::Union{String, Vector{String}},
             <optional keyword arguments>) -> Dynema.Dynema_struct`

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

- `parallel`: Runs in parallel via `Distrbuted.jl` with `pmap`
- `H0`: Null hypothesis value. By default `0`
- `imposenull`: Logical inditicating whether to impose a null in the bootstrap data generating process (DGP). If true, a score test 
is applied, otherwise a wald test is used
- B: Number of bootstrap iterations to apply. By default 39999 iterations at apply to achieve a p-value of 5 x 10^-5 for a two-tail test. 
For adaptive bootstrapping, the number of iterations might be specified as a vector indicating the number of iterations
to perform in each step
- ptype: Type of bootstrap p-value to return (:equaltail, :symmetric)
- rboot: Whether to return the bootstrap distributions for each variant. Useful for direct assessment of bootstrap statistics and
internal debugging
- rng: Random number generator
- pos: A numeric value specifying a genomic location for each genetic variant. Stored in final output for convenience
- gene: Name of the gene being tested. Stored in final output for convenience
- chr: Chromosome position of gene being tested. Stored in final output for convenience.
"""




function map_locus(f::FormulaTerm; pheno::AbstractVector, geno::Union{AbstractDataFrame, AbstractVector}, 
                    meta::AbstractDataFrame, groups::Union{AbstractDataFrame, AbstractVector}, termtest::Union{String, Vector{String}}, 
                    parallel = false,
                    H0::Float64 = Float64(0), imposenull::Bool = true,
                    B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    ptype::Symbol = :equaltail, rboot = false, rng::AbstractRNG = StableRNG(66), 
                    pos::Union{Nothing, Vector{Int64}, Vector{Float64}} = nothing,
                    gene::Union{Nothing, String} = nothing,
                    chr::Union{Nothing, String, Int} = nothing)

    # --------------------- Vectorize arguments if necessary --------------------- #

    # Vectorize genotype
    geno = geno isa AbstractVector ? DataFrame(G = geno) : geno
    
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

    # ------------ Run association for each SNP in input genotype data ----------- #

    t0 = time()

    results = if parallel

        @showprogress pmap(eachcol(geno)) do snp

            safe_map_snp(snp; f = f,  d = design, groups = groups, R = R, r = r, imposenull = imposenull, 
                    B = B, ptype = ptype, rboot = rboot, rng = rng)

        end

    else
        @showprogress map(eachcol(geno)) do snp
            
            safe_map_snp(snp; f = f,  d = design, groups = groups, R = R, r = r, imposenull = imposenull, 
                    B = B, ptype = ptype, rboot = rboot, rng = rng)
        end

    end
    t1 = time()
    timewait = t1 - t0

    # ------------------------ Collect summary statistics ------------------------ #

    failed_snps = isnothing.(results)

    snp_names = if any(failed_snps)
        
        failed_snps_names = names(geno)[failed_snps]
        println(Crayon(foreground = :yellow), "The following SNPs failed:\n$(join(failed_snps_names, '\n'))")
        println(Crayon(foreground = :red), "Removing failed SNPs from output...")
        results = results[Not(failed_snps)]
        
        names(geno)[Not(failed_snps)]

    else

        names(geno)

    end


    summ_stats = rboot ? reduce(vcat, [first(x) for x in results]) : vcat(results...)
    insertcols!(summ_stats, 1, :snp => snp_names)
    
    # ---------------- Collect bootstrap distribution if necessary --------------- #

    boot_dist = rboot ? [last(x) for x in results] : []


    # --------------------------- Create Dynema object --------------------------- #
    res = DynemaModel(f, termtest, nrow(design), length(unique(groups[:, 1])), summ_stats, 
                        B, boot_dist, timewait, imposenull, pos, gene, chr)


    return(res)


end


function map_snp(snp::AbstractVector; f::FormulaTerm, d::AbstractDataFrame,
                    groups::Matrix, R::BitMatrix, r::Vector{Float64}, imposenull::Bool = true, B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    ptype::Symbol = :equaltail, rboot = true, rng::AbstractRNG = StableRNG(66))

        
        # ------------- Add expression and genotype data to model matrix ------------- #

        design = deepcopy(d)
        # snp = geno[:, 1]
        design[!, :G] = Float64.(snp)

        # ---------------------------------- Fit GLM --------------------------------- #
        
        m = glm(f, design, Poisson(), LogLink())

        # ---------------------- Extract predictors and response --------------------- #

        X = modelmatrix(m)
        y = response(m)

        if imposenull
            
            # Impose null on model
            v = vec(any(R, dims=1))
            f0 = FormulaTerm(f.lhs, f.rhs[.!vec(v)])
            m0 = glm(f0, design, Poisson(), LogLink())

            # Calculate betas and variance-covariance matrix of the unrestricted model
            # at the solution of the restricted model.
            # Scores will be calculated the same way too
            betas0 = coef(m0)
            betas = insert_zeros(vec(v), betas0)

            μ̂ = fitted(m0)
            A = (X' * (μ̂ .* X)) \ I

        else

            betas = coef(m)
            μ̂ = fitted(m)
            A = vcov(m)

        end
        

        # ------------------------------- Build scores ------------------------------- #

        scores = (y .- μ̂) .* X

        # ------------------------Calculate analytical p-value------------------------ #


        p_analytical = if imposenull

               scoretest(R, r; 
                            resp = y, 
                            scores = scores,
                            beta = betas,
                            A = A,
                            clustid = groups, 
                            ml = true,
                            scorebs = true,
                            small = false)
        
        else

                waldtest(R, r; 
                            resp = y, 
                            scores = scores,
                            beta = betas,
                            A = A,
                            clustid = groups, 
                            ml = true,
                            scorebs = true,
                            small = false)
        end

        # --------------------- Run first round of bootstrapping --------------------- #


        test = wildboottest(R, r; 
                            resp = y, 
                            scores = scores,
                            beta = betas,
                            A = A,
                            clustid = groups, 
                            ml = true,
                            scorebs = true,
                            imposenull = imposenull,
                            small = false,
                            rng = rng, ptype = ptype, reps = B[1] - 1)

        # Extract results
        statistic = teststat(test)
        boot = dist(test)[1, :]
        stattype = test.stattype
        counts = pass(statistic, boot, stattype)

        # ---------------- Remaining rounds of bootstrapping if needed --------------- #
        
        # Compute final p-value 
        pval = compute_pvalue(statistic, boot, stattype)

        # ------------------- Extract betas from unrestricted model ------------------ #

        betas = DataFrame(transpose(coef(m)), coefnames(m))

        # ------------------------------ Gather results ------------------------------ #
        
        res = DataFrame(stattype => statistic)
        res[!, :p] = [pval]
        res[!, :p_analytical] = [p_analytical.p]
        res = hcat(res, betas)

        # ---------------- Return bootstrap distributions if required ---------------- #

        if rboot
           res = (res, boot)
        end

        return res

end


function safe_map_snp(snp; kwargs...)
    try
        map_snp(snp; kwargs...)
    catch err
        bt = catch_backtrace()
        @warn "map_snp failed for a SNP" exception = (err, bt)
        nothing
    end
end
