function map_locus(f::FormulaTerm; pheno::AbstractVector, geno::Union{AbstractDataFrame, AbstractVector}, meta::AbstractDataFrame, 
                    groups::Union{AbstractDataFrame, AbstractVector}, termtest::Union{String, Vector{String}}, H0::Float64 = Float64(0), 
                    imposenull::Bool = true,
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

    # ------------- Validate dimensionality of input data structures ------------- #

    verify_nobs_map_locus(geno, pheno, meta, groups)

    # ---------------------------- Add phenotype data ---------------------------- #

    design = deepcopy(meta)

    # -------------------- Set R and r for hypothesis testing -------------------- #

   

    # Define R matrix
    terms = termnames(f)[2]
    R = falses(1, length(terms))
    i_terms = [first(findall(bt .== terms)) for bt in termtest]

    if length(i_terms) != length(termtest)
        throw("Some termtest variables are not found in the formula. Verify all bootstrap terms tested are included")
    else
        R[1, i_terms] .= true
    end

    # r (by default zero)
    r = [H0]

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
    results = @showprogress pmap(eachcol(geno)) do snp
        
        map_snp(snp; f = f,  d = design, groups = groups, R = R, r = r, imposenull = imposenull, 
                B = B, ptype = ptype, rboot = rboot, rng = rng)
    end
    t1 = time()
    timewait = t1 - t0

    # ------------------------ Collect summary statistics ------------------------ #

    summ_stats = rboot ? reduce(vcat, [first(x) for x in results]) : vcat(results...)
    insertcols!(summ_stats, 1, :snp => names(geno))
    
    # ---------------- Collect bootstrap distribution if necessary --------------- #

    boot_dist = rboot ? [last(x) for x in results] : []


    # --------------------------- Create Dynema object --------------------------- #
    res = DynemaModel(f, termtest, nrow(design), length(unique(groups[:, 1])), summ_stats, B, boot_dist, timewait, pos, gene, chr)


    return(res)


end


function map_snp(snp::AbstractVector; f::FormulaTerm, d::AbstractDataFrame,
                    groups::Matrix, R::BitMatrix, r::Vector{Float64}, imposenull::Bool = true, B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    ptype::Symbol = :equaltail, rboot = true, rng::AbstractRNG = StableRNG(66))

        
        # ------------- Add expression and genotype data to model matrix ------------- #

        design = deepcopy(d)
        design[!, :G] = Float64.(snp)

        # ---------------------------------- Fit GLM --------------------------------- #
        
        m = glm(f, design, Poisson(), LogLink())

        # ------ Extract naive statistics (useful for diagnostics and debuggin) ------ #

        p_naive = coeftable(m).cols[4][vec(R)]
        p_naive = DataFrame(transpose(p_naive), "p_naive - " .* coefnames(m)[vec(R)])

        # ---------------------- Extract predictors and response --------------------- #

        X = modelmatrix(m)
        y = response(m)

        if imposenull
            
            # Impose null on model
            f0 = FormulaTerm(f.lhs, f.rhs[.!vec(R)])
            m0 = glm(f0, design, Poisson(), LogLink())

            # Calculate betas and variance-covariance matrix of the unrestricted model
            # at the solution of the restricted model.
            # Scores will be calculated the same way too
            betas0 = coef(m0)
            betas = insert_zeros(vec(R), betas0)

            μ̂ = fitted(m0)
            A = (X' * (μ̂ .* X)) \ I

        else

            betas = coef(m)
            μ̂ = fitted(m)
            A = vcov(m)

        end
        

        # ------------------------------- Build scores ------------------------------- #

        scores = (y .- μ̂) .* X

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
        counts = pass(statistic, boot)

        # ---------------- Remaining rounds of bootstrapping if needed --------------- #

        if counts <= 20

            for j in 2:length(B)
                
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
                            rng = rng, ptype = ptype, 
                            reps = B[j])
                
                
                boot_i = dist(test)[1, :]
                boot = vcat(boot, boot_i)
                counts = pass(statistic, boot)
                
                if counts > 20
                    break
                end
            
            end

        end

        # Compute final p-value 
        pval = compute_pvalue(statistic, boot)

    
        # -------------------------- Calculate CRVE p-value -------------------------- #

        vcov_cr = vcov(CR0(groups), m)
        zr = coef(m) ./ sqrt.(diag(vcov_cr))
        zr = zr[vec(R)]
        pval_CRVE = [2 .*(1 .- cdf.(Normal(), abs.(x))) for x in zr]
        pval_CRVE = DataFrame(transpose(pval_CRVE), "p_analytical - " .* coefnames(m)[vec(R)])


        # ------------------- Extract betas from unrestricted model ------------------ #

        betas = DataFrame(transpose(coef(m)), coefnames(m))

        # ------------------------------ Gather results ------------------------------ #
        
        res = DataFrame(stat = statistic)
        res = hcat(betas, res)
        res = hcat(res, pval_CRVE)
        res[!, :p_boot] = [pval]
        res = hcat(res, p_naive)

        # ---------------- Return bootstrap distributions if required ---------------- #

        if rboot
           res = (res, boot)
        end

        return res

end