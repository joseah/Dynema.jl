function map_locus(f::FormulaTerm; pheno::AbstractVector, geno::AbstractDataFrame,  meta::AbstractDataFrame, 
                    groups::AbstractDataFrame, bterm::String, B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    r = [0],  ptype::Symbol = :equaltail, rboot = false, rng::AbstractRNG = MersenneTwister(66), 
                    pos::Union{Nothing, Vector{Int64}, Vector{Float64}} = nothing,
                    gene::Union{Nothing, String} = nothing,
                    chr::Union{Nothing, String, Int} = nothing)
    
    # ------------- Validate dimensionality of input data structures ------------- #

    verify_nobs_map_locus(geno, pheno, meta, groups)

    # ---------------------------- Add phenotype data ---------------------------- #

    design = deepcopy(meta)

    # ------------------------------ Verify formula ------------------------------ #

    terms = termnames(f)[2]
        
    R = terms .== bterm

    if !any(R)
        throw("Term '$bterm' not found in formula. Verify the bootstrap term in included in the formula")
    else
        R = reshape(R, 1, length(R))
    end

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
        
        map_snp(snp; f = f,  d = design, groups = groups, bterm = bterm,
                B = B, R = R, r = r, ptype = ptype, rboot = rboot, rng = rng)
    end
    t1 = time()
    timewait = t1 - t0

    # ------------------------ Collect summary statistics ------------------------ #

    summ_stats = rboot ? reduce(vcat, [first(x) for x in results]) : vcat(results...)
    insertcols!(summ_stats, 1, :snp => names(geno))
    
    # ---------------- Collect bootstrap distribution if necessary --------------- #

    boot_dist = rboot ? [last(x) for x in results] : []


    # --------------------------- Create Dynema object --------------------------- #
    res = DynemaModel(f, bterm, nrow(design), length(unique(groups[:, 1])), summ_stats, B, boot_dist, timewait, pos, gene, chr)


    return(res)


end


function map_snp(snp::AbstractVector; f::FormulaTerm, d::AbstractDataFrame,
                    groups::Matrix, bterm::String, B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    R::BitMatrix, r = [0],  ptype::Symbol = :equaltail, rboot = true, rng::AbstractRNG = MersenneTwister(66))

        
        # Add expression and genotype data to covariates
        design = deepcopy(d)
        design[!, :G] = Float64.(snp)

        # Fit GLM
        m = glm(f, design, Poisson(), LogLink())


        # ------------------------------- Build scores ------------------------------- #
        
        p_naive_analytical = last(coeftable(m).cols[4][vec(R)])
        X = modelmatrix(m)
        μ̂ = fitted(m)
        y = response(m)
        A = vcov(m)
        betas = coef(m)
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
                            imposenull = false,
                            rng = rng, ptype = ptype, reps = B[1] - 1)
        



        # Extract results
        stat = teststat(test)
        boot = dist(test)[1, :]
        counts = pass(stat, boot)
        b = test.b

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
                            imposenull = false,
                            rng = rng, ptype = ptype, 
                            reps = B[j])
                
                
                boot_i = dist(test)[1, :]
                boot = vcat(boot, boot_i)
                counts = pass(stat, boot)
                
                if counts > 20
                    break
                end
            
            end

        end

        # Compute final p-value 
        pval = compute_pvalue(stat, boot)

        # Calculate analytical cluster robust p-value
        pval_CRVE = 2 .*(1 .- cdf.(Normal(), abs.(stat)))

        # Gather results
        res = (DataFrame(coef = b, stat = stat, p = pval, p_analytical = pval_CRVE, p_naive_analytical = p_naive_analytical))

        if rboot
           res = (res, boot)
        end

        return res

end