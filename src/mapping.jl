function map_locus(f::FormulaTerm; pheno::Vector{Float64}, geno::AbstractDataFrame,  metadata::AbstractDataFrame, 
                    groups::AbstractDataFrame, term::String, B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    r = [0], ptype::Symbol = :equaltail, rboot = false, rng::AbstractRNG = MersenneTwister(66))
    
    # ------------- Validate dimensionality of input data structures ------------- #

    verify_nobs_map_locus(geno, pheno, metadata, groups)

    # ---------------------------- Add phenotype data ---------------------------- #

    design = deepcopy(metadata)

    # ------------ Run association for each SNP in input genotype data ----------- #

    results = @showprogress pmap(eachcol(geno)) do snp
        
        map_snp(snp; f = f, pheno = pheno, metadata = design, groups = groups, term = term,
                B = B, r = r, ptype = ptype, rboot = rboot, rng = rng)
    end

    # ------------------------ Collect summary statistics ------------------------ #

    summ_stats = rboot ? reduce(vcat, [first(x) for x in results]) : vcat(results...)
    insertcols!(summ_stats, 1, :snp => names(geno))
    
    # ---------------- Collect bootstrap distribution if necessary --------------- #

    boot_dist = rboot ? [last(x) for x in results] : []


    # --------------------------- Create Dynema object --------------------------- #
    res = DynemaModel(f, term, nrow(design), length(unique(groups[:, 1])), summ_stats, B, boot_dist)


    return(res)


end





#function map_snp(snp::AbstractVector; f::FormulaTerm;  pheno::Vector{Float64}, data::AbstractDataFrame,
#                    groups::AbstractDataFrame, term::String, B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
#                    r = [0], ptype::Symbol = :equaltail, dist::Bool = false, rboot = true, rng::AbstractRNG = MersenneTwister(66))

#function map_snp(snp::AbstractVector; kargs...)

function map_snp(snp::AbstractVector; f::FormulaTerm, pheno::Vector{Float64}, metadata::AbstractDataFrame,
                    groups::AbstractDataFrame, term::String, B::Vector{Int64} = [200, 200, 1600, 2000, 16000, 20000], 
                    r = [0], ptype::Symbol = :equaltail, rboot = true, rng::AbstractRNG = MersenneTwister(66))

        # Add expression and genotype data to covariates
        design = deepcopy(metadata)
        design[!, :E] = pheno
        design[!, :G] = Float64.(snp)

        # Model data based on formula
        f = apply_schema(f, schema(f, design))
        _ , predexog = modelcols(f, design)

        # Recode grouping variables
        for i in 1:ncol(groups)
            groups[!, i]  = levelcode.(CategoricalArray(groups[:, i]))
        end
        groups = Matrix(groups)    
        

        terms = termnames(f)[2]
        
        R = terms .== term

        if !any(R)
            throw("Term '$term' not found in formula. Verify the term in included in the formula")
        else
            R = reshape(R, 1, length(R))
        end

        
        # --------------------- Run first round of bootstrapping --------------------- #
        B_total = B[1] - 1
        test = wildboottest(R, r; resp = pheno, predexog = predexog, clustid = groups, 
                            rng = rng, ptype = ptype, reps = B_total)
        
        # Extract results
        stat = teststat(test)
        boot = dist(test)[1, :]
        counts = pass(stat, boot)
        #stat_type = stattype(test)
        b = test.b

        if counts <= 20

            for j in 1:length(B)
                
                test = wildboottest(R, r; resp = pheno, predexog = predexog, clustid = groups, 
                                    rng = rng, ptype = ptype, reps = B[j])
                boot_i = dist(test)[1, :]
                boot = vcat(boot, boot_i)
                total_B = length(boot)
                counts = pass(stat, boot)
                
                if counts > 20
                    break
                end
            
            end

        end

        # Compute final p-value 
        pval = compute_pvalue(stat, boot, B_total)
        
        # Gather results
        res = (DataFrame(coef = b, stat = stat, p = pval))

        if rboot
           res = (res, boot)
        end

        return res

end