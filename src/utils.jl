function printlnv(args...; verbose = true)
    if verbose
        println(args...)
    end
end


function verify_nobs_map_locus(geno::AbstractDataFrame, pheno::AbstractVector, 
                                meta::AbstractDataFrame, groups::AbstractDataFrame)

    geno_size = nrow(geno)
    pheno_size = length(pheno)
    meta_size = nrow(meta)
    groups_size = nrow(groups)



   n_obs = Dict(:snp => geno_size, :expr => pheno_size, 
                 :meta => meta_size, :groups => groups_size)

   n_obs_vals  = collect(values(n_obs))

   if(!all(x -> x == first(n_obs_vals), n_obs_vals))
       
       println("Verify number of observatios of input files")
       println("Number of observations per input data\n")
       for (key, value) in pairs(n_obs)
           println(key, " = ", value)
           
       end
       throw("Number of observations in input data do not match between input files")
       
   end


end