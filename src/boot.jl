"""
    boot_snp(rng, f, data, geno, snp_index, n)

Fits linear mixed model `n` times for a given SNP using sampling with replacement
at the observation level.

Parameters:

- `rng``: Random number generator
- `f``: Formula
- `data``: Dataframe containing response and covariates
- `geno``: Dataframe containing genotype data:
    + First column: donor id
    + Remaining columns: Each SNPs with respective allele dosages

    Both `data` and `geno` should be coded at the observational level.

- `n`: Number of bootstrap replicates


"""


function boot_snp(rng::AbstractRNG, f::FormulaTerm, data::AbstractDataFrame, 
                  geno::AbstractDataFrame, snp_index::Integer, n::Integer, by_snp::Bool=true)

    # Create copy of covariate and response data
    set = deepcopy(data)

    # Add genotype information for snp defined by the index "geno_index"
    # snp_index should start from 2, as column 1 from genotype data 
    # corresponds to the donor id
    set[!, :G] = geno[:, snp_index]

    # Fit `n_boot` models using sampling with replacement
    boot = boot_model(rng, set, f, n, by_snp)

    # Gather all betas from all bootstrap fits
    boot = reduce(vcat, boot)
    
    # Return bootstrapping results
    return boot

end


"""
    boot_locus(rng, form, md, geno, snp_set, n_boot)

Performs locus-wide eQTL mapping for a set of SNPs.

- `rng``: Random number generator
- `f``: Formula
- `data``: Dataframe containing response and covariates
- `geno``: Dataframe containing genotype data:
    + First column: donor id
    + Remaining columns: Each SNPs with respective allele dosages

    Both `data` and `geno` should be coded at the observational level.
- `snp_set``: A boolean vector specifying which SNPs to test for
- `n`: Number of bootstrap replicates



"""

function boot_locus(rng::AbstractRNG, f::FormulaTerm, data::AbstractDataFrame, 
    geno::AbstractDataFrame, snp_set, n::Integer, by_snp::Bool=true)

    println("Bootstrap n = $n")

    # Create symbol for current bootstrap iteration
    boot_symbol = Symbol("n_" * string(n))


    if by_snp
        # Test each SNP specified in `snp_set`
        boot_res = @showprogress pmap(snp_set) do snp_index
            boot = boot_snp(rng, f, data, geno, snp_index + 1, n, by_snp = true)
        end

    else
        # Test each SNP specified in `snp_set`
        boot_res = @showprogress map(snp_set) do snp_index
            boot = boot_snp(rng, f, data, geno, snp_index + 1, n, by_snp = false)
        end

    end

    # Gather results and return
    snp_names = names(geno)[snp_set .+ 1]
    boot_res = DataFrame(snp = snp_names, boot_n = boot_res)
    rename!(boot_res, :boot_n => boot_symbol)

    return boot_res

end


"""
    pass_boot(rng, form, md, geno, snp_set, n_boot)


- `rng``: Random number generator
- `f``: Formula
- `data``: Dataframe containing response and covariates
- `geno``: Dataframe containing genotype data:
    + First column: donor id
    + Remaining columns: Each SNPs with respective allele dosages

    Both `data` and `geno` should be coded at the observational level.
- `snp_set``: A boolean vector specifying which SNPs to test for
- `n`: Number of bootstrap replicates



"""

function pass_boot(res::AbstractDataFrame, n::Integer, boot_sizes::Vector{Int64}, target_params::Vector{Symbol})
    
    # Create bootstrap iteration symbols
    i = n .>=  boot_sizes
    boot_cumm_symbols  = @. Symbol("n_" * string(boot_sizes[i]))
    
    # For each SNP, gather all bootstrapping result across all
    # iterations. Count number of non-zero values for each 
    # parameter of interest and return a logical vector
    # specifying which parameters there's evidence of non-zero
    # values for.

    pass = map(eachrow(res)) do snp
        # Collect result across bootstrap iterations
        boot = reduce(vcat, snp[boot_cumm_symbols])

        # If only one parameter is of interest, create array
        # of patameters of length one
        if typeof(target_params) == Symbol
            target_params = [target_params]
        end

        # Subset fixed effect parameters of interest
        boot = boot[:, target_params]

        # Check if there are non-zero values for each parameter
        param_status = map(count_nonzeros, eachcol(boot))

        # Return logical vector specifying which parameters pass
        # condition
        any(param_status .== 0)
    end

    return pass

end




function map_locus(f::FormulaTerm, data::AbstractDataFrame, 
                    geno::AbstractDataFrame, n::Vector{Int64}, target_params::Vector{Symbol})

    # Create indivator vector to specify which SNPs to test for
    i = ncol(geno) - 1
    pass = 1:i

    # Extract gene names ancd create Dataframe to store all results
    snp_names = names(geno)[2:end]
    res = DataFrame(snp = snp_names)

    # For each bootstrap round, bootstrap model across all specified SNPs.
    # The list of SNPs will change depending on whether they keep being 
    # significant after each bootstrap round.

    for n_i in n
        current_boot = Symbol("n_" * string(n_i))
        # Performs locus-wide eQTL mapping for the specified set of SNPs.
        if n_i < 10_000
            res_boot = boot_locus(MersenneTwister(n_i), f, data, geno, pass, n_i, by_snp = trye)
        else
            res_boot = boot_locus(MersenneTwister(n_i), f, data, geno, pass, n_i, by_snp = false)
        end
        # Add current results to global results
        res = leftjoin(res, res_boot, on = :snp, order = :left)
        # Substitute missing values for non-significant SNPs with emoty dataframes 
        res[!, current_boot] = replace(res[:, current_boot], missing => DataFrame())
        # Test which SNPs will be tested for in the next bootstrapping iteration
        pass = findall(pass_boot(res, n_i, n, target_params))
        println("N. SNPs pass = ", length(pass))
        println("-------------")
    end

    return res

end