"""
    boot_model(rng, form, data, geno, geno_index, target_params, n)

Fits linear mixed model `n` times for a given SNP.

Parameters


"""


function boot_snp(rng::AbstractRNG, f::FormulaTerm, data::AbstractDataFrame, 
                  geno::AbstractDataFrame, snp_index::Integer, target_params::Vector{Symbol}, 
                  n::Integer)

    # Create copy of covariate and response data
    set = deepcopy(data)

    # Add genotype information for snp defined by the index "geno_index"
    # snp_index should start from 2, as column 1 from genotype data 
    # corresponds to the donor id
    set[!, :G] = geno[:, snp_index]

    # Fit `n_boot` models using sampling with replacement
    boot = boot_model(rng, set, f, n)

    # Gather all betas from all bootstrap fits
    boot = reduce(vcat, boot)
    
    # Return bootstrapping results for parameters of interest
    if typeof(target_params) == Symbol
        boot = DataFrame([boot[:, target_params]], [target_params])
    else
        boot = boot[:, target_params]
    end

    return boot

end