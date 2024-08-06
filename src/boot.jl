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
                  geno::AbstractDataFrame, snp_index::Integer, n::Integer)

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
    
    # Return bootstrapping results
    return boot

end