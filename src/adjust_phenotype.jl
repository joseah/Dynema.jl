
"""
    adjust_phenotype(pheno::AbstractDataFrame, covariates::AbstractDataFrame)

Adjust phenotype information for biological and technical covariates using a 
    Poisson Generalized Linear Model.


# Arguments

`pheno` Phenotype information for each cell.

`covariates` Biological and technical covariates to regress out from phenotye information.


# Return

Vector with adjusted phenotype information (deviance residuals).

"""


function adjust_phenotype(pheno::AbstractDataFrame, covariates::AbstractDataFrame)

    fe_names = names(covariates)
    pheno_covs = hcat(pheno, covariates)

    f = term(:C) ~ sum(term.(fe_names))

    fit = GLM.glm(f, pheno_covs, Poisson())
    y = pheno_covs.C
    μ = predict(fit)
    res = @. sign(y - μ) * sqrt(devresid(Poisson(), y, μ))
    
    return(res)

end