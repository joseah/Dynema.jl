
"""
    adjust_phenotype(pheno::AbstractDataFrame, covariates::AbstractDataFrame)

Adjust phenotype information for biological and technical covariates using a 
    Poisson Generalized Linear Model.


# Arguments

`pheno` Phenotype information for each cell.

`covariates` Biological and/or technical covariates to regress out from the phenotye information.


# Return

Vector with adjusted phenotype information (deviance residuals).

"""


function adjust_phenotype(pheno::Vector{Int64}, covariates::AbstractDataFrame)

    fe_names = names(covariates)
    covariates[!, :C] = pheno

    f = term(:C) ~ sum(term.(fe_names))

    fit = GLM.glm(f, covariates, Poisson())
    y = covariates.C
    μ = predict(fit)
    # Clamp devresid result to avoid negative values due to floating-point precision errors
    res = @. sign(y - μ) * sqrt(max(devresid(Poisson(), y, μ), 0)) 
    
    return(res)

end