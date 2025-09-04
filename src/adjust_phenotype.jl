
"""
    adjust_phenotype(f::FormulaTerm, pheno::AbstractVector, covariates::AbstractDataFrame)

Adjust phenotype information for biological and technical covariates using 
    a Poisson Generalized Linear Model


# Arguments

`f` Formula object describing residualization strategy

`pheno` Phenotype information for each cell.

`covariates` Biological and/or technical covariates to regress out from the phenotype information.


# Return

Vector with adjusted phenotype information (raw residuals).

"""


function adjust_phenotype(f::FormulaTerm, pheno::AbstractVector, covariates::AbstractDataFrame)

    if(length(pheno) != nrow(covariates)) throw("Phenotype and covariates have different number of observarions") end

    schema_f = schema(f, covariates)
    f_res = apply_schema(f, schema_f)
    _, predexog = modelcols(f_res, covariates)

    m = GLM.glm(predexog, pheno, Poisson())
    resid = pheno - predict(m)

    return resid

end