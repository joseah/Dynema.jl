
"""
    adjust_phenotype(f::FormulaTerm, pheno::AbstractVector, covariates::AbstractDataFrame)

Adjust phenotype information for biological and technical covariates using 
    ordinary least squares (OLS)


# Arguments

`f` Formula object describing residualization strategy

`pheno` Phenotype information for each cell.

`covariates` Biological and/or technical covariates to regress out from the phenotype information.


# Return

Vector with adjusted phenotype information (OLS residuals).

"""


function adjust_phenotype(f::FormulaTerm, pheno::AbstractVector, covariates::AbstractDataFrame)

    if(length(pheno) != nrow(covariates)) throw("Phenotype and covariates have different number of observarions") end

    schema_f = schema(f, covariates)
    f_res = apply_schema(f, schema_f)
    _, predexog = modelcols(f_res, covariates)

    return ols_residuals(pheno, predexog)

end



function ols_residuals(y::AbstractVector, X::Matrix{Float64})

    β = Vector{Float64}(undef, size(X, 2))
    r = Vector{Float64}(undef, size(X, 1))
    β .= X \ y
    mul!(r, X, β)
    @. r = y - r
    return r

end