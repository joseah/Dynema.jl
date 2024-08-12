
function fit_model(f::FormulaTerm, resample::AbstractDataFrame)

    mdl = LinearMixedModel(f, resample)
    mdl.optsum.ftol_rel = 1e-8
    boot_model = fit!(mdl)
    DataFrame([boot_model.βs])

end



function boot_model(rng::AbstractRNG, md::AbstractDataFrame, f::FormulaTerm, n::Integer)

    res_boot = map(1:n) do i 
    index = sample_index(rng, nrow(md))
    fit_model(f, md[index, :])

    return res_boot

end