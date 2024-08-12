
function fit_model(f::FormulaTerm, resample::AbstractDataFrame)

    mdl = LinearMixedModel(f, resample)
    mdl.optsum.ftol_rel = 1e-8
    boot_model = fit!(mdl)
    DataFrame([boot_model.βs])

end



function boot_model(rng::AbstractRNG, md::AbstractDataFrame, f::FormulaTerm, n::Integer, type::Symbol; cluster::Symbol)

    res_boot = map(1:n) do i 
    index = resample(rng, md, type; cluster = cluster)
        fit_model(f, md[index, :])
    end

    return res_boot

end

