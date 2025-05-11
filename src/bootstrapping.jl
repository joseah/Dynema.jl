
function resample(rng::AbstractRNG, data::AbstractDataFrame)
    
    boot_id = Vector{Int64}(undef, nrow(data))
    sample!(rng, 1:nrow(data), boot_id; replace=true)

    return(boot_id)

end

function fit_model(f::FormulaTerm, resample::AbstractDataFrame)

    mdl = LinearMixedModel(f, resample)
    mdl.optsum.ftol_rel = 1e-8
    fit!(mdl)
    return(DataFrame([mdl.βs]))

end


function boot_model(rng::AbstractRNG, md::AbstractDataFrame, f::FormulaTerm, n::Integer)

    res_boot = map(1:n) do i 
        index = resample(rng, md)
        @views fit_model(f, md[index, :])
    end

    return res_boot

end
