

function sample_index(rng::AbstractRNG, n::Integer)

    index = collect(1:n)
    index = sample(rng, index, length(index); replace=true)
    return index

end


function resample(rng::AbstractRNG, data::AbstractDataFrame)
    
    boot_id = sample_index(rng, nrow(data))
    return(boot_id)

end

function fit_model(f::FormulaTerm, resample::AbstractDataFrame)

    mdl = LinearMixedModel(f, resample)
    mdl.optsum.ftol_rel = 1e-8
    boot_model = fit!(mdl)
    DataFrame([boot_model.βs])

end


function boot_model(rng::AbstractRNG, md::AbstractDataFrame, f::FormulaTerm, n::Integer)

    res_boot = map(1:n) do i 
    index = resample(rng, md)
        fit_model(f, md[index, :])
    end

    return res_boot

end


function count_nonzeros(x)
    minimum([sum(x .> 0), sum(x .< 0)])
end
