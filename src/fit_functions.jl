
function fit_model(form, resample)

    mdl = LinearMixedModel(form, resample)
    mdl.optsum.ftol_rel = 1e-8
    boot_model = fit!(mdl)
    DataFrame([boot_model.βs])

end


function sample_index(rng, n_cells)

    index = collect(1:n_cells)
    index = sample(rng, index, length(index); replace=true)
    return index

end


function fit_lmm(rng, md, form, n; batch_size = 30)

    res_boot = pmap(1:n, batch_size = batch_size) do i 
    index = sample_index(rng, nrow(md))
        fit_model(form, md[index, :])
    end

    return res_boot

end

export fit_model
export sample_index
export fit_lmm