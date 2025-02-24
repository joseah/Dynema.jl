
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



function sample_index(rng::AbstractRNG, n::Integer)

    index = collect(1:n)
    index = sample(rng, index, length(index); replace=true)
    return index

end


function resample(rng::AbstractRNG, data::AbstractDataFrame, type::Symbol; cluster::Symbol)

    if type == :obs

        boot_id = sample_index(rng, nrow(data))

    else

        # Extract cluster info for each observation
        cluster_values = data[:, cluster]
        # Determine unique cluster labels
        cluster_levels = unique(cluster_values)
        # Resample clusters
        resample_levels = sample(rng, cluster_levels, length(cluster_levels); replace=true)

        
        cluster_dict= Dict()

        # Subset observations per cluster
        for level in resample_levels

            idx = findall(cluster_values .== level)
            cluster_dict[level] = idx
        
        end


        if type == :cluster
            boot_id = reduce(vcat, values(cluster_dict))
        
        elseif type == :two_stage

            cluster_dict_resampled= Dict()
            for key in keys(cluster_dict)
                orig_idx = cluster_dict[key]

                i = sample_index(rng, length(orig_idx))
                idx = orig_idx[i]
                cluster_dict_resampled[key] = idx
            end

            boot_id = reduce(vcat, values(cluster_dict_resampled))

        end

    end


    return(boot_id)

end