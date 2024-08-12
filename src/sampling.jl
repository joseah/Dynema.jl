function sample_index(rng::AbstractRNG, n::Integer)

    index = collect(1:n)
    index = sample(rng, index, length(index); replace=true)
    return index

end


function resample(rng::AbstractRNG, data::AbstractDataFrame)
    index = collect(1:n)
    index = sample(rng, index, length(index); replace=true)
    return index
end


function resample(rng::AbstractRNG, data::AbstractDataFrame, type::Symbol = :two_stage; cluster::String)

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

    return(boot_id)

end