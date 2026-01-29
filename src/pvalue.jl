
function pass(stat, boot, stattype)

    if stattype == "z"
        minimum([sum(stat .> boot), sum(stat .< boot)])
    elseif stattype == "χ²"
        sum(boot .>= stat)
    end
    
end

function compute_pvalue(stat, boot, stattype)
    
    pass_obs = pass(stat, boot, stattype)

    if stattype == "z"
        
        pass_obs == 0 ? 2 / (length(boot) + 1) : 2 * pass_obs / length(boot)

    elseif stattype == "χ²"

        pass_obs == 0 ? 1 / (length(boot) + 1) : pass_obs / length(boot)
    
    end

end

function approximate_pvalue(x)

    p = map(x) do d
        p = cdf(Normal(Float64(d), std(x) / length(x)), 0)

        if mean(x) < 0
            p = 1 .- p
        end

        return p / length(x)
    end

    p = 2 * sum(p)
    
    if p == 0
        p = floatmin(Float64)
    end

    return(p)

end