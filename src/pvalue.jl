
function pass(stat, boot)
    minimum([sum(stat .> boot), sum(stat .< boot)])
end

function compute_pvalue(stat, boot)

    pass_obs = pass(stat, boot)
    pval = pass_obs == 0 ? 2 / (length(boot) + 1) : 2 * pass_obs / length(boot)
    return pval

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