
function calculate_bootstrap_pvalue(x)
    min_zero = minimum([sum(x .> 0), sum(x .< 0)])

    if min_zero == 0
        pval = 2 / length(x)
    else
        pval = 2 * min_zero / length(x)
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

function count_nonzeros(x)
    minimum([sum(x .> 0), sum(x .< 0)])
end