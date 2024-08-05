function calculate_pvalue(x)
    2 * minimum([sum(x .> 0) + 1, sum(x .< 0) + 1]) / (length(x) + 1)
end



function count_nonzeros(x)
    minimum([sum(x .> 0), sum(x .< 0)])
end
