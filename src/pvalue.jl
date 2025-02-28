function interp_pval(q)

    R = length(q)
    tstar = sort(q)
   zero = searchsortedlast(tstar, 0)

    if(zero == 0 | zero == R) 
        pval = 2/R
    else
        pval = 2*min(zero/R, (R-zero)/R)
        
    end

    return(pval)

end
   

function basic_p(obs, boot; null = 0)

    interp_pval(2*obs .- boot .- null)

end


function calculate_pvalue(x)
    2 * minimum([sum(x .> 0) + 1, sum(x .< 0) + 1]) / (length(x) + 1)
end


