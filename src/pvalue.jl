
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



function crvetest(R, r; resp::AbstractVector, scores::AbstractMatrix, betas::AbstractVector, 
                    A::AbstractMatrix, clustid::AbstractMatrix, imposenull::Bool, small::Bool = false)

    p_analytical = if imposenull

        scoretest(R, r; 
                    resp = resp, 
                    scores = scores,
                    beta = betas,
                    A = A,
                    clustid = clustid, 
                    ml = true,
                    scorebs = true,
                    small = small)

    else

        waldtest(R, r; 
                    resp = resp, 
                    scores = scores,
                    beta = betas,
                    A = A,
                    clustid = clustid, 
                    ml = true,
                    scorebs = true,
                    small = small)
    end


    return p_analytical

end