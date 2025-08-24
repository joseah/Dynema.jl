
function aggregator(pds::AbstractDataFrame, 
    data::Union{AbstractDataFrame, AbstractMatrix},
    f::Function)

    # Iterate over disks
    data_aggr = @showprogress map(1:nrow(pds)) do i
        # Extract cells from ith disk (skip first two rows)
        pd = collect(pds[i, 3:end])
        # Remove missing values (cells assigned randomly to other disks``)
        pd = filter(!ismissing, collect(pd))
        # Aggregate data across 
        # all cells for each gene contained in the ith disk and store results
        cell2diskaggr(pd, data, f)

    end

    return data_aggr

end



function cell2diskaggr(pd::AbstractVector, 
                        data::Union{AbstractDataFrame, AbstractMatrix}, 
                        f::Function)

   res = if isa(data, AbstractMatrix)
    try
        collect(f(data[:, pd]))
    catch e
        throw("Verify that cells in target data were used for poisson disk generation")
    end

   else
        subdata = length(pd) > 1 ? data[pd, :] : data[[pd], :]
        [f(col) for col in eachcol(subdata)]
   end

   return res

end




function aggregate_expr(pds::AbstractDataFrame, expr::NamedArray, f::Function = x -> sum(x, dims = 2))

    # Aggregate gene expression data
    data_aggr = aggregator(pds, expr, f)

    # Gather results across all disks
    data_aggr = reduce(hcat, data_aggr)

    # Create named array with gene and disks names 
    data_aggr = NamedArray(data_aggr, 
                            names = (names(expr, 1), pds.disk_id), 
                            dimnames = ("features", "disks")) 

    return data_aggr

end



function aggregate_meta(pds::AbstractDataFrame, meta::AbstractDataFrame, f::Function = mean)

    # Create vector of vectors to store results
    data_aggr  = aggregator(pds, meta, f)

    # Gather results across all disks
    data_aggr = DataFrame(permutedims(reduce(hcat, data_aggr)), names(meta))

    return data_aggr

end


function uniq(x::AbstractArray)

    u = unique(x)
    if length(u) == 1
        u = first(u)
    else
        throw("Column contains more than one unique value")
    end


    return u

end