module npeQTL


using MixedModels
using StatsBase
using Distributed
using DataFrames

# Write your package code here.
export fit_model
export sample_index
export fit_lmm
export calculate_pvalue
export count_nonzeros

include("fit.jl")
include("sampling.jl")
include("pvalue.jl")

end
