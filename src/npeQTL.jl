module npeQTL


using MixedModels
using StatsBase
using Distributed
using DataFrames
using Random
using StatsModels

# Write your package code here.
export boot_model
export sample_index
export fit_model
export calculate_pvalue
export count_nonzeros
export boot_snp

include("fit.jl")
include("sampling.jl")
include("pvalue.jl")
include("boot.jl")

end
