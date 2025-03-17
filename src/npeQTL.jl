module npeQTL

using MixedModels
using StatsBase
using Distributed
using DataFrames
using Random
using CategoricalArrays
using GLM
using StatsModels
using LinearAlgebra
using StaticArrays
using NearestNeighborDescent
using Distances

export calculate_poisson_disks
export map_locus
export boot_snp
export calculate_pvalue
export count_nonzeros
export adjust_phenotype
export basic_p
export interp_pval


include("adjust_phenotype.jl")
include("bootstrapping.jl")
include("mapping.jl")
include("pvalue.jl")
include("poisson_disk_sampling.jl")

end
