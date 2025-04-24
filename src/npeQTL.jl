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
using Distributions
using PrettyTables

export calculate_poisson_disks
export map_locus
export adjust_phenotype


include("adjust_phenotype.jl")
include("bootstrapping.jl")
include("mapping.jl")
include("pvalue.jl")
include("poisson_disk_sampling.jl")
include("bootstrapmap.jl")

end
