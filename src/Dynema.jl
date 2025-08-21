module Dynema

using WildBootTests
using StatsBase
using Distributed
using DataFrames
using Random
using CategoricalArrays
using StatsModels
using LinearAlgebra
using StaticArrays
using NearestNeighborDescent
using Distances
using Distributions
using PrettyTables
using ProgressMeter

export map_locus
export adjust_phenotype
export f, fterm, ncells, ndonors, stats, bstat, coefs, snps, B, bootdist


include("adjust_phenotype.jl")
include("mapping.jl")
include("pvalue.jl")
include("Dynema_struct.jl")
include("utils.jl")

end
