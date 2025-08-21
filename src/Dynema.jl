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
export f, bterm, ncells, ndonors, sumstats, bstats, coefs, pvalues, snps, B, bootdists


include("adjust_phenotype.jl")
include("mapping.jl")
include("pvalue.jl")
include("Dynema_struct.jl")
include("utils.jl")

end
