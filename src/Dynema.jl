module Dynema

using WildBootTests
using GLM
using StatsBase
using Distributed
using DataFrames
using Random
using CategoricalArrays
using StatsModels
using LinearAlgebra
using StaticArrays
using NamedArrays
using NearestNeighborDescent
using Distances
using Distributions
using PrettyTables
using ProgressMeter

export map_locus
export adjust_phenotype
export f, bterm, ncells, ndonors, sumstats, bstats, coefs, pvalues, snps, B, bootdists
export pos, genename, genechr
export setpos!, setgene!, setchr!
export calculate_poisson_disks, aggregate_expr, aggregate_meta, uniq


include("Dynema_struct.jl")
include("adjust_phenotype.jl")
include("mapping.jl")
include("pvalue.jl")
include("poisson_disk_sampling.jl")
include("aggregation.jl")
include("utils.jl")

end
