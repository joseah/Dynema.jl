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
using StableRNGs
using CovarianceMatrices

export map_locus
export get_f, get_termtest, get_ncell, get_ndonor, get_summary, get_stat, get_coef, get_p, get_snp, get_B, get_bootdists
export get_pos, get_gene, get_chr
export set_pos!, set_gene!, set_chr!
export calculate_poisson_disks, aggregate_expr, aggregate_meta, uniq

include("Dynema_struct.jl")
include("mapping.jl")
include("pvalue.jl")
include("poisson_disk_sampling.jl")
include("aggregation.jl")
include("utils.jl")

end
