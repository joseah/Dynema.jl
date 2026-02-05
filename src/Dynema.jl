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

export map_locus
export expand_geno
export get_f, get_termtest, get_ncell, get_ndonor, get_summary
export get_stat, get_p, get_snp, get_B, get_bootdists, get_time
export get_stattype, get_testtype, get_boot
export get_pos, get_gene, get_chr
export set_pos!, set_gene!, set_chr!


include("Dynema_struct.jl")
include("ExpandedGeno.jl")
include("mapping.jl")
include("pvalue.jl")
include("utils.jl")

end
