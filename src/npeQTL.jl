module npeQTL

using MixedModels
using StatsBase
using Distributed
using DataFrames
using Random
using StatsModels
using LinearAlgebra
using MixedModels: fixef!, getθ!
using StaticArrays
using ProgressMeter

export nonparametricbootstrap
export resample!
export residuals_from_blups
export inflation_factor

export map_locus_interactions
export boot_snp
export calculate_pvalue
export count_nonzeros

export adjust_phenotype


include("adjust_phenotype.jl")
include("bootstrapping.jl")
include("mapping.jl")

end
