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
using Distributions
using PrettyTables
using ProgressMeter
using StableRNGs
using CSV
using CodecZlib
using Statistics
using htslib_jll

export map_locus
export expand_geno, expand_genotypes
export get_f, get_termtest, get_ncell, get_ndonor, get_summary
export get_stat, get_p, get_variant, get_B, get_bootdists, get_time
export get_stattype, get_testtype, get_boot
export get_pos, get_gene, get_chr
export set_pos!, set_gene!, set_chr!
export extract_geno_dataframe
export extract_gene_expression, resolve_mtx_triplet


include("Dynema_struct.jl")
include("ExpandedGeno.jl")
include("mapping.jl")
include("pvalue.jl")
include("utils.jl")
include("bootstrap.jl")
include("vcf_genotypes.jl")
include("matrix_market_expression.jl")

end
