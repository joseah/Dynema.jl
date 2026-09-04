module Dynema

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
export extract_gene_expression, resolve_mtx_triplet, prepare_gene_expression

# WildBootTests.jl is an *optional* (weak) dependency: Dynema's own direct
# CRVE score test covers standard analytical mapping, and only score
# bootstrapping (`boot = true`) and a few fallback/self-check paths need the
# library. The DynemaWildBootTestsExt extension (loaded automatically when
# WildBootTests is installed and `using WildBootTests` has run alongside
# Dynema) provides the real `crvetest`/`scorebootstrap` methods and flips
# this flag; without it, those paths error with install instructions (see
# the stubs in pvalue.jl/bootstrap.jl) and the fast path skips its one-time
# library self-check.
const HAS_WILDBOOTTESTS = Ref(false)


include("Dynema_struct.jl")
include("ExpandedGeno.jl")
include("mapping.jl")
include("pvalue.jl")
include("utils.jl")
include("bootstrap.jl")
include("vcf_genotypes.jl")
include("matrix_market_expression.jl")

end
