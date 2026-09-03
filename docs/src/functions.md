```@meta
CurrentModule = Dynema
```

# API Reference

This page documents every function exported by Dynema.jl.

## Mapping

The core entry point for single-cell eQTL mapping.

```@docs
map_locus
```

## Genotype containers

A memory-efficient, lazily-indexed representation of donor-level genotype
data expanded to the single-cell level.

```@docs
expand_geno
expand_genotypes
```

## Model results

`map_locus` returns a `Dynema.DynemaModel`. The accessors below extract
individual fields without needing to know the internal struct layout.

```@docs
get_f
get_termtest
get_ncell
get_ndonor
get_summary
get_stat
get_p
get_variant
get_B
get_bootdists
get_time
get_stattype
get_testtype
get_boot
get_pos
get_gene
get_chr
```

## Model mutators

Positional/gene/chromosome metadata can be attached to a fitted model after
the fact (for example, once you know where a variant maps in the genome).

```@docs
set_pos!
set_gene!
set_chr!
```

## Command-line extraction helpers

The extraction behind the [`dynema-map`/`dynema-extract-geno` command-line
tools](tutorials/command_line.md) (`--vcf`/`--expr-prefix`) is implemented
as regular exported functions, callable directly from Julia without going
through either CLI wrapper.

```@docs
extract_geno_dataframe
extract_gene_expression
resolve_mtx_triplet
prepare_gene_expression
```
