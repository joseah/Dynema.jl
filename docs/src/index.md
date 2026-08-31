```@meta
CurrentModule = Dynema
```

# Dynema.jl

*Dynamic single-cell eQTL mapping.*

Dynema tests for genetic associations with single-cell gene expression,
including associations that change across a continuous or categorical cell
state (a *dynamic* or genotype × cell-state eQTL). It models single-cell
expression counts as a function of SNP genotype dosage, cell state(s), and
genotype × cell-state interaction(s) with a Poisson generalized linear
model (GLM), while controlling for sample- and single-cell-level covariates.

Cluster-robust variance estimators (CRVE) account for the non-independence
of cells drawn from the same donor. Because parametric p-values from a GLM
can be poorly calibrated in this setting, Dynema avoids relying on
distributional assumptions by instead computing p-values with a
**score bootstrap**: single-cell-level contributions to the score function
are bootstrapped with wild weights applied at the donor level, and p-values
are derived from the empirical distribution of the resulting bootstrapped,
standardized coefficients.

## Installation

Dynema.jl is not yet registered in the General registry. Install it
directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/joseah/Dynema.jl")
```

## Getting started

- [Set up simulated data](tutorials/simulation.md) downloads a simulated
  single-cell dataset and prepares it for mapping.
- [Mapping main eQTL effects](tutorials/main_effect.md) shows how to test a
  context-independent (main) eQTL effect.
- [Mapping interaction eQTL effects](tutorials/interaction.md) shows how to
  test a genotype × cell-state interaction term, single- or multi-context.
- [API Reference](functions.md) documents every exported function.


## Getting help

Please open an issue on [GitHub](https://github.com/joseah/Dynema.jl/issues)
for bug reports or feature requests.
