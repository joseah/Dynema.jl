```@meta
CurrentModule = Dynema
```

# Dynema.jl

See our pre-print using this [link](https://www.biorxiv.org/content/10.64898/2026.08.25.747138v1)

*Dynema* (Dynamic eQTL mapping in single cells) is a fast and calibrated method to map eQTL effects genome-wide at true single-cell resolution. We argue that an eQTL effect is different for each individual cell. Dynema can decompose such effects into: i) context-independent (main), ii) context-dependent (interaction), and iii) total (main and interaction).

## Installation

You can install it directly from GitHub:

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
