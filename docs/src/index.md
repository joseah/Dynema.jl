```@meta
CurrentModule = Dynema
```

# Dynema.jl

See our pre-print using this [link](https://www.biorxiv.org/content/10.64898/2026.08.25.747138v1)

*Dynema* (Dynamic eQTL mapping in single cells) is a fast and calibrated method to map eQTL effects genome-wide at true single-cell resolution. We argue that an eQTL effect is different for each individual cell. Dynema can decompose such effects into: i) context-independent (main), ii) context-dependent (interaction), and iii) total (main and interaction).

**No Julia coding required:** Dynema ships a command-line interface that runs the full mapping workflow -- straight from a VCF and a Matrix Market count matrix (or plain TSV/CSV files) to a summary statistics table -- entirely from the terminal. See [CLI: Command-line overview](tutorials/command_line.md) to get started. A Julia API is also available for anyone who wants to call Dynema directly from their own scripts or pipelines -- see the [API Reference](functions.md).

## Getting started

- [Installing Julia](tutorials/installing_julia.md) walks through getting
  Julia set up if you don't have it yet.
- [CLI: Command-line interface overview](tutorials/command_line.md) familiarizes with 
  running *Dynema* and its parameters from a VCF and a Matrix Market count 
  matrix -- or plain TSV/CSV files -- without writing any Julia code.
- [CLI: Main effect](tutorials/cli_main_effect.md),
  [CLI: Interaction effect](tutorials/cli_interaction_effect.md), and
  [CLI: Total effect](tutorials/cli_total_effect.md) each run one of the
  single-cell eQTL effects using a small demo dataset, entirely from the command 
  line.
- [API Reference](functions.md) documents every exported function, for
  anyone calling Dynema directly from Julia.


## Getting help

Please open an issue on [GitHub](https://github.com/joseah/Dynema.jl/issues)
for bug reports or feature requests.
