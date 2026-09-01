<p align="center">
  <img src="docs/src/assets/logo.svg" alt="Dynema.jl logo" width="220">
</p>

# Dynema (Dynamic eQTL mapping for single cells)

[![Build Status](https://github.com/joseah/Dynema.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/joseah/Dynema.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/joseah/Dynema.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://joseah.github.io/Dynema.jl/dev/)
[![Docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://joseah.github.io/Dynema.jl/stable/)

*Dynema* is a method to map single-cell eQTL effects at real cellular resolution.

*Dynema*'s generalized framework enables testing complex regulatory effects including:

- **Main effects**: A standard eQTL effect, independent of any context
- **Interaction effects**: eQTL effects that dependent/vary depending on one (**single-context**) or more (**multi-context**) contexts
- **Total effects**: Joint effect of main and interaction eQTL components. This effect captures any genetic signal driven by either main or interaction eQTL effects


*Dynema* scales to genome-wide analysis, accounts for repeated measurements (multiple cells per donor), and provides calibrated p-values by using cluster-robust variance estimators (CRVEs). Additionally, it provides robust inferences in extreme scenarios including (1) small number of donors or (2) extreme imbalanced number cells per donor, by leveraging the score bootstrapping, implemented in [WildBootTests.jl](https://github.com/droodman/WildBootTests.jl).


# Installation

*Dynema* can be installed in Julia as follows:

```julia
using Pkg
Pkg.add(url = "https://github.com/joseah/Dynema.jl")
```

For more details and tutorials, see the [documentation website](https://joseah.github.io/Dynema.jl/dev/).

# Command-line usage

You don't need to write any Julia to run *Dynema*: `bin/dynema-map` maps a single gene against one or more genetic variants directly from plain TSV/CSV files.

```bash
./bin/dynema-map \
  --expr      simulated_cell_expression.txt \
  --meta      simulated_cell_expression.txt \
  --geno      simulated_donor_genotypes.txt \
  --gene      E \
  --covariates "" \
  --contexts  C1,C2,C3 \
  --test      interaction \
  --interaction-with C1 \
  --out       gxc1_results.tsv
```

Requires only a working Julia installation (>= 1.6) -- the first run bootstraps its own environment automatically (installs `Dynema` plus `ArgParse`/`CSV`/`DataFrames`/`StatsModels`/`CategoricalArrays`). If Julia isn't on your `PATH` as `julia`, invoke it directly instead: `julia --project=bin bin/dynema_map.jl [options]`.

Input files:

- `--expr`: single-cell gene expression counts (cell-id column + one column per gene, or a single expression column). Alternatively, pass `--expr-prefix PREFIX --gene GENE` to pull one gene straight out of a Matrix Market count matrix (`PREFIX.mtx`, `PREFIX.features`, `PREFIX.barcodes`, each tried gzipped then plain) -- see "Extracting a gene from a Matrix Market file" below. `--expr` and `--expr-prefix` are mutually exclusive.
- `--meta`: single-cell metadata (cell id, donor id, cell-state contexts, and any covariates), one row per cell.
- `--geno`: donor-level genotype dosages (donor-id column + one column per variant); expanded to the single-cell level automatically. Alternatively, pass `--vcf genotypes.vcf.gz` instead of `--geno` to extract the tested gene's cis-window genotypes on the fly from a tabix-indexed VCF (see below) -- `--geno` and `--vcf` are mutually exclusive.

Run `./bin/dynema-map --help` for the full list of options, including `--test main|interaction|total`, `--boot` for score-bootstrapped p-values, and `--wald` to use a Wald test instead of the default score/Lagrange-multiplier test.

Every run also saves a full transcript -- the exact command invoked, plus everything printed to the console (progress, warnings, the results table) -- to a log file, `--log` (default: `--out` with its extension swapped for `.log`, e.g. `--out CTSS_results.tsv` -> `CTSS_results.log`). `dynema-extract-geno` does the same.

## Extracting genotypes from a VCF

`--vcf` on `dynema-map` uses the same on-the-fly extraction described below, so most of the time you don't need `dynema-extract-geno` as a separate step at all:

```bash
./bin/dynema-map \
  --expr simulated_cell_expression.txt --meta simulated_cell_expression.txt --gene E \
  --vcf genotypes.vcf.gz --tss-file genes.tsv \
  --window 250000 --contexts C1,C2,C3 --test interaction --interaction-with C1 \
  --out gxc1_results.tsv
```

`dynema-map`'s `--gene` doubles as the gene id to look up in `--tss-file` (or pass `--chr`/`--tss` directly instead). `dynema-extract-geno` below remains useful when you want the extracted matrix saved to a file -- e.g. to reuse across multiple `dynema-map` runs on the same gene, or to inspect/QC it directly -- rather than re-querying the VCF every time.

Real genotype data usually lives in a single genome-wide, bgzipped and tabix-indexed VCF (`.vcf.gz`) rather than a pre-extracted per-gene matrix. `bin/dynema-extract-geno` pulls the cis-window around one gene's TSS out of such a VCF and writes it in the donor x variant format `dynema-map --geno` expects, converting genotype dosages (`DS`) or genotype probabilities (`GP`) -- whichever the VCF provides -- to allele dosage. Hard-call genotypes (`GT`) are not read.

```bash
./bin/dynema-extract-geno \
  --vcf genotypes.vcf.gz \
  --chr chr17 --tss 43044295 --window 250000 \
  --field auto \
  --out BRCA1_geno.tsv
```

`--field auto` (the default) prefers `GP` over `DS` when a variant has both, falling back to `DS` if `GP` isn't present; pass `--field GP` or `--field DS` to force one. The TSS can also be looked up from an annotation file with `--tss-file genes.tsv --gene BRCA1` instead of `--chr`/`--tss`. See `./bin/dynema-extract-geno --help` for sample-id remapping (`--samples`), MAF/missingness filtering, and other options.

The VCF must already be tabix-indexed (a `.tbi`/`.csi` file next to it). `tabix` itself doesn't need to be installed on your system: this script depends on [`htslib_jll`](https://github.com/JuliaBinaryWrappers/htslib_jll.jl), which Julia downloads automatically as a package artifact (for essentially any platform) the same way it downloads `ArgParse`/`CSV`/etc. If you don't have `tabix` installed separately either, you can index a VCF with that same bundled copy:

```bash
julia --project=bin -e 'using htslib_jll; htslib_jll.bgzip() do exe; run(`$exe -c genotypes.vcf`) end' > genotypes.vcf.gz
julia --project=bin -e 'using htslib_jll; htslib_jll.tabix() do exe; run(`$exe -p vcf genotypes.vcf.gz`) end'
```

(or, more conveniently, install `htslib` system-wide via `conda install -c bioconda htslib` or `brew install htslib` and just run `bgzip`/`tabix` directly.)

A synthetic VCF for exercising this end-to-end (built from the tutorial's own `simulated_donor_genotypes.txt`, with deliberately mismatched sample ids and a handful of edge cases) lives in `test/fixtures/mock_vcf/`, along with a walkthrough of the exact expected output in `test/fixtures/mock_vcf/NOTES.md`.

## Extracting a gene from a Matrix Market file

Single-cell count matrices are often exported as Matrix Market (`.mtx`, optionally gzipped) with separate features (gene names, one per row) and barcodes (cell ids, one per column) files -- the standard Seurat/scanpy/10x layout. `--expr-prefix` on `dynema-map` reads just the tested gene out of one of these directly -- give it the shared filename prefix and it resolves and validates all three companion files for you (trying `PREFIX.mtx.gz`/`PREFIX.features.gz`/`PREFIX.barcodes.gz` first, then the plain, ungzipped names):

```bash
./bin/dynema-map \
  --expr-prefix expr_0.05 \
  --gene CTSS \
  --meta your_meta.tsv \
  --vcf genotypes.vcf.gz --chr chr1 --tss 150765778 --window 250000 \
  --contexts C1,C2,C3 --test interaction --interaction-with C1 \
  --out CTSS_results.tsv
```

That example expects `expr_0.05.mtx.gz`, `expr_0.05.features.gz`, and `expr_0.05.barcodes.gz` (or their non-gzipped equivalents) to exist; if any are missing, `dynema-map` errors up front naming exactly which file(s) it looked for and where, rather than failing deep into a (potentially minutes-long) matrix scan.

Worth understanding what this does and doesn't buy you: Matrix Market coordinate files store nonzero `(row, col, value)` triples usually sorted by column, not by gene, and carry no separate index -- there's no tabix-style trick that lets Dynema jump straight to one gene's data. So `--expr-prefix` still reads through the entire file once. What it avoids is loading the *whole matrix into memory* the way a generic Matrix Market reader would: it streams the file and keeps only the (typically sparse, so few) values belonging to the requested gene, discarding everything else as it goes. For a matrix with hundreds of millions of nonzero entries, expect a full pass to take low-single-digit minutes dominated by decompression and line parsing, but with memory use that stays roughly constant rather than scaling with the size of the whole matrix.

