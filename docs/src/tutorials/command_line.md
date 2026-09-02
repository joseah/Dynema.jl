# Command-line interface (CLI)

See [Installing Julia](installing_julia.md) if you don't have Julia
installed yet.

**You don't need to write any Julia to run *Dynema*!**. The `bin/` folder of the
[Dynema.jl repository](https://github.com/joseah/Dynema.jl) ships a
self-contained command-line tool:

- **`dynema-map`**: This tool maps a single gene against one or 
  more genetic variant -- reads gene expression, single-cell metadata, 
  and genotypes, runs [`map_locus`](@ref), and writes a
  summary statistics table.

This script bootstraps its own Julia environment automatically the first
time it's run (installing `Dynema` plus a handful of CLI-only
dependencies), so all you need installed is Julia (>= 1.6)
itself. 


```bash
curl -L https://github.com/joseah/Dynema.jl/archive/refs/heads/main.tar.gz | tar -xz
mv Dynema.jl-main Dynema.jl
cd Dynema.jl
./bin/dynema-map --help
```

If `julia` isn't on your `PATH`, invoke the underlying script directly
instead: `julia --project=bin bin/dynema_map.jl [options]`.

## What `dynema-map` needs

`dynema-map` inputs:

- **Gene expression**: this can be provided as Matrix Market format
  (`PREFIX.mtx`, `PREFIX.features`, `PREFIX.barcodes`–— commonly exported by 
  Seurat/scanpy/10x pipelines), or a plain text file (TSV/CSV with a cell-id column 
  and one column per gene).
- **Genotypes**: either a VCF file (tabix-indexed) to extract the tested gene's cis-window 
  variants,  or a plain text file (a pre-extracted donor-level dosage table: one donor-id
  column, one column per variant).
- **Metadata**: plain text file with cell id, donor id, cell-state
  contexts, and any donor or single-cell covariates, one row per cell.
- **Type of eQTL effect**: either main, interaction, or total.

## Running Dynema with VCF and a Matrix Market Format matrix

Real single-cell eQTL studies usually start from a single genome-wide,
bgzipped and tabix-indexed VCF (`.vcf.gz` + `.vcf.gz.tbi`/`.csi`) rather
than a pre-extracted per-gene genotype table, and from a Matrix Market
count export (`.mtx[.gz]` + `.features[.gz]` + `.barcodes[.gz]`). `dynema-map` 
can read both directly. 

Here's a simple example just to illustrate the parameter usage (do not run):

```bash
./bin/dynema-map \
  --expr-prefix "$input/expr" \
  --gene CTSS \
  --meta "$input/meta.tsv" \
  --vcf "$input/genotypes.vcf.gz" \
  --tss-file "$input/tss.tsv" \
  --window 500000 \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --contexts C1,C2,C3 \
  --interaction-with C1,C2,C3 \
  --donor-col donor_id \
  --cell-id-col cell_id \
  --effect interaction \
  --out CTSS_multi-interaction.tsv
```

Arguments:

- `--expr-prefix` specifies the basename of the matrix market files (e.g. [expr].mtx.gz, 
  [expr].barcodes.gz, [expr].features.gz)
- `--gene`: The name of the gene tested. The name provided is used to extract its
  expression from the expressin file inputs and must match a gene in [expr].features.gz. 
  This step might take some minutes the larger the expression data is. Expression data
  can be split by gene or in batches to load faster. Otherwise, we can use a plain-text file
  as shown in the next section.
- `--vcf`: VCF file (*.vcf or *.vcf.gx) with **genotype dosages, not hard calls.** 
  The VCF must already carry per-sample genotype dosages (`DS`) and/or genotype 
  probabilities (`GP`) -- as produced by standard imputation pipelines (Minimac, IMPUTE2/5,
  Beagle, etc.). `--field auto` (the default) prefers `GP` over `DS` when a
  variant has both; hard-call genotypes (`GT`) are never read. 
  Must be accompanied by a tabix index file (*.tbi).
- `--tss-file`: TSV/CSV gene-to-TSS lookup table (columns `gene_id`, `chr`,
  `tss` position) to look up `--gene`'s chromosome and transcription start site from --
  a cis-window is built around the looked-up position. Alternative to passing
  `--chr`/`--tss` directly.
- `--window`: *cis* region. By default 1Mb.
- `--meta`: plain text file with cell id, donor id, cell-state contexts, and any donor or 
  single-cell covariates, one row per cell.
- `covariates`: list of covariates. Must be comma-separated and should match column names in 
  metadata file.
- `donor-col`: Column name in metadata file including the donor ids. This must match the VCF donor ids too.
- `cell-col`: Column name in metadata file including the cell ids/barcodes. This must match the barcodes 
  provided in the expression data.
- `--contexts`: list of contexts to consider. Must be comma-separated and should match column names in 
  metadata file.
- `--effect`: Type of single-cell eQTL to test: main, interaction, or total
- `--interaction-with`: indicates the interactions we want test for. Only necessary if a total or interaction
effect is assessed. 
- `--out`: output file with summary statistics.


## Running Dynema plain text files

If you already have gene expression and genotypes as pre-extracted plain
TSV/CSV tables, pass them with
`--expr`/`--geno` instead of `--expr-prefix`/`--vcf` -- everything else
about `--meta`/`--covariates`/`--contexts`/`--effect`/`--interaction-with`
works exactly the same way.

Here's a simple example just to illustrate the parameter usage (do not run):

```bash
./bin/dynema-map \
  --expr "$input/expr.tsv" \
  --gene CTSS \
  --meta "$input/meta.tsv" \
  --geno "$input/genotypes.tsv" \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --contexts C1,C2,C3 \
  --interaction-with C1,C2,C3 \
  --donor-col donor_id \
  --cell-id-col cell_id \
  --effect interaction \
  --out CTSS_multi-interaction.tsv
```

## Learn more

Run `./bin/dynema-map --help` for the complete list of options (adaptive score bootstrapping
with `--boot`, distributing variants across worker processes with `--parallel`, and more).


## Examples: main, interaction, and total effects

Three separate tutorials each walk through one of `--effect`'s modes:

- [Main effect (context-independent)](cli_main_effect.md)
- [Interaction effect (context-dependent)](cli_interaction_effect.md)
- [Total effect (main + interaction)](cli_total_effect.md)

Take a look at the section below to understand parameter specification for Dynema.

