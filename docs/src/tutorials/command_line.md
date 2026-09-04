# Command-line interface (CLI)

See [Installing Julia](installing_julia.md) if you don't have Julia
installed yet.

**You don't need to write any Julia to run *Dynema*!**. The `bin/` folder of the
[Dynema.jl repository](https://github.com/joseah/Dynema.jl) ships a
self-contained command-line tool:

- **`dynema-map`**: This tool maps one or more genes (a batch, defined by
  the rows of `--bed`) against their cis genetic variants -- reads gene
  expression, single-cell metadata, and genotypes, runs [`map_locus`](@ref)
  per gene, and writes one summary statistics table per gene.
- **`dynema-prepare-expr`** (optional, run once per expression matrix):
  builds a gene-major index (`.dgx`) next to a Matrix Market export, after
  which `dynema-map` loads any single gene in milliseconds instead of
  scanning the whole multi-GB matrix per run.
- **`dynema-prepare-bed`** (optional, run once per study): builds
  `dynema-map`'s bed-like gene file(s) from a GTF annotation -- keeping only
  genes present in the expression data's features file, picking an
  unambiguous identifier per gene, and optionally splitting into fixed-size
  chunks for HPC batching.

This script bootstraps its own Julia environment automatically the first
time it's run (installing `Dynema` plus a handful of CLI-only
dependencies), so all you need installed is Julia (>= 1.9)
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

- **Gene expression** (recommended: Matrix Market): the standard sparse export
  (`PREFIX.mtx`, `PREFIX.features`, `PREFIX.barcodes` -- as produced by
  Seurat/scanpy/10x pipelines), containing raw counts. Run `dynema-prepare-expr`
  on it once to build a `.dgx` index, and every later `dynema-map` run loads its
  gene in milliseconds. A plain text file (TSV/CSV with a cell-id column and one
  column per gene) is also accepted as an alternative for small or pre-extracted
  data.
- **Genotypes** (recommended: VCF): a bgzipped, tabix-indexed VCF with genotype
  dosages, from which the tested gene's cis-window variants are extracted in
  about a second. A plain text file (a pre-extracted donor-level dosage table:
  one donor-id column, one column per variant -- e.g. the output of
  `dynema-extract-geno`) is also accepted.
- **The gene(s) to map**: a bed-like plain text file (`--bed`) with columns
  chr, start, end, gene, strand -- one data row per gene. It serves two
  purposes at once: its gene column (a name/symbol or gene id) specifies
  *what* to map, and its positions/strand give each TSS, derived
  FastQTL-style -- the gene's start position on the plus strand, its end
  position on the minus strand. A multi-row file defines a batch: the genes
  are mapped one after another in the same run, each writing its own output
  file. (A small curated table, not a full GTF -- reading it stays instant.)
- **Metadata**: plain text file with cell id, donor id, cell-state
  contexts, and any donor or single-cell covariates, one row per cell.
- **Type of eQTL effect**: either main, interaction, or total.

The recommended route is Matrix Market expression (indexed once with
`dynema-prepare-expr`) plus a genome-wide VCF and a bed-like gene annotation:
all inputs are then compressed, indexed, and random-access, so per-gene runs
pay seconds of I/O, not minutes.

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
  --meta "$input/meta.tsv" \
  --vcf "$input/genotypes.vcf.gz" \
  --bed "$input/CTSS.bed" \
  --window 500000 \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --interaction-with C1,C2,C3 \
  --donor-col donor_id \
  --cell-id-col cell_id \
  --effect interaction \
  --out interaction
```

Arguments:

- `--expr-prefix` specifies the basename of the matrix market files (e.g. [expr].mtx.gz, 
  [expr].barcodes.gz, [expr].features.gz). If a [expr].dgx index built by
  `dynema-prepare-expr` sits next to them, the gene loads from it in
  milliseconds; otherwise the matrix is scanned once per run (minutes for
  multi-GB files -- index it!).
- `--vcf`: VCF file (*.vcf or *.vcf.gx) with **genotype dosages, not hard calls.** 
  The VCF must already carry per-sample genotype dosages (`DS`) and/or genotype 
  probabilities (`GP`) -- as produced by standard imputation pipelines (Minimac, IMPUTE2/5,
  Beagle, etc.). `--field auto` (the default) prefers `GP` over `DS` when a
  variant has both; hard-call genotypes (`GT`) are never read. 
  Must be accompanied by a tabix index file (*.tbi).
- `--bed`: bed-like file (plain or gzipped) specifying the gene(s) to map:
  one data row per gene with columns chr, start, end, gene, strand (a
  standard 6-column BED with a score column also works; header/`#` lines are
  skipped). Each row's gene column -- a gene name/symbol (`CTSS`) or a gene
  id (`ENSG00000163131`) -- names that gene: with a single-column features file
  (e.g. a Seurat export) it must match that column exactly; with a 10x
  features file (gene_id, gene_name, modality) it is searched against
  gene_name (column 2) first and, failing that, against gene_id (column 1) --
  no match is an error (Ensembl id version suffixes are ignored). The TSS is
  derived the same way FastQTL does: the gene's start position on the plus
  strand, its end position on the minus strand. The annotation's chromosome
  naming must match the VCF's (e.g. `chr1` vs `1`). A cis-window is built
  around the derived TSS.
- `--window`: *cis* region half-width in bp around the TSS. By default 500000 (0.5 Mb).
- `--meta`: plain text file with cell id, donor id, cell-state contexts, and any donor or 
  single-cell covariates, one row per cell.
- `covariates`: list of covariates. Must be comma-separated and should match column names in 
  metadata file.
- `donor-col`: Column name in metadata file including the donor ids. This must match the VCF donor ids too.
- `cell-col`: Column name in metadata file including the cell ids/barcodes. This must match the barcodes 
  provided in the expression data.
- `--effect`: Type(s) of single-cell eQTL to test, comma-separated: any of
  `main`, `interaction`, and `total` (e.g. `--effect main,interaction`).
  Multiple effects share each gene's data extraction and write separate
  files, `<out>_<effect>_<gene>.tsv`. In a multi-effect run, `main` uses the
  classic model without G × context terms.
- `--interaction-with`: comma-separated context column(s) in the metadata file to
  test G × context interactions for; their main effects are added to the model
  automatically. Required when `--effect` includes interaction/total. Contexts
  you only want to adjust for (without an interaction) belong in `--covariates`.
- `--out`: output *prefix*. Each gene writes its own summary statistics
  table named `<out>_<gene>.tsv` (or `<out><gene>.tsv` if the prefix ends in
  `/`; directories are created as needed; with no `--out`, just
  `<gene>.tsv`). This makes it easy to distinguish analyses of the same
  genes: `--out main` and `--out interaction` give `main_CTSS.tsv` and
  `interaction_CTSS.tsv`.


## Alternative: plain text files

For small datasets, quick checks, or genotypes pre-extracted with
`dynema-extract-geno`, gene expression and genotypes can also be passed as
plain TSV/CSV tables with `--expr`/`--geno` instead of
`--expr-prefix`/`--vcf` (no `--bed` needed: a pre-extracted genotype table
already fixes its own variants). Everything else about
`--meta`/`--covariates`/`--effect`/`--interaction-with` works
exactly the same way. For full-size data, prefer the indexed Matrix Market +
VCF route above.

Here's a simple example just to illustrate the parameter usage (do not run):

```bash
./bin/dynema-map \
  --expr "$input/expr.tsv" \
  --meta "$input/meta.tsv" \
  --geno "$input/genotypes.tsv" \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --interaction-with C1,C2,C3 \
  --donor-col donor_id \
  --cell-id-col cell_id \
  --effect interaction \
  --out interaction
```

## Mapping many genes: batches

Dynema maps genome-wide studies the same way FastQTL does: split the gene
list into chunks and let your HPC scheduler parallelize over chunks. Each
`dynema-map` invocation takes one bed-like chunk file and maps its genes
sequentially -- metadata, plain-text expression, the model formula, and (via
the `.dgx` index) the expression matrix are all loaded once per invocation,
so per-gene overhead stays low. A gene that fails (e.g. absent from the
features file) is reported and skipped without sinking the rest of the
batch, and each gene writes its own `<out>_<gene>.tsv`. For example, with a
genome-wide study:

```bash
# once: build 100-gene chunk beds from your GTF, matched to your features file
./bin/dynema-prepare-bed --gtf gencode.gtf.gz --features expr.features.gz \
  --chunk-size 100 --out beds/chunk

# then submit one job per chunk, e.g. (SLURM):
#   ./bin/dynema-map --bed beds/chunk_001.bed --expr-prefix expr --vcf genotypes.vcf.gz \
#     ... --out main --log main_chunk_001.log
```

Each invocation also writes `<out>_summary.tsv` -- one row per gene with its
lead variant and statistics -- so concatenating the chunk summaries gives the
study-wide top-associations table without parsing the per-gene files.

## Learn more

Run `./bin/dynema-map --help` for the complete list of options (adaptive score bootstrapping
with `--boot`, distributing variants across worker processes with `--parallel`, and more).


## Examples: main, interaction, and total effects

Three separate tutorials each walk through one of `--effect`'s modes:

- [Main effect (context-independent)](cli_main_effect.md)
- [Interaction effect (context-dependent)](cli_interaction_effect.md)
- [Total effect (main + interaction)](cli_total_effect.md)

Take a look at the section below to understand parameter specification for Dynema.

