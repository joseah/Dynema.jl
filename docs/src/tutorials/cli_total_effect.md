# Total effect (main + interaction)

A **total effect** tests a variant's main effect `G` together with its
`G & context` interaction term(s) in a single joint test -- useful when
you want to know whether a variant has *any* effect on gene expression,
whether constant across cells or context-dependent.

## Demo dataset

These commands run end to end against a small, fully public demo dataset
hosted in the [Dynema_datasets](https://github.com/joseah/Dynema_datasets)
repository -- a fully simulated dataset (100 donors, ~4000 cells, 100
variants, gene `CTSS`, three cell-state contexts, plus a few synthetic
covariates). Download it once:

```bash
input=demo_data
mkdir -p "$input"

# Downloads the whole data/cli_demo/ directory in one shot -- including
# the VCF's pre-built tabix index -- via the repository's tarball, rather
# than fetching each file individually:
curl -L https://github.com/joseah/Dynema_datasets/archive/refs/heads/main.tar.gz | \
  tar -xz -C "$input" --strip-components=3 "Dynema_datasets-main/data/cli_demo"
```

## Run it

`--effect total` tests `G` and all three `G & context` interaction terms
(named via `--interaction-with`) together:

```bash
input=demo_data

./bin/dynema-map \
  --expr-prefix "$input/expr" \
  --meta "$input/meta.tsv" \
  --vcf "$input/genotypes.vcf.gz" \
  --bed "$input/CTSS.bed" \
  --window 500000 \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --contexts C1,C2,C3 \
  --interaction-with C1,C2,C3 \
  --donor-col donor_id \
  --cell-id-col cell_id \
  --effect total \
  --out CTSS_total.tsv
```

This is the command-line equivalent of building the formula
`@formula(0 ~ 1 + G + C1 + C2 + C3 + <covariates...> + G & C1 + G & C2 + G & C3)`
by hand and calling [`map_locus`](@ref) with
`termtest = ["G", "G & C1", "G & C2", "G & C3"]`.

This run reads the same gene, variants, and cells as [Main
effect](cli_main_effect.md) and [Interaction
effect](cli_interaction_effect.md).
