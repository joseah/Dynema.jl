# Main effect (context-independent)

A traditional eQTL mapping focuses on testing a variant's **main effect**
on gene expression -- constant and **independent** of any context(s).

## Demo dataset

These commands run end to end against a demo dataset
hosted in the [Dynema_datasets](https://github.com/joseah/Dynema_datasets)
repository. Download it once:

```bash
input=demo_data
mkdir -p "$input"

curl -L https://github.com/joseah/Dynema_datasets/archive/refs/heads/main.tar.gz | \
  tar -xz -C "$input" --strip-components=3 "Dynema_datasets-main/data/cli_demo"
```


The `meta.tsv` has one row per cell, with columns:
1. `cell_id`: the cell id or barcode
2. `donor_id`: donor that each cell belongs to
3. `C1`: cell state 1 (cytotoxicity)
4. `C2`: cell state 2 (Treg/activation)
5. `C3`: cell state 3 (central memory)
6. `scaled_age`, `sex`: donor-level covariates (age is mean-0/variance-1 scaled across donors)
7. `scaled_log_nUMI`, `percent_mito`: single-cell-level covariates
8. `gPC1`-`gPC5`: donor-level genetic principal components
9. `ePC1`-`ePC5`: single-cell-level expression principal components


## Run it

`--effect main` tests the main effect alone, adjusted for the covariates and contexts (no interactions):

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
  --donor-col donor_id \
  --cell-id-col cell_id \
  --effect main \
  --out main

```

> This is the command-line equivalent of building the formula
> `@formula(0 ~ 1 + G + C1 + C2 + C3 + <covariates...>)` by hand and calling
> [`map_locus`](@ref) with `termtest = "G"`.

To also test the `G & context` interaction terms on the same data instead,
see [Interaction effect](cli_interaction_effect.md) or [Total
effect](cli_total_effect.md).
