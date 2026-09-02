# Interaction effect (context-dependent)

Besides a variant's **main effect**, Dynema can map context-dependent
eQTLs by testing **interaction effect(s)** -- whether a variant's effect on
gene expression differs across cell states (single-context, e.g. `G & C1`)
or several jointly (multi-context, e.g. `G & C1 + G & C2 + G & C3`). 

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

## Run it

`--interaction-with` is used to provide names of the context(s) we want to test interaction(s) for. 
Here, we provide three cell contexts (`C1`, `C2`, `C3`) to test multi-context interaction effect.

`--effect interaction` specifies that this is an interaction test.

```bash
input=demo_data


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

To test a single context's interaction instead (e.g. just `G & C1`), pass
`--interaction-with C1` on its own.

This is the command-line equivalent of building the formula
`@formula(0 ~ 1 + G + C1 + C2 + C3 + <covariates...> + G & C1 + G & C2 + G & C3)`
by hand and calling [`map_locus`](@ref) with
`termtest = ["G & C1", "G & C2", "G & C3"]` -- see [Interaction
effect](interaction.md) for what single- vs. multi-context interaction
terms mean and how to interpret them.

To also test `G`'s main effect alongside these interactions on the same
data, see [Total effect](cli_total_effect.md); for just the main effect on
its own, see [Main effect](cli_main_effect.md).
