```@meta
CurrentModule = Dynema
```

# Download simulated data for demo

To demonstrate *Dynema*'s capabilities, let's download a simulated dataset
of T cells from `https://github.com/joseah/Dynema_datasets` -- the same
dataset used by the [command-line interface tutorials](../tutorials/command_line.md).
There, gene `CTSS`'s expression and genotypes for the variants around it
are read directly from a Matrix Market export and a VCF by `dynema-map`;
here, we do the same extraction ourselves from Julia, using the exported
functions that power `--expr-prefix`/`--vcf` under the hood.

Let's first import the single-cell metadata:

```@example setup

using Dynema, CSV, DataFrames, StatsModels, HTTP

base = "https://raw.githubusercontent.com/joseah/Dynema_datasets/refs/heads/main/data/cli_demo"

resp = HTTP.get("$base/meta.tsv")
meta = CSV.read(resp.body, DataFrame)

nothing # hide
```

The `meta` DataFrame has one row per cell, with columns:
1. `cell_id`: the cell id or barcode
2. `donor_id`: donor that each cell belongs to
3. `C1`: cell state 1 (cytotoxicity)
4. `C2`: cell state 2 (Treg/activation)
5. `C3`: cell state 3 (central memory)
6. `scaled_age`, `sex`: donor-level covariates (age is mean-0/variance-1 scaled across donors)
7. `scaled_log_nUMI`, `percent_mito`: single-cell-level covariates
8. `gPC1`-`gPC5`: donor-level genetic principal components
9. `ePC1`-`ePC5`: single-cell-level expression principal components

## Gene expression from a Matrix Market export

The dataset also ships gene `CTSS`'s expression as a Matrix Market export
(`.mtx.gz` + `.features.gz` + `.barcodes.gz`) -- as commonly produced by
Seurat/scanpy/10x pipelines for a whole count matrix, here trimmed to a
single gene. [`extract_gene_expression`](@ref) streams through the `.mtx`
file, pulling out just the requested gene without ever loading the whole
matrix into memory:

```@example setup

expr = mktempdir() do dir
    for suffix in ("mtx.gz", "features.gz", "barcodes.gz")
        open(joinpath(dir, "expr.$suffix"), "w") do io
            write(io, HTTP.get("$base/expr.$suffix").body)
        end
    end
    extract_gene_expression(
        mtx = joinpath(dir, "expr.mtx.gz"),
        features = joinpath(dir, "expr.features.gz"),
        barcodes = joinpath(dir, "expr.barcodes.gz"),
        gene = "CTSS",
        verbose = false,
    ).expr.CTSS
end

nothing # hide
```

`expr` is a vector containing the expression of `CTSS` for every cell in
`meta`, in the same order.

## Genotypes from a VCF

Genotypes for the variants around `CTSS` are similarly shipped as a
bgzipped, tabix-indexed VCF (`.vcf.gz` + a pre-built `.vcf.gz.tbi`).
[`extract_geno_dataframe`](@ref) queries the cis-window around a TSS with
`tabix` -- bundled by
[`htslib_jll`](https://github.com/JuliaBinaryWrappers/htslib_jll.jl), a
Julia package artifact, so no system `tabix`/`bcftools` install is needed
-- and converts each variant's dosage field to a donor x variant table:

```@example setup

geno = mktempdir() do dir
    vcf_path = joinpath(dir, "genotypes.vcf.gz")
    for suffix in ("", ".tbi")
        open(vcf_path * suffix, "w") do io
            write(io, HTTP.get("$base/genotypes.vcf.gz$suffix").body)
        end
    end

    extract_geno_dataframe(
        vcf = vcf_path,
        chr = "chr1",
        tss = 1011067,
        window = 1_000_000,
        donor_col = "donor_id",
        verbose = false,
    ).geno
end

nothing # hide
```

`geno` is a DataFrame with columns:

1. `donor_id`: the genotyped donor
2. `rs_snp_*`: allele dosage for the variants found in the queried window

Dynema requires all predictors at the single-cell level. In the case of genotying data, `Dynema` can create a memory effecient genotype representation to store and expand the genotype data to the single-cell level (i.e. assign genotypes to each donor's cell). For this task we use the `expand_genotypes` function with the following parameters:

1. A genotype matrix object with each column containing the allele dosage for single variant
2. A vector of donor ids, in the same order as the rows in the genotype matrix
3. A vector of donor ids for each cell/barcode. This vector should be in the same order as the metadata and gene expression at the single-cell level.
4. A vector with hte variant names (e.g. rs ids). Optional but strongly encourage to keep track of the variant names

```@example setup

donor_ids = geno.donor_id

# Drop donor_id column from genotype data
select!(geno, Not(:donor_id))

# Convert genotype data to matrix and 
ex_geno = expand_genotypes(Matrix(geno), donor_ids, meta.donor_id, names(geno))
nothing # hide
```

We have now all inputs in shape to use Dynema! Let's proceed to the next sections :).

## Note: plain text files instead

If your own data is already pre-extracted as plain tables rather than a
VCF/Matrix Market export -- or you'd just rather skip the extraction step
for a quick test -- the same dataset is also available as plain TSVs, and
loading it is simpler (no `htslib_jll`/tabix step, no streaming parse):

```julia
resp = HTTP.get("$base/expr.tsv")
expr = CSV.read(resp.body, DataFrame).CTSS

resp = HTTP.get("$base/genotypes.tsv")
geno = CSV.read(resp.body, DataFrame)
```

`expr` and `geno` come out in exactly the same shape either way, so
everything from `expand_genotypes` onward is identical regardless of which
path you took to get there. This isn't the priority path for these
tutorials (real single-cell eQTL data more often starts from a VCF and a
Matrix Market export, so that's what's shown above), but it's there if you
need it.
