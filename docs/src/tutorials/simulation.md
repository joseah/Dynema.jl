```@meta
CurrentModule = Dynema
```

# Download simulated data for demo

To demonstrate *Dynema*'s capabilities, let's download a simulated dataset of T cells from `https://github.com/joseah/Dynema_datasets`.

Let's first import the gene expression and metadata:

```@example setup

using Dynema, CSV, DataFrames, StatsModels, HTTP

url = "https://raw.githubusercontent.com/joseah/Dynema_datasets/refs/heads/main/simulated_cell_expression.txt"
resp = HTTP.get(url)
data = CSV.read(resp.body, DataFrame)
meta = data[:, [:cell_id, :C1, :C2, :C3, :donor_id]]
expr = data.E

nothing # hide
```


The `meta` DataFrame contains six columns:
1. `cell_id`: The cell id or barcode
3. `C1`: Cell state 1 (cytotoxicity)
4. `C2`: Cell state 2 (Treg/activation)
5. `C3`: Cell state 3 (central memory)
6. `donor_id`: Donor that each cell belongs to

`expr` is a vector containing the expression of a single gene


Let's now import the genotype data:

```@example setup

url = "https://raw.githubusercontent.com/joseah/Dynema_datasets/refs/heads/main/simulated_donor_genotypes.txt"
resp = HTTP.get(url)
geno = CSV.read(resp.body, DataFrame)

nothing # hide
```


`geno` is a DataFrame with columns:

1. `donor_id`: The genotyped donor
2. `rs_snp_*`: simulated allele dosage for 100 variants


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
