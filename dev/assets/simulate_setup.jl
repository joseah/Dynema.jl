# Shared data-loading code for the tutorials. Included (not run standalone)
# via a `@setup` block on any page that needs `meta`, `expr`, and `ex_geno`
# without re-showing this code (it's already shown, block by block, on the
# "Set up simulated data" tutorial page).

using Dynema, CSV, DataFrames, StatsModels, HTTP

base = "https://raw.githubusercontent.com/joseah/Dynema_datasets/refs/heads/main/data/cli_demo"

resp = HTTP.get("$base/meta.tsv")
meta = CSV.read(resp.body, DataFrame)

# Gene expression: download the Matrix Market export and pull out CTSS with
# `extract_gene_expression` -- the same function `dynema-map --expr-prefix`
# uses under the hood (see the command-line interface tutorials).
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

# Genotypes: download the VCF (plus its pre-built tabix index -- .tbi
# files are binary and don't round-trip well through git diffs, but
# they're small and static here, so it's shipped alongside the VCF rather
# than rebuilt on every run) and extract with `extract_geno_dataframe` --
# the same function `dynema-map --vcf` uses under the hood.
ex_geno = mktempdir() do dir
    vcf_path = joinpath(dir, "genotypes.vcf.gz")
    for suffix in ("", ".tbi")
        open(vcf_path * suffix, "w") do io
            write(io, HTTP.get("$base/genotypes.vcf.gz$suffix").body)
        end
    end

    r = extract_geno_dataframe(
        vcf = vcf_path,
        chr = "chr1",
        tss = 1011067,
        window = 1_000_000,
        donor_col = "donor_id",
        verbose = false,
    )
    donor_ids = r.geno.donor_id
    select!(r.geno, Not(:donor_id))
    expand_genotypes(Matrix(r.geno), donor_ids, meta.donor_id, names(r.geno))
end
