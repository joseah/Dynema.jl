# Shared data-loading code for the tutorials. Included (not run standalone)
# via a `@setup` block on any page that needs `meta`, `expr`, and `ex_geno`
# without re-showing this code (it's already shown, block by block, on the
# "Set up simulated data" tutorial page).

using Dynema, CSV, DataFrames, StatsModels, HTTP

url = "https://raw.githubusercontent.com/joseah/Dynema_datasets/refs/heads/main/simulated_cell_expression.txt"
resp = HTTP.get(url)
data = CSV.read(resp.body, DataFrame)
meta = data[:, [:cell_id, :C1, :C2, :C3, :donor_id]]
expr = data.E

url = "https://raw.githubusercontent.com/joseah/Dynema_datasets/refs/heads/main/simulated_donor_genotypes.txt"
resp = HTTP.get(url)
geno = CSV.read(resp.body, DataFrame)

donor_ids = geno.donor_id
select!(geno, Not(:donor_id))
ex_geno = expand_genotypes(Matrix(geno), donor_ids, meta.donor_id, names(geno))
