# ---------------------------------------------------------------------------- #
#                           Define DynemaModel struct                          #
# ---------------------------------------------------------------------------- #

struct DynemaModel

    f::FormulaTerm
    bterm::String
    ncells::Int
    ndonors::Int
    stats::DataFrame
    B::Vector{Int64}
    bootdist::AbstractVector

end

# ---------------------------------------------------------------------------- #
#                               Define accessors                               #
# ---------------------------------------------------------------------------- #

"""

`f(::Dynema.DynemaModel)`

Extract formula used for a DynemaModel
"""

f(m::DynemaModel) = m.f

"""

`bterm(::Dynema.DynemaModel)`

Extract formula used for a DynemaModel
"""

bterm(m::DynemaModel) = m.bterm

"""

`ncells(::Dynema.DynemaModel)`

Extract number of cells used for a DynemaModel
"""

ncells(m::DynemaModel) = m.ncells

"""

`ndonors(::DynemaModel)`

Extract number of donors/individuals for a DynemaModel
"""

ndonors(m::DynemaModel) = m.ndonors

"""

`stats(::Dynema.DynemaModel)`

Extract all summary statistics for a DynemaModel
"""

stats(m::DynemaModel) = m.stats


"""

`bstat(::Dynema.DynemaModel)`

Extract bootstrapepd statistic for a DynemaModel
"""

bstat(m::DynemaModel) = m.stats.bstat

"""

`coefs(::Dynema.DynemaModel)`

Extract OLS beta coefficients for all SNPS for the tested bootstrapped 'bterm' with a DynemaModel
"""

coefs(m::DynemaModel) = m.stats.b

"""

`pvalues(::Dynema.DynemaModel)`

Extract empirical p-values for a DynemaModel
"""

pvalues(m::DynemaModel) = m.stats.p


"""

`snps(::Dynema.DynemaModel)`

Extract SNP/genetic variant names provided as column names in genotypying data from a DynemaModel
"""

snps(m::DynemaModel) = m.stats.snp


"""

`B(::Dynema.DynemaModel)`

Extract number of bootstrap iterations applied iteratively 
for a DynemaModel
"""

B(m::DynemaModel) = m.B


"""

`bootdist(::Dynema.DynemaModel)`

Extract bootstrap stat distributions for each SNP for a DynemaModel
"""

bootdist(m::DynemaModel) = m.bootdist


# ---------------------------------------------------------------------------- #
#                                Define printing                               #
# ---------------------------------------------------------------------------- #

function Base.show(io::IO, ::MIME"text/plain", m::DynemaModel)

    print(Crayon(foreground = :light_yellow, bold = true), "\nDynamic non-parametric eQTL mapping (Dynema) model\n\n")
    print(Crayon(foreground = :green), "Wild cluster bootstrap via WildBootTests.jl\n\n")
    print(Crayon(foreground = :blue), f(m), "\n\n")


    print(Crayon(reset = true, bold = true), "Term tested   = ")
    println(Crayon(foreground = :red, bold = true), bterm(m))

    
    print(Crayon(reset = true, bold = true), "N. bootstraps = ")
    println(Crayon(foreground = :red, bold = true), "$(sum(B(m)))")

    
    print(Crayon(reset = true, bold = true), "N. SNPs       = ")
    println(Crayon(foreground = :red, bold = true), "$(nrow(stats(m)))")

    
    print(Crayon(reset = true, bold = true), "N. cells      = ")
    println(Crayon(foreground = :red, bold = true), "$(ncells(m))")
    
    print(Crayon(reset = true, bold = true), "N. donors     = ")
    println(Crayon(foreground = :red, bold = true), "$(ndonors(m))")


    if nrow(stats(m)) >= 10
        
        glance = first(sort(stats(m), [order(:p), order(:stat, by = abs, rev = true)]), 10)
        push!(glance, fill("...", ncol(stats(m))), promote = true)

    else

        glance = stats(m)

    end

    println(Crayon(reset = true), "\nResults")
    pretty_table(glance, header = (names(glance)))
    println("** smallest p-value = $(2/sum(B(m))); report as p < $(2/sum(B(m)))")
    
    
end
