# ---------------------------------------------------------------------------- #
#                           Define DynemaModel struct                          #
# ---------------------------------------------------------------------------- #

struct DynemaModel

    f::FormulaTerm
    term::String
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

`term(::Dynema.DynemaModel)`

Extract formula used for a DynemaModel
"""

term(m::DynemaModel) = m.term

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

Extract summary statistics for a DynemaModel
"""

stats(m::DynemaModel) = m.stats

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
    println(Crayon(foreground = :red, bold = true), term(m))

    
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
