# ---------------------------------------------------------------------------- #
#                           Define DynemaModel struct                          #
# ---------------------------------------------------------------------------- #

mutable struct DynemaModel

    const f::FormulaTerm
    const termtest::Union{String, Vector{String}}
    const ncells::Int
    const ndonors::Int
    const summary::DataFrame
    const B::Vector{Int64}
    const bootdists::AbstractVector
    const time::Float64
    pos::Union{Nothing, Vector{Real}}
    gene::Union{Nothing, String}
    chr::Union{Nothing, String, Int}

end

# ---------------------------------------------------------------------------- #
#                               Define accessors                               #
# ---------------------------------------------------------------------------- #

"""

`get_f(::Dynema.DynemaModel)`

Extract formula used for a DynemaModel
"""

get_f(m::DynemaModel) = m.f

"""

`get_termtest(::Dynema.DynemaModel)`

Extract formula used for a DynemaModel
"""

get_termtest(m::DynemaModel) = m.termtest

"""

`get_ncell(::Dynema.DynemaModel)`

Extract number of cells used for a DynemaModel
"""

get_ncell(m::DynemaModel) = m.ncells

"""

`get_ndonor(::DynemaModel)`

Extract number of donors/individuals for a DynemaModel
"""

get_ndonor(m::DynemaModel) = m.ndonors

"""

`get_summary(::Dynema.DynemaModel)`

Extract all summary statistics for a DynemaModel
"""

get_summary(m::DynemaModel) = m.summary


"""

`get_stat(::Dynema.DynemaModel)`

Extract bootstrapepd statistic for a DynemaModel
"""

get_stat(m::DynemaModel) = m.summary.stat


"""

`get_p(::Dynema.DynemaModel)`

Extract empirical p-values for a DynemaModel
"""

get_p(m::DynemaModel) = m.summary.p_boot


"""

`get_snp(::Dynema.DynemaModel)`

Extract SNP/genetic variant names provided as column names in genotypying data from a DynemaModel
"""

get_snp(m::DynemaModel) = m.summary.snp


"""

`get_B(::Dynema.DynemaModel)`

Extract number of bootstrap iterations applied iteratively 
for a DynemaModel
"""

get_B(m::DynemaModel) = m.B


"""

`get_bootdists(::Dynema.DynemaModel)`

Extract bootstrap stat distributions for each SNP for a DynemaModel
"""

get_bootdists(m::DynemaModel) = m.bootdists


"""

`get_time(::Dynema.DynemaModel)`

Extract total elapsed time in seconds
"""

get_time(m::DynemaModel) = m.time

"""

`get_pos(::Dynema.DynemaModel)`

Extract genomic position for each SNP
"""

get_pos(m::DynemaModel) = m.pos


"""

`get_gene(::Dynema.DynemaModel)`

Extract name for tested gene
"""

get_gene(m::DynemaModel) = m.gene



"""

`get_chr(::Dynema.DynemaModel)`

Extract name for tested gene
"""

get_chr(m::DynemaModel) = m.chr


"""

`set_pos(::Dynema.DynemaModel)`

Sets positions for all SNPs/genetic variants for a DynemaModel
"""


function set_pos!(m::DynemaModel, pos::Union{Nothing, Vector{Int64}, Vector{Float64}})
    m.pos = pos
    return m
end



"""

`set_gene(::Dynema.DynemaModel)`

Sets gene name for a DynemaModel
"""


function set_gene!(m::DynemaModel, gene::Union{Nothing, String})
    m.gene = gene
    return m
end

"""

`set_chr(::Dynema.DynemaModel)`

Sets chromosome name for gene tested
"""


function set_chr!(m::DynemaModel, chr::Union{Nothing, String, Int})
    m.chr = chr
    return m
end



# ---------------------------------------------------------------------------- #
#                                Define printing                               #
# ---------------------------------------------------------------------------- #

function Base.show(io::IO, ::MIME"text/plain", m::DynemaModel)

    print(Crayon(foreground = :light_yellow, bold = true), "\nDynamic eQTL mapping (Dynema) model\n\n")
    print(Crayon(foreground = :green), "Score bootstrap via WildBootTests.jl\n\n")
    print(Crayon(foreground = :blue), get_f(m), "\n\n")

    if !isnothing(get_gene(m))
            print(Crayon(reset = true, bold = true), "Gene name    = ")
            println(Crayon(foreground = :green, bold = true), get_gene(m))
    end

    if !isnothing(get_chr(m))
            print(Crayon(reset = true, bold = true), "Gene chr.    = ")
            println(Crayon(foreground = :green, bold = true), get_chr(m))
    end


    print(Crayon(reset = true, bold = true), "Term(s) tested   = ")
    println(Crayon(foreground = :red, bold = true), get_termtest(m))

    
    print(Crayon(reset = true, bold = true), "N. bootstraps = ")
    println(Crayon(foreground = :red, bold = true), "$(sum(get_B(m)))")

    
    print(Crayon(reset = true, bold = true), "N. SNPs       = ")
    println(Crayon(foreground = :red, bold = true), "$(nrow(get_summary(m)))")

    
    print(Crayon(reset = true, bold = true), "N. cells      = ")
    println(Crayon(foreground = :red, bold = true), "$(get_ncell(m))")
    
    print(Crayon(reset = true, bold = true), "N. donors     = ")
    println(Crayon(foreground = :red, bold = true), "$(get_ndonor(m))")

    summ = get_summary(m)[:, .!occursin.("p_naive", names(get_summary(m)))]

    if nrow(summ) >= 10
        
        glance = first(sort(summ, [order(:p), order(:stat, by = abs, rev = true)]), 10)
        push!(glance, fill("...", ncol(summ)), promote = true)

    else

        glance = summ

    end

    println(Crayon(reset = true), "\nResults")
    pretty_table(glance, header = (names(glance)))
    println("** smallest p-value computed= $(2/sum(get_B(m))); report as p < $(2/sum(get_B(m)))\n")
    
    print(Crayon(reset = true, bold = true), "Computation time = ")
    println(Crayon(foreground = :green, bold = true), "$(round(get_time(m) / 60, sigdigits = 4)) mins.")

end
