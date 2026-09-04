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
    const imposenull::Bool
    const boot::Bool
    const stattype::String
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
get_p(m::DynemaModel) = m.summary.p

"""
`get_variant(::Dynema.DynemaModel)`

Extract variant names provided as column names in genotypying data from a DynemaModel
"""
get_variant(m::DynemaModel) = m.summary.variant

"""
`get_B(::Dynema.DynemaModel)`

Extract number of bootstrap iterations applied iteratively
for a DynemaModel
"""
get_B(m::DynemaModel) = m.B

"""
`get_bootdists(::Dynema.DynemaModel)`

Extract bootstrap stat distributions for each variant
"""
get_bootdists(m::DynemaModel) = m.bootdists

"""
`get_time(::Dynema.DynemaModel)`

Extract total elapsed time in seconds
"""
get_time(m::DynemaModel) = m.time

"""
`get_stattype(::Dynema.DynemaModel)`

Extract statistic type (z or χ²)
"""
get_stattype(m::DynemaModel) = m.stattype

"""
`get_testtype(::Dynema.DynemaModel)`

Whether a score or wald test were performed
"""
get_testtype(m::DynemaModel) = m.imposenull ? "Score/lagrange multiplier" : "Wald"

"""
`get_boot(::Dynema.DynemaModel)`

Determine whether bootstrapping was perfomed
"""
get_boot(m::DynemaModel) = m.boot

"""
`get_pos(::Dynema.DynemaModel)`

Extract genomic position for each variant
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
`set_pos!(::Dynema.DynemaModel)`

Sets positions for all variants for a DynemaModel
"""
function set_pos!(m::DynemaModel, pos::Union{Nothing, Vector{Int64}, Vector{Float64}})
    m.pos = pos
    return m
end

"""
`set_gene!(::Dynema.DynemaModel)`

Sets gene name for a DynemaModel
"""
function set_gene!(m::DynemaModel, gene::Union{Nothing, String})
    m.gene = gene
    return m
end

"""
`set_chr!(::Dynema.DynemaModel)`

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

    print(io, Crayon(foreground = :light_yellow, bold = true), "\nDynamic eQTL mapping (Dynema) model\n\n")
    print(io, Crayon(foreground = :blue), get_f(m), "\n\n")

    if !isnothing(get_gene(m))
            print(io, Crayon(reset = true, bold = true), "Gene name    = ")
            println(io, Crayon(foreground = :green, bold = true), get_gene(m))
    end

    if !isnothing(get_chr(m))
            print(io, Crayon(reset = true, bold = true), "Gene chr.    = ")
            println(io, Crayon(foreground = :green, bold = true), get_chr(m))
    end


    print(io, Crayon(reset = true, bold = true), "Term(s)   = ")
    println(io, Crayon(foreground = :red, bold = true), get_termtest(m))


    if m.boot

    print(io, Crayon(reset = true, bold = true), "N. bootstraps = ")
    println(io, Crayon(foreground = :green, bold = true), "$(sum(get_B(m)))")

    end


    print(io, Crayon(reset = true, bold = true), "N. variants   = ")
    println(io, Crayon(foreground = :green, bold = true), "$(nrow(get_summary(m)))")


    print(io, Crayon(reset = true, bold = true), "N. cells      = ")
    println(io, Crayon(foreground = :green, bold = true), "$(get_ncell(m))")

    print(io, Crayon(reset = true, bold = true), "N. donors     = ")
    println(io, Crayon(foreground = :green, bold = true), "$(get_ndonor(m))")

    print(io, Crayon(reset = true, bold = true), "Test type     = ")
    println(io, Crayon(foreground = :green, bold = true), "$(get_testtype(m))")

    summ = get_summary(m)

    if nrow(summ) >= 10

        glance = first(sort(summ, [order(:p), order(get_stattype(m), by = abs, rev = true)]), 10)
        push!(glance, fill("...", ncol(summ)), promote = true)

    else

        glance = summ

    end

    # With many covariates/contexts, `summ` can have dozens of columns --
    # too wide to print without either wrapping mid-row (illegible) or
    # relying on terminal-width auto-cropping (unreliable once stdout is
    # redirected, e.g. for --log). Instead, always show variant/stat/p plus
    # whichever column(s) were actually tested (`termtest`); everything else
    # (intercept, covariates, untested contexts) is still in `summ`/--out,
    # just not in this console preview.
    tt_cols = get_termtest(m) isa AbstractVector ? get_termtest(m) : [get_termtest(m)]
    key_cols = filter(c -> c in names(glance), unique(vcat(["variant", get_stattype(m), "p", "p_boot", "p_boot_approx"], tt_cols)))

    if length(key_cols) < ncol(glance)

        n_hidden = ncol(glance) - length(key_cols)
        shown = glance[:, key_cols]
        println(io, Crayon(reset = true), "\nResults (variant/stat/p/tested term(s); ",
                "$n_hidden more column(s) -- covariates, untested contexts -- in the full results/--out file)")
        pretty_table(io, shown, header = (names(shown)))

    else

        println(io, Crayon(reset = true), "\nResults")
        pretty_table(io, glance, header = (names(glance)))

    end

    if m.boot
        println(io, "** empirical p_boot is floored at $(2/sum(get_B(m))) (report as p < $(2/sum(get_B(m)))); " *
                    "p_boot_approx extrapolates below it via a FastQTL-style beta approximation " *
                    "of the bootstrap distribution\n")
    end

    print(io, Crayon(reset = true, bold = true), "Computation time = ")
    println(io, Crayon(foreground = :green, bold = true), "$(round(get_time(m) / 60, sigdigits = 4)) mins.")

    # Reset terminal styling: without this, the green/bold state above leaks
    # into whatever the caller prints next (e.g. the CLI's following
    # bullets/sections).
    print(io, Crayon(reset = true))

end
