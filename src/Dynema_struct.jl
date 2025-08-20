struct DynemaModel

    f::FormulaTerm
    term::String
    ncells::Int
    ndonors::Int
    summ_stats::DataFrame
    B::Vector{Int64}
    boot_dist::AbstractVector

end


function Base.show(io::IO, ::MIME"text/plain", m::DynemaModel)

    print(Crayon(foreground = :light_yellow, bold = true), "\nDynamic non-parametric eQTL mapping (Dynema) model\n\n")
    print(Crayon(foreground = :green), "Wild cluster bootstrap via WildBootTests.jl\n\n")
    print(Crayon(foreground = :blue), m.f, "\n\n")


    print(Crayon(reset = true, bold = true), "Term tested   = ")
    println(Crayon(foreground = :red, bold = true), m.term)

    
    print(Crayon(reset = true, bold = true), "N. bootstraps = ")
    println(Crayon(foreground = :red, bold = true), "$(sum(m.B))")

    
    print(Crayon(reset = true, bold = true), "N. SNPs       = ")
    println(Crayon(foreground = :red, bold = true), "$(nrow(m.summ_stats))")

    
    print(Crayon(reset = true, bold = true), "N. cells      = ")
    println(Crayon(foreground = :red, bold = true), "$(m.ncells)")
    
    print(Crayon(reset = true, bold = true), "N. donors     = ")
    println(Crayon(foreground = :red, bold = true), "$(m.ndonors)")


    if nrow(m.summ_stats) >= 6
        glance = first(m.summ_stats, 3)
        push!(glance, fill("...", ncol(m.summ_stats)), promote = true)
        glance = vcat(glance, last(m.summ_stats, 3))
        

    else
        glance = m.summ_stats
    end



    println(Crayon(reset = true), "\nResults")
    pretty_table(glance, header = (names(glance)))
    println("** smallest p-value = $(2/sum(B)); report as p < $(2/sum(B))")
    
    
end