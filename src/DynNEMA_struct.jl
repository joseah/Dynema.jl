struct DynNEMAModel

    coefs::DataFrame
    p_boot::DataFrame
    p_approx::DataFrame
    p_analytical::DataFrame
    std_analytical::DataFrame
    var_comp::DataFrame
    boot::Vector{DataFrame}
    boot_iter::Vector{Int64}
    formula::FormulaTerm
    snps::AbstractVector
    contexts::AbstractVector

end


function Base.show(io::IO, ::MIME"text/plain", bm::DynNEMAModel)
    
    
    n_snps = length(bm.snps)
    n_contexts = length(bm.contexts)
    
    print(Crayon(foreground = :light_yellow, bold = true), "\nDynamic non-parametric eQTL mapping (DynNEMA) model\n\n")

    if length(bm.boot_iter) == 1 && bm.boot_iter[1] == 0
        info = DataFrame(:c => ["SNPs", "Contexts"], :N => [n_snps, n_contexts])
        print(Crayon(foreground = :red, bold = true), "No bootstrapping: Parametric results only\n\n")
    else
        n_boot = sum(bm.boot_iter)
        info = DataFrame(:c => ["SNPs", "Contexts", "Max. bootstrap iters."], :N => [n_snps, n_contexts, n_boot])

    end


    print(Crayon(reset = true))
    println(bm.formula, "\n")
    pretty_table(info, header = (["Parameters", "N."]))
    
    
end