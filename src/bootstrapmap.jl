struct BootstrapMap

    coefs::DataFrame
    p_boot::DataFrame
    p_approx::DataFrame
    p_analytical::DataFrame
    std_analytical::DataFrame
    var_comp::DataFrame
    boot::DataFrame
    boot_iter::Vector{Int64}
    formula::FormulaTerm
    snps::AbstractVector
    contexts::AbstractVector

end


function Base.show(io::IO, ::MIME"text/plain", bm::npeQTL.BootstrapMap)
    
    
    n_snps = length(bm.snps)
    n_contexts = length(bm.contexts)
    n_boot = sum(bm.boot_iter)
    info = DataFrame(:c => ["SNPs", "Contexts", "Max. bootstrap iters."], :N => [n_snps, n_contexts, n_boot])
    print(Crayon(foreground = :light_yellow, bold = true), "\nDynamic non-parametric eQTL mapping (DynNEMA) model\n\n")
    print(Crayon(reset = true))
    println(bm.formula, "\n")
    pretty_table(info, header = (["Parameters", "N."]))
    
    
end