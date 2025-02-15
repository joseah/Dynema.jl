
function map_locus_interactions(snps::AbstractMatrix, expr::AbstractVector, 
    contexts::AbstractMatrix, donor::AbstractVector, batch::AbstractVector;
    n::Vector{Int64} = [100, 400, 500, 4000, 5000], main_effect::Bool = false)
    
    

   # ------------- Validate dimensionality of input data structures ------------- #

    snp_size = size(snps)
    expr_len = length(expr)
    contexts_size = size(contexts)
    donor_len = length(donor)
    batch_len = length(batch)


    n_obs = Dict(:snp => snp_size[1], :expr => expr_len, 
                 :contexts => contexts_size[1], :donor => donor_len, 
                 :batch => batch_len)
   
   n_obs_vals  = collect(values(n_obs))

   if(!all(x -> x == n_obs_vals[1], n_obs_vals))
       
       println("Verify number of observatios of input files")
       println("Number of observations per input data\n")
       for (key, value) in pairs(n_obs)
           println(key, " = ", value)
           
       end
       throw("Number of observations in input data do not match between input files")
       
   end


   # ----------------- Assign dummy variable names for modelling ---------------- #

   # Add context information
   context_names = Symbol.("C" .* string.(1:contexts_size[2]))
   contexts = DataFrame(contexts, context_names)

   # Add snp information
   snps = DataFrame(snps, :auto)


   design = DataFrame(E = expr, donor = CategoricalArray(donor), batch = CategoricalArray(batch))
   design = hcat(design, contexts)

   # Add contexts
   for col_name in names(contexts)
       design[!, col_name] = contexts[!, col_name]
   end


   # Create formula
   f = term(:E) ~ term(:G) * sum(term.(context_names)) + 
       (term(1) | term(:donor)) + (term(1) | term(:batch))

   boot_terms = Symbol.("G & " .* string.(context_names))


   if(main_effect)
    boot_terms = vcat(:G, boot_terms)
   end
   
   
   # Add genotype information for SNp of interest

   res = map(eachcol(snps)) do snp
       boot_snp(f, snp, design, boot_terms, n)
   end
   

   res = Dict(:coefs => reduce(vcat, [x[:coefs] for x in res]), 
              :p     => reduce(vcat, [x[:p] for x in res]))


   return(res)


end






function boot_snp(f::FormulaTerm, snp::Vector{Float64}, data::AbstractDataFrame, 
    boot_terms::Vector{Symbol}, boot_sizes::Vector{Int64})

    # ------------------------------- Add genotypes ------------------------------ #
    set = deepcopy(data)

    set[!, :G] = snp

    # ------------------- Compute betas and variance components ------------------ #

    # Fit model
    model = fit(MixedModel, f, set)
    
    # Extract betas and variance components
    coefs = DataFrame([model.betas])
    
    sigmas = vcat([x[1] for x in model.sigmas], model.sigma)
    var_values = sigmas.^2
    sigmas_names = vcat(collect(string.(keys(model.sigmas))), "residual")
    sigmas_names = Symbol.("var_" .* sigmas_names)

    # Merge betas and variance components
    re_var = DataFrame(transpose(var_values), sigmas_names)
    coefs = hcat(coefs, re_var)

    # ----------------------- Start adaptive bootstrapping ----------------------- #

    boot = DataFrame()

    for n in boot_sizes
       
        boot_i = nonparametricbootstrap(MersenneTwister(n), n, model; progress = false)
        boot_i = DataFrame(boot_i.β)
        boot_i = unstack(boot_i, :coefname, :β)
        boot_i = boot_i[! , Not(:iter)]

        # Combine with prebious results
        boot = vcat(boot, boot_i)
        # Return logical vector specifying which parameters pass
        # condition
        param_status = map(count_nonzeros, eachcol(boot[:, boot_terms]))
        # Stop bootstrapping if more than 20 values cross zero
        if(all(param_status .> 20))
            break
        end

    end

    boot = boot[:, boot_terms]
    p_vals = map(calculate_pvalue, eachcol(boot))
    p_vals = DataFrame(transpose(p_vals), names(boot))

    res = Dict(:coefs => coefs, :p => p_vals)

    res[:coefs]
    res[:p]

    return(res)

end




function calculate_pvalue(x)
    2 * minimum([sum(x .> 0) + 1, sum(x .< 0) + 1]) / (length(x) + 1)
end



function count_nonzeros(x)
    minimum([sum(x .> 0), sum(x .< 0)])
end
