```@meta
CurrentModule = Dynema
```

# Total effect (main + interaction)

A **total effect** tests a variant's main effect `G` together with its
`G & context` interaction term(s) in a single joint test -- useful when
you want to know whether a variant has *any* effect on gene expression,
whether constant across cells or context-dependent. It combines the 
[main effect](main_effect.md) and [interaction effect](interaction.md) tests 
from the previous two tutorials into one `termtest`.

```@setup total_effect
include("src/assets/simulate_setup.jl")
```

As with the interaction formulas, we test all three contexts jointly here
(a multi-context interaction), together with `G`'s main effect:

```@example total_effect
f_total = @formula(0 ~ 1 + G + C1 + C2 + C3 +
                        scaled_age + sex + scaled_log_nUMI + percent_mito +
                        gPC1 + gPC2 + gPC3 + gPC4 + gPC5 +
                        ePC1 + ePC2 + ePC3 + ePC4 + ePC5 +
                        G & C1 + G & C2 + G & C3)
```

This is the same formula as the multi-context interaction test in
[Interaction effect](interaction.md) -- what changes is `termtest`, which
now names `G` *and* every interaction term:

```@example total_effect
res_total = map_locus(f_total;
        pheno = expr,
        geno = ex_geno,
        meta = meta,
        groups = meta.donor_id,
        termtest = ["G", "G & C1", "G & C2", "G & C3"])
res_total
```

Let's extract the summary statistics:

```@example total_effect
summ_total = get_summary(res_total)
first(summ_total, 5)
```

This is the same covariates the [command-line interface](cli_total_effect.md)
passes via `--covariates`/`--contexts`/`--interaction-with`, and the same
`termtest` `--effect total` builds under the hood.
