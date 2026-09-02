```@meta
CurrentModule = Dynema
```

# Main effect (context-independent)


A traditional eQTL mapping focuses on testing a variant's **main effect** on gene expression. This effect is constant and **indepedendent** of any context(s). In the case of single-cell data, it reflects an average or baseline genetic effect across the whole cell population. 

`Dynema` can map main eQTL effects. For demonstration, let's use the simulation data introduced in the [simulation section](simulation.md).

```@setup main_effect
include("src/assets/simulate_setup.jl")
```

First, let's define a formula to specify our eQTL mappin strategy:


```@example main_effect
f_main = @formula(0 ~ 1 + G + C1 + C2 + C3 +
                       scaled_age + sex + scaled_log_nUMI + percent_mito +
                       gPC1 + gPC2 + gPC3 + gPC4 + gPC5 +
                       ePC1 + ePC2 + ePC3 + ePC4 + ePC5)
```

To map a main effect, the minimal terms required are:

- `1`: intercept
- `G`: genetic main effect

> Note: 
> 1. `Dynema` exlusively reserves the use of the term `G` to encode the genetic variant term. When using Dynema, please do not name any column as `G` as it might trigger some errors.
> 2. The response coding does not matter, for convention we use `0 ~`.

Additionally, we can include covariates to control for biological and experimental sources of confounding. For Dynema's framework, we need to encode all variables at the single-cell-level (i.e. a value for every single cell) -- a donor-level covariate like `scaled_age` just repeats the same value across all of that donor's cells.

In this example, we add cell state covariate terms for `C1-3` to account for differential expression confounding, plus donor-level covariates (`scaled_age`, `sex`, `gPC1-5`) and single-cell-level covariates (`scaled_log_nUMI`, `percent_mito`, `ePC1-5`) -- the same covariates the [command-line interface](../tutorials/cli_main_effect.md) passes via `--covariates`.


`Dynema`' main function is `map_locus`. Because this is a generalized function, we can map any type eQTL effect with it by just changing the formula. Let's provide the following arguments to `map_locus`:

1. `f`: formula
2. `pheno`: gene expression at single cell level (counts
3. `geno`: allele dosage or probabilities. Ideally processed using `expand_genotypes()` for memory efficiency
4. `meta`: Single-cell level metadata including donor and single-cell covariates, and contexts (e.g. cell states)
5. `groups`: donor ids for each cell
6. `termtest`: term that we want to test for. In this case, we are
interested in testing the main effect `G`. `termtest` must name a term that appears in the formula's right-hand


Let's run map_locus with the previous arguments:


```@example main_effect
res_main = map_locus(f_main;
        pheno = expr,
        geno = ex_geno,
        meta = meta,
        groups = meta.donor_id,
        termtest = "G")
res_main
```

