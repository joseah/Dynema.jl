```@meta
CurrentModule = Dynema
```

# Mapping interaction eQTL effects (context-depedent)


Besides an eQTL's **main effect**, Dynema can map context-depedent eQTLs by testing **interaction effect(s)**. In this example, we can map interaction effects with each of the cell states with a **single-context interaction** test or all cell states of interests jointly with a **multi-contetext interaction** test.


```@setup interaction
include("src/assets/simulate_setup.jl")
```

**Single-context interaction**

Let's consider the question: is the effect of an eQTL different between non-cytotoxic and cytotoxic T cells? We can answer this question by testing an interaction with C1 (T cell cytotoxicity). In this case, C1 is a single-cell level context variable with mean= 0 and variance = 1, where positive values reflect more cytotoxicity.


Let's define the formula:

```@example interaction
f_GXC1 = @formula(0 ~ 1 + G + C1 + C2 + C3 + G & C1)
```

Notice that the interaction term **G x C1** is denoted as `G & CV1`. 

Next, we can simply run `map_locus` providing the interaction formula and specifying the term we want to test, in this case `G & CV1`.


```@example interaction
res_gxc1 = map_locus(f_GXC1;
        pheno = expr,
        geno = ex_geno,
        meta = meta,
        groups = meta.donor_id,
        termtest = ["G & C1"])
res_gxc1
```

Let's extract the summary statistics and write into a text file (optional)

```@example interaction
summ_gxc1 = get_summary(res_gxc1)
CSV.write("summ_gxc1.tsv", summ_gxc1, delim = "\t")
```


Assuming input variants are ordered by genomic location, let's make a quick
locus plot:


```@example interaction
using Plots

p_gxc1 = scatter(1:nrow(summ_gxc1), -log10.(summ_gxc1.p), 
        xlabel = "Genomic location", ylabel = " -log(p)", 
        title = "G x C1 interaction", legend = false)
```

---


Let's consider now mapping single-context interaction with Treg/activation status


```@example interaction
f_GXC2 = @formula(0 ~ 1 + G + C1 + C2 + C3 + G & C2)

res_gxc2 = map_locus(f_GXC2;
        pheno = expr,
        geno = ex_geno,
        meta = meta,
        groups = meta.donor_id,
        termtest = ["G & C2"])
res_gxc2
```


Let's make a second locus plot for G x C1:


```@example interaction
summ_gxc2 = get_summary(res_gxc2)

p_gxc2 = scatter(1:nrow(summ_gxc2), -log10.(summ_gxc2.p), 
        xlabel = "Genomic location", ylabel = " -log(p)", 
        title = "G x C2 interaction", legend = false, color = "orange")
```


Let's compare results


```@example interaction
plot(p_gxc1, p_gxc2)
```


**Multi-context interaction**


We can also test multiple contexts at once as follows:



```@example interaction
f_int = @formula(0 ~ 1 + G + C1 + C2 + C3 + G & C1 + G & C2 + G & C3)

res_int = map_locus(f_int;
        pheno = expr,
        geno = ex_geno,
        meta = meta,
        groups = meta.donor_id,
        termtest = ["G & C1", "G & C2", "G & C3"])
res_int
```



```@example interaction
summ_int = get_summary(res_int)

p_int = scatter(1:nrow(summ_int), -log10.(summ_int.p), 
        xlabel = "Genomic location", ylabel = " -log(p)", 
        title = "Multi-context interaction", legend = false, color = "green")
```


The multi-context interaction captures both G x C1 and G x C2 signals.








