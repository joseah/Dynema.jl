# Extra: Reusing an extracted genotype matrix: `dynema-extract-geno`

- **`dynema-extract-geno`** (optional): pulls a donor x variant dosage matrix
  for one gene's cis-window out of a genome-wide VCF and writes it to a file, for
  when you want to reuse the same extracted genotypes across multiple
  `dynema-map` runs (or inspect/QC them directly) instead of re-querying
  the VCF every time.

`--vcf` on `dynema-map` extracts on the fly for one run; `dynema-extract-geno`
does the same underlying extraction but writes the result to a file, for
reuse across multiple runs on the same gene or direct inspection/QC:

```bash
./bin/dynema-extract-geno \
  --vcf genotypes.vcf.gz \
  --bed BRCA1.bed --window 250000 \
  --field auto \
  --out BRCA1_geno.tsv
```

The single-gene bed-like file (one data row; columns chr, start, end, gene,
strand) specifies which gene's cis-window to extract; the TSS is derived the
same way FastQTL does it: the gene's start position on the plus strand, its
end position on the minus strand.

The resulting file is in exactly the shape `dynema-map --geno` expects. See
`./bin/dynema-extract-geno --help` for sample-id remapping (`--samples`)
and MAF/missingness filtering (`--maf`, `--max-missing`).
