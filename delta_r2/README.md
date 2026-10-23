# delta_r2

For each sample, we fit two linear models across genes, with all variables z-scored within the sample:

```
null:  observed ~ gene_mean
full:  observed ~ gene_mean + model
ΔR² = R²full − R²null
```

The gene-mean baselines here always come from other cohorts, because a cohort's own leave-one-out baseline would be circular.

| File | Description |
|------|-------------|
| `01_delta_r2_discovery.R` | LBP LIV brain (GTEx, ROSMAP, MSBB baselines) and LBP LIV blood (GTEx, COVID-19). Draws Fig. 4, writes the Supp. Table 1 summary |
| `02_delta_r2_replication.R` | ROSMAP (GTEx, MSBB, LBP PM) and COVID-19 (GTEx, LBP LIV blood). Writes the Supp. Table 2 summary |
| `03_delta_r2_reference_genome.R` | Repeats 01 and 02 with reference-genome predictions in place of personalized ones, on the same samples, genes and baselines. Draws Supp. Figs 5–6 |

All three use the shared gene panel from `cross_gene_correlation` (13,061 brain / 10,562 blood genes). The medians at the end of `01` and `02` are the values in Supp. Tables 1 and 2.
