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
| `01_delta_r2_discovery.R` | LBP LIV brain (GTEx, ROSMAP, MSBB, LBP PM and FANTOM5 baselines) and LBP LIV blood (GTEx, COVID-19, FANTOM5). Writes the tables for Fig. 2 and Supp. Fig. 5 and its rows of Supp. Table 3 |
| `02_delta_r2_replication.R` | ROSMAP (GTEx, MSBB, LBP PM, FANTOM5) and COVID-19 (GTEx, LBP LIV blood, FANTOM5). Writes the table for Supp. Fig. 5 and its rows of Supp. Table 3 |

Both use the shared gene panel from `cross_gene_correlation` (13,061 brain / 10,562 blood genes). The figures are drawn by `figures/Fig2.R` and `figures/SuppFig5.R`.
