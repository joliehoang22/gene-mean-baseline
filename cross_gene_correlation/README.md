# cross_gene_correlation

For each sample, these scripts compute the cross-gene Spearman and Pearson correlations between observed expression and each arm (every model track and every gene-mean baseline). Within a cohort, all arms use the same genes and samples.

| File | Description |
|------|-------------|
| `01_score_all_cohorts.R` | Scores every arm in all four cohorts, after dropping samples whose expression barely correlates with their own cohort's leave-one-out gene mean (two COVID-19 samples). Writes `arm_scores_allcohorts_qcexact.csv` |
| `02_score_personal_vs_reference.R` | Adds the reference-genome arm. Writes `arm_scores_allcohorts_withref.csv` |
| `03_wilcoxon_training_genemean.R` | Paired Wilcoxon tests of each model against the gene-mean baseline from its own training data (Supp. Fig. 3) |
| `04_paired_stats.R` | Paired t-tests quoted in the text, at participant level (a participant's samples averaged), plus every model vs every outside baseline |

A sample is kept only if it has genotypes and a usable prediction from all four models, so every model is scored on the same samples. The genes are those present in every gene-mean source, resulting in 13,061 brain genes and 10,562 blood genes.
