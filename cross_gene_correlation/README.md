# cross_gene_correlation

For each sample, these scripts compute the cross-gene Spearman and Pearson correlations between observed expression and each arm (every model track and every gene-mean baseline). Within a cohort, all arms use the same genes and samples.

| File | Description |
|------|-------------|
| `01_score_all_cohorts.R` | Scores every model track and the gene-mean baselines in Fig. 1, after dropping samples whose expression barely correlates with their own cohort's leave-one-out gene mean (two COVID-19 samples). Writes `arm_scores_allcohorts_qcexact.csv` |
| `02_score_personal_vs_reference.R` | Scores every model track with personalized and with reference-genome predictions. Writes `arm_scores_allcohorts_withref.csv` |
| `03_paired_stats.R` | The 48 paired t-tests in the text and Supp. Table 1 (every model vs every gene-mean baseline in each cohort), at participant level (a participant's samples averaged) |

A sample is kept only if it has genotypes and a usable prediction from all four models, so every model is scored on the same samples. The genes are those present in every gene-mean source, FANTOM5 included, resulting in 13,061 brain genes and 10,562 blood genes. 
