# gene_mean_baseline

These scripts compute the observed expression for each cohort and the gene-mean baselines. Run them in numeric order.

| File | Description |
|------|-------------|
| `01_expression_brain.R` | log10(TPM + 1) for GTEx, MSBB, ROSMAP, LBP LIV and LBP PM brain; expressed-gene universe. Writes `Null/mats.rds`, `Null/genes.rds` |
| `02_expression_blood.R` | Same for GTEx whole blood, LBP LIV blood and COVID-19. Writes `Null/mats_blood.rds`, `Null/genes_blood.rds` |
| `03_fantom5_gene_means.R` | FANTOM5 CAGE gene-mean baseline: mean over the 7 brain and 4 blood CAGE libraries behind the model tracks |

We quantify gene expression in transcript per million (TPM) because CAGE counts one tag per transcript, so it does not depend on transcript length, and TPM removes transcript length from RNA-seq counts. 
