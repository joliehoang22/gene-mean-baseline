# gene_mean_baseline

These scripts prepare the measured gene expression that the model predictions are compared with. Scripts 01 and 02 put every cohort's RNA-seq data on the same scale and keep the genes found in all cohorts and expressed in GTEx (TPM ≥ 1 in at least 25% of GTEx samples). A cohort's gene-mean baseline is each gene's average expression across that cohort's samples; the scoring scripts calculate it from these files. Script 03 calculates the FANTOM5 gene-mean baseline.

| File | Description |
|------|-------------|
| `01_expression_brain.R` | Compute the log10(TPM + 1)-transformed gene expression for brain samples from GTEx, MSBB, ROSMAP, LBP LIV, and LBP PM |
| `02_expression_blood.R` | Compute the log10(TPM + 1)-transformed gene expression for whole blood samples from GTEx, LBP LIV blood and Mount Sinai COVID-19 |
| `03_fantom5_gene_means.R` | Compute the gene-mean baseline from the FANTOM5 CAGE dataset: mean over the 7 brain and 4 blood CAGE libraries behind the model tracks |

We quantify gene expression in transcript per million (TPM) because CAGE counts one tag per transcript, so it does not depend on transcript length, and TPM removes transcript length from RNA-seq counts. 
