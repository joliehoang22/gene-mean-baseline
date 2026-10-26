# gene_mean_baseline

These scripts compute the gene-mean baselines that the model predictions are compared against.

| File | Description |
|------|-------------|
| `01_expression_brain.R` | Compute the log10(TPM + 1)-transformed gene expression for brain samples from GTEx, MSBB, ROSMAP, LBP LIV, and LBP PM |
| `02_expression_blood.R` | Compute the log10(TPM + 1)-transformed gene expression for whole blood samples from GTEx, LBP LIV blood and Mount Sinai COVID-19 |
| `03_fantom5_gene_means.R` | Compute the gene-mean baseline from the FANTOM5 CAGE dataset: mean over the 7 brain and 4 blood CAGE libraries behind the model tracks |

We quantify gene expression in transcript per million (TPM) because CAGE counts one tag per transcript, so it does not depend on transcript length, and TPM removes transcript length from RNA-seq counts. 
