# figures

There is one script per figure. Each reads tables written by `cross_gene_correlation` or `delta_r2` (Supp. Fig. 4 reads the expression matrices) and saves a PDF and a PNG to `$PLOT_DIR`.

| File | Figure |
|------|--------|
| `Fig1_SuppFig1-2.R` | Fig. 1 (best track per model), Supp. Figs 1–2 (every CAGE track; Spearman and Pearson) |
| `Fig2.R` | Fig. 2: ΔR² in LBP brain with five gene-mean baselines |
| `Fig3.R` | Fig. 3: personalized minus reference genome |
| `SuppFig3.R` | Supp. Fig. 3: each model minus the gene-mean baseline from its own training data |
| `SuppFig4.R` | Supp. Fig. 4: gene-mean baseline from each dataset scored against every other dataset |
| `SuppFig5.R` | Supp. Fig. 5: ΔR² in LBP blood, ROSMAP brain and COVID-19 blood |

Supp. Tables 2 and 3 are the summary tables written by `delta_r2/01_delta_r2_discovery.R` and `delta_r2/02_delta_r2_replication.R`.
