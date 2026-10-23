# figures

There is one script per figure. Each reads tables written by `cross_gene_correlation` or `delta_r2` and saves a PDF and a PNG to `$PLOT_DIR`.

| File | Figure |
|------|--------|
| `Fig1_SuppFig1-2.R` | Fig. 1 (best track per model), Supp. Figs 1–2 (every CAGE track; Spearman and Pearson) |
| `Fig2.R` | Fig. 2: each model minus the gene-mean baseline from its own training data |
| `Fig3.R` | Fig. 3: personalized minus reference genome |
| `SuppFig3.R` | Supp. Fig. 3: gene-mean baseline from each dataset scored against every other dataset |
| `SuppFig4.R` | Supp. Fig. 4: ΔR² in the replication cohorts |

Fig. 4 is drawn by `delta_r2/01_delta_r2_discovery.R`, and Supp. Figs 5–6 by `delta_r2/03_delta_r2_reference_genome.R`.
