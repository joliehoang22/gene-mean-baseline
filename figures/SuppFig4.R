#!/usr/bin/env Rscript
## Supp. Fig. 4: replication dR2 histograms, replotted from the table written by
## delta_r2/02_delta_r2_replication.R (no parquet re-read).
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUT_DIR  <- file.path(BASE, "evaluation", "complementarity")

res  <- read.csv(file.path(OUT_DIR, "dR2_external_genemeans_replication_qc.csv"),         stringsAsFactors = FALSE)
summ <- read.csv(file.path(OUT_DIR, "dR2_external_genemeans_replication_qc_summary.csv"), stringsAsFactors = FALSE)
for (cc in unique(res$cohort))
  cat(sprintf("%-14s n = %d samples\n", cc, length(unique(res$sample[res$cohort == cc]))))

res$model  <- factor(res$model,  levels = c("Enformer","Basenji2","ExPecto","Xpresso"))
res$cohort <- factor(res$cohort, levels = c("ROSMAP Brain","COVID Blood"))
GM_LEV  <- c("GTEx","MSBB","LBP (PM)","GTEx Blood","LBP Blood")
gm_cols <- c(GTEx = "#7048E8", MSBB = "#F59F00", `LBP (PM)` = "#0CA678",
             `GTEx Blood` = "#7048E8", `LBP Blood` = "#E03131")
res$gene_mean <- factor(res$gene_mean, levels = GM_LEV)

med <- summ
med$model     <- factor(med$model,     levels = levels(res$model))
med$cohort    <- factor(med$cohort,    levels = levels(res$cohort))
med$gene_mean <- factor(med$gene_mean, levels = GM_LEV)
med <- med %>% group_by(cohort, model) %>% arrange(gene_mean, .by_group = TRUE) %>%
  mutate(vjust = 1.7 + 1.6 * (row_number() - 1)) %>% ungroup()

p <- ggplot(res, aes(dR2, fill = gene_mean)) +
  geom_histogram(bins = 45, position = "identity", alpha = 0.55, colour = NA) +
  geom_vline(data = med, aes(xintercept = dR2, colour = gene_mean),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = med,
            aes(x = Inf, y = Inf,
                label = sprintf("%s  %.2f%%", as.character(gene_mean), pct_added),
                colour = gene_mean),
            hjust = 1.08, vjust = med$vjust, size = 3, show.legend = FALSE) +
  scale_fill_manual(values = gm_cols, drop = TRUE, name = NULL) +
  scale_colour_manual(values = gm_cols, drop = TRUE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(cohort ~ model, scales = "free_y") +
  labs(title = "Extra cross-gene variance explained by adding the model to an external gene mean",
       subtitle = expression("replication cohorts, shared gene panel; one value per sample; dashed lines and corner labels are the median " * Delta * italic(R)^2),
       x = expression("Additional variance explained (" * Delta * italic(R)^2 * ")"),
       y = "samples") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOT_DIR, "delta_rsquared_external_genemeans_replication_qc.pdf"), p, width = 13, height = 7)
ggsave(file.path(PLOT_DIR, "delta_rsquared_external_genemeans_replication_qc.png"), p, width = 13, height = 7, dpi = 300)
cat("saved delta_rsquared_external_genemeans_replication_qc.pdf/.png\n")
