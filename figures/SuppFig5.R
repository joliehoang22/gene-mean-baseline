#!/usr/bin/env Rscript
## Supp. Fig. 5: ΔR² in LBP blood, ROSMAP brain and COVID-19 blood. Replots the per-sample
## tables written by delta_r2/01_delta_r2_discovery.R (LBP blood rows) and
## delta_r2/02_delta_r2_replication.R (no regressions refit). Gene-mean baselines:
##   LBP blood:    GTEx, COVID-19, FANTOM5
##   ROSMAP brain: LBP PM, MSBB, GTEx, FANTOM5
##   COVID blood:  GTEx, LBP Blood, FANTOM5
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUT_DIR <- file.path(BASE, "evaluation", "complementarity")
stem    <- file.path(PLOT_DIR, "delta_rsquared_LBPblood_ROSMAP_COVID")

disc <- read.csv(file.path(OUT_DIR, "dR2_external_genemeans_shared.csv"),         stringsAsFactors = FALSE)
disc <- disc[disc$tissue == "Blood", ]; disc$cohort <- "LBP-Blood"
repl <- read.csv(file.path(OUT_DIR, "dR2_external_genemeans_replication_qc.csv"), stringsAsFactors = FALSE)
repl$cohort <- c(`ROSMAP Brain` = "ROSMAP-Brain", `COVID Blood` = "COVID-Blood")[repl$cohort]
cols <- c("cohort", "model", "gene_mean", "sample", "dR2")
res  <- rbind(disc[, cols], repl[, cols])
res$gene_mean[res$gene_mean == "FANTOM5 CAGE"] <- "FANTOM5"

MODELS  <- c("Enformer", "Basenji2", "ExPecto", "Xpresso")
COH_LEV <- c("LBP-Blood", "ROSMAP-Brain", "COVID-Blood")
GM_LEV  <- c("LBP PM", "MSBB", "GTEx", "LBP Blood", "COVID-19", "FANTOM5")   # also the order of the corner labels
gm_cols <- c(`LBP PM` = "#1971C2", MSBB = "#F59F00", GTEx = "#7048E8",
             `LBP Blood` = "#E03131", `COVID-19` = "#0CA678", FANTOM5 = "#868E96")
stopifnot(!anyNA(res$cohort), all(res$gene_mean %in% GM_LEV), all(res$model %in% MODELS))

summ <- res %>% group_by(cohort, model, gene_mean) %>%
  summarise(dR2 = median(dR2), pct_added = 100 * dR2, n = n(), .groups = "drop")
for (co in COH_LEV)
  cat(sprintf("  %-13s %3d samples | %s\n", co, length(unique(res$sample[res$cohort == co])),
              paste(intersect(GM_LEV, unique(res$gene_mean[res$cohort == co])), collapse = ", ")))
cat("\n=== median dR2 (%) ===\n")
print(as.data.frame(tidyr::pivot_wider(summ[, c("cohort", "gene_mean", "model", "pct_added")],
                                       names_from = model, values_from = pct_added)), digits = 3)

res$model     <- factor(res$model,  levels = MODELS)
res$cohort    <- factor(res$cohort, levels = COH_LEV)
res$gene_mean <- factor(res$gene_mean, levels = GM_LEV)
med <- summ
med$model     <- factor(med$model,  levels = MODELS)
med$cohort    <- factor(med$cohort, levels = COH_LEV)
med$gene_mean <- factor(med$gene_mean, levels = GM_LEV)
med <- med %>% group_by(cohort, model) %>% arrange(gene_mean, .by_group = TRUE) %>%
  mutate(vjust = 1.7 + 1.6 * (row_number() - 1)) %>% ungroup()

p <- ggplot(res, aes(dR2, fill = gene_mean)) +
  geom_histogram(bins = 45, position = "identity", alpha = 0.55, colour = NA) +
  geom_vline(data = med, aes(xintercept = dR2, colour = gene_mean),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = med,
            aes(x = Inf, y = Inf, label = sprintf("%s  %.2f%%", as.character(gene_mean), pct_added),
                colour = gene_mean),
            hjust = 1.08, vjust = med$vjust, size = 3, show.legend = FALSE) +
  scale_fill_manual(values = gm_cols, drop = TRUE, name = NULL) +
  scale_colour_manual(values = gm_cols, drop = TRUE) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(cohort ~ model, scales = "free_y") +
  labs(x = expression("Additional variance explained (" * Delta * italic(R)^2 * ")"), y = "samples") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(paste0(stem, ".pdf"), p, width = 13, height = 9.5)
ggsave(paste0(stem, ".png"), p, width = 13, height = 9.5, dpi = 300)
cat("\nsaved", paste0(stem, ".{pdf,png}"), "\n")
