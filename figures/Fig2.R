#!/usr/bin/env Rscript
## Fig. 2: ΔR² in LBP brain for each model and each gene-mean baseline. Rows are the five
## baselines, coloured by how their RNA was prepared compared with LBP brain (Ribo-Zero);
## columns are the models, with headers shaded by training data (CAGE darker, RNA-seq lighter).
## Replots the per-sample table written by delta_r2/01_delta_r2_discovery.R (no regressions refit).
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(grid)})

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUT_DIR <- file.path(BASE, "evaluation", "complementarity")
stem    <- file.path(PLOT_DIR, "delta_rsquared_LBPbrain_5refs")

res <- read.csv(file.path(OUT_DIR, "dR2_external_genemeans_shared.csv"), stringsAsFactors = FALSE)
res <- res[res$tissue == "Brain", ]

MODELS <- c("Enformer", "Basenji2", "ExPecto", "Xpresso")
TRAIN  <- c(Enformer = "CAGE", Basenji2 = "CAGE", ExPecto = "RNA-seq", Xpresso = "RNA-seq")
REFS   <- c("GTEx", "ROSMAP", "MSBB", "LBP PM", "FANTOM5 CAGE")
PREP   <- setNames(c("poly-A", "poly-A", "Ribo-Zero", "Ribo-Zero", "CAGE"), REFS)
MATCH  <- c(`Ribo-Zero` = "Same library prep (Ribo-Zero)", `poly-A` = "Different library prep (poly-A)",
            CAGE = "Different assay (CAGE)")
MATCH_COLS <- setNames(c("#0072B2", "#D55E00", "#009E73"), MATCH)
stopifnot(all(res$gene_mean %in% REFS), all(res$model %in% MODELS))
cat(sprintf("%d samples x %d models x %d references\n", length(unique(res$sample)), length(MODELS), length(REFS)))

summ <- res %>% group_by(model, gene_mean) %>%
  summarise(dR2 = median(dR2), pct_added = 100 * dR2, n = n(), .groups = "drop")
print(as.data.frame(tidyr::pivot_wider(summ[, c("gene_mean", "model", "pct_added")],
                                       names_from = model, values_from = pct_added)), digits = 3)

lab_row <- function(g) sprintf("%s\n(%s)", g, PREP[g])
lab_col <- function(m) sprintf("%s (trained on %s)", m, TRAIN[m])
dress <- function(d) { d$row <- factor(lab_row(d$gene_mean), levels = lab_row(REFS))
                       d$col <- factor(lab_col(d$model),     levels = lab_col(MODELS))
                       d$match <- factor(MATCH[PREP[d$gene_mean]], levels = MATCH); d }
res <- dress(res); summ <- dress(summ)

p <- ggplot(res, aes(dR2, fill = match)) +
  geom_histogram(bins = 45, alpha = 0.8, colour = NA) +
  geom_vline(data = summ, aes(xintercept = dR2, colour = match), linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = summ, aes(x = Inf, y = Inf, label = sprintf("median %.2f%%", pct_added), colour = match),
            hjust = 1.08, vjust = 1.7, size = 3.2, show.legend = FALSE) +
  scale_fill_manual(values = MATCH_COLS, drop = FALSE, name = "Reference vs. study cohort (LBP living brain, Ribo-Zero):") +
  scale_colour_manual(values = MATCH_COLS, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(row ~ col, scales = "free_y") +
  labs(x = expression("Additional variance explained (" * Delta * italic(R)^2 * ")"), y = "samples") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour = "grey40"))

## shade the model headers by training data (CAGE darker, RNA-seq lighter); reference labels white
recolor <- function(gr, fill) {
  if (inherits(gr, "rect")) { gr$gp$fill <- fill; return(gr) }
  if (inherits(gr, "gtable")) { for (j in seq_along(gr$grobs)) gr$grobs[[j]] <- recolor(gr$grobs[[j]], fill); return(gr) }
  if (inherits(gr, "gTree"))  { for (j in seq_along(gr$children)) gr$children[[j]] <- recolor(gr$children[[j]], fill) }
  gr
}
g   <- ggplotGrob(p)
top <- which(grepl("^strip-t", g$layout$name)); top <- top[order(g$layout$l[top])]
for (k in seq_along(top)) g$grobs[[top[k]]] <- recolor(g$grobs[[top[k]]], if (TRAIN[MODELS[k]] == "CAGE") "grey70" else "grey95")
for (k in which(grepl("^strip-r", g$layout$name))) g$grobs[[k]] <- recolor(g$grobs[[k]], "white")

pdf(paste0(stem, ".pdf"), width = 14, height = 11); grid.draw(g); invisible(dev.off())
png(paste0(stem, ".png"), width = 14, height = 11, units = "in", res = 300); grid.draw(g); invisible(dev.off())
cat("\nsaved", paste0(stem, ".{pdf,png}"), "\n")
