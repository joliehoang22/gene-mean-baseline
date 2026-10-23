#!/usr/bin/env Rscript
## Fig. 2: per-sample paired difference between each model and the gene-mean
## baseline from its own training data (model minus baseline). Paired Wilcoxon,
## Bonferroni-adjusted. Reads cross_gene_correlation/01_score_all_cohorts.R output.
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
SRC <- file.path(BASE, "evaluation", "allcohorts", "arm_scores_allcohorts_qcexact.csv")
d   <- read.csv(SRC, stringsAsFactors = FALSE)

COHORT_LEV <- c("LBP-Brain","LBP-Blood","ROSMAP-Brain","COVID-Blood")
COHORT_LEV <- COHORT_LEV[COHORT_LEV %in% unique(d$cohort)]
fam_cols   <- c(Enformer = "#0072B2", Basenji2 = "#009E73", ExPecto = "#D55E00")

PAIRS <- list(
  list(model="Enformer", gm_brain="Gene-mean: FANTOM5 CAGE", gm_blood="Gene-mean: FANTOM5 CAGE", short="FANTOM5 CAGE"),
  list(model="Basenji2", gm_brain="Gene-mean: FANTOM5 CAGE", gm_blood="Gene-mean: FANTOM5 CAGE", short="FANTOM5 CAGE"),
  list(model="ExPecto",  gm_brain="Gene-mean: GTEx",         gm_blood="Gene-mean: GTEx Blood",   short="GTEx"))
PAIR_LEV <- c("Enformer  -  FANTOM5 CAGE", "Basenji2  -  FANTOM5 CAGE", "ExPecto  -  GTEx")

best_track <- function(x, mdl) {
  z <- x[x$model == mdl, ]; if (!nrow(z)) return(NA_character_)
  a <- aggregate(spearman ~ track, z, function(v) median(v, na.rm = TRUE))
  a$track[which.max(a$spearman)]
}
fmt_p <- function(p) ifelse(p < 1e-16, "p < 1e-16",
                     ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p)))

build <- function(metric) {
  diffs <- list(); stats <- list()
  for (co in COHORT_LEV) {
    x <- d[d$cohort == co, ]; is_blood <- grepl("Blood", co)
    for (pi in seq_along(PAIRS)) {
      p <- PAIRS[[pi]]
      bt <- best_track(x, p$model); if (is.na(bt)) next
      gmn <- if (is_blood) p$gm_blood else p$gm_brain
      m <- x[x$model == p$model & x$track == bt, c("sample", metric)]
      g <- x[x$track == gmn,                     c("sample", metric)]
      names(m)[2] <- "model"; names(g)[2] <- "genemean"
      j <- merge(m, g, by = "sample")
      j <- j[is.finite(j$model) & is.finite(j$genemean), ]
      if (nrow(j) < 5) next
      lab <- PAIR_LEV[pi]
      diffs[[paste(co,p$model)]] <- data.frame(cohort=co, pair=lab, model=p$model,
                                               diff=j$model - j$genemean, stringsAsFactors=FALSE)
      w <- wilcox.test(j$model, j$genemean, paired = TRUE)
      stats[[paste(co,p$model)]] <- data.frame(cohort=co, pair=lab, model=p$model,
        n=nrow(j), med=median(j$model - j$genemean), p=w$p.value, stringsAsFactors=FALSE)
    }
  }
  D <- do.call(rbind, diffs); S <- do.call(rbind, stats)
  if (is.null(D)) { cat("no pairings for", metric, "\n"); return(invisible()) }
  S$p_bonf <- pmin(1, S$p * nrow(S))
  S$txt <- sprintf("median %+.3f   %s", S$med, fmt_p(S$p_bonf))

  D$pair   <- factor(D$pair,   levels = rev(PAIR_LEV))
  S$pair   <- factor(S$pair,   levels = rev(PAIR_LEV))
  D$cohort <- factor(D$cohort, levels = COHORT_LEV)
  S$cohort <- factor(S$cohort, levels = COHORT_LEV)
  D$model  <- factor(D$model,  levels = names(fam_cols))
  S$model  <- factor(S$model,  levels = names(fam_cols))

  rng <- range(D$diff, na.rm = TRUE); pad <- diff(rng) * 0.55
  stat_word <- if (metric == "spearman") "Spearman" else "Pearson"
  ylab_expr <- if (metric == "spearman")
      expression("Difference in cross-gene Spearman " * rho * "  (model - gene mean)")
    else
      expression("Difference in cross-gene Pearson " * italic(r) * "  (model - gene mean)")

  p <- ggplot(D, aes(pair, diff, fill = model)) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_violin(scale = "width", width = 0.85, alpha = 0.55, colour = NA, trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.size = 0.2, alpha = 0.95,
                 colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
    geom_text(data = S, aes(pair, y = rng[2] + pad * 0.10, label = txt),
              hjust = 0, size = 2.8, colour = "grey25", inherit.aes = FALSE) +
    scale_fill_manual(values = fam_cols, drop = TRUE, name = NULL) +
    coord_flip(ylim = c(rng[1] - pad * 0.05, rng[2] + pad)) +
    facet_wrap(~ cohort, ncol = 2) +
    labs(title = sprintf("Does each model beat the gene mean built from its own training data? (%s)", stat_word),
         subtitle = "per-sample paired difference (model - matched gene mean); values below 0 mean the gene mean wins.\nEnformer and Basenji2 were trained on FANTOM5 CAGE; ExPecto on GTEx. Paired Wilcoxon, Bonferroni-adjusted.",
         x = NULL, y = ylab_expr) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top",
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

  fn <- sprintf("paired_diff_model_vs_training_genemean_qc_%s", metric)
  ggsave(file.path(PLOT_DIR, paste0(fn, ".pdf")), p, width = 13, height = 7)
  ggsave(file.path(PLOT_DIR, paste0(fn, ".png")), p, width = 13, height = 7, dpi = 300)
  cat("\n=== ", stat_word, " ===\n", sep = "")
  print(as.data.frame(S[, c("cohort","pair","n","med","p","p_bonf")]), row.names = FALSE, digits = 3)
  cat("saved: ", fn, ".pdf/.png\n", sep = "")
  invisible(S)
}

for (mt in c("spearman","pearson")) build(mt)
