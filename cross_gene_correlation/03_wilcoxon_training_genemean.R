#!/usr/bin/env Rscript
## Paired Wilcoxon tests behind Supp. Fig. 3: each model against the gene-mean baseline
## built from its own training data (FANTOM5 CAGE for Enformer and Basenji2,
## GTEx for ExPecto), per cohort. Writes the test table for the supplement.
suppressPackageStartupMessages({library(dplyr)})
source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
SRC  <- file.path(BASE, "evaluation", "allcohorts", "arm_scores_allcohorts_qcexact.csv")
d    <- read.csv(SRC, stringsAsFactors = FALSE)

## the gene-mean baseline built from each model's training data
PAIRS <- list(
  list(model = "Enformer", gm_brain = "Gene-mean: FANTOM5 CAGE", gm_blood = "Gene-mean: FANTOM5 CAGE",
       why = "Enformer trained on FANTOM5 CAGE"),
  list(model = "Basenji2", gm_brain = "Gene-mean: FANTOM5 CAGE", gm_blood = "Gene-mean: FANTOM5 CAGE",
       why = "Basenji2 trained on FANTOM5 CAGE"),
  list(model = "ExPecto",  gm_brain = "Gene-mean: GTEx",         gm_blood = "Gene-mean: GTEx Blood",
       why = "ExPecto trained on GTEx RNA-seq"))

best_track <- function(x, mdl) {
  z <- x[x$model == mdl, ]
  if (!nrow(z)) return(NA_character_)
  a <- aggregate(spearman ~ track, z, function(v) median(v, na.rm = TRUE))
  a$track[which.max(a$spearman)]
}

res <- list()
for (co in unique(d$cohort)) {
  x <- d[d$cohort == co, ]
  is_blood <- grepl("Blood", co)
  for (p in PAIRS) {
    bt <- best_track(x, p$model); if (is.na(bt)) next
    gmn <- if (is_blood) p$gm_blood else p$gm_brain
    for (metric in c("spearman","pearson")) {
      m  <- x[x$model == p$model & x$track == bt, c("sample", metric)]
      g  <- x[x$track == gmn,                     c("sample", metric)]
      names(m)[2] <- "model"; names(g)[2] <- "genemean"
      j  <- merge(m, g, by = "sample")
      j  <- j[is.finite(j$model) & is.finite(j$genemean), ]
      if (nrow(j) < 5) next
      w  <- wilcox.test(j$model, j$genemean, paired = TRUE)
      res[[paste(co,p$model,metric)]] <- data.frame(
        cohort = co, model = p$model, model_track = bt, gene_mean = gmn,
        metric = metric, n = nrow(j),
        median_model = median(j$model), median_genemean = median(j$genemean),
        median_diff = median(j$model - j$genemean),
        pct_model_higher = 100 * mean(j$model > j$genemean),
        W = unname(w$statistic), p_value = w$p.value,
        rationale = p$why, stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, res); rownames(out) <- NULL
out$p_bonferroni <- pmin(1, out$p_value * sum(out$metric == "spearman"))
out <- out[order(out$metric, out$cohort, out$model), ]

cat("\n=== paired Wilcoxon: model vs its training-data gene mean ===\n")
print(out[out$metric == "spearman",
          c("cohort","model","gene_mean","n","median_model","median_genemean",
            "median_diff","pct_model_higher","p_value","p_bonferroni")],
      row.names = FALSE, digits = 3)

write.csv(out, file.path(BASE, "evaluation", "allcohorts",
                         "wilcoxon_model_vs_training_genemean_qc.csv"), row.names = FALSE)
