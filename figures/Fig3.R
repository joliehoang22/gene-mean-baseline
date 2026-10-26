#!/usr/bin/env Rscript
## Fig. 3: paired personalized-vs-reference difference, 4 models x 4 cohorts.
## Reads cross_gene_correlation/02_score_personal_vs_reference.R output (has a `mode` column).
suppressPackageStartupMessages({library(ggplot2)})
source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
SRC <- file.path(BASE, "evaluation", "allcohorts", "arm_scores_allcohorts_withref.csv")
stopifnot(file.exists(SRC))

d <- read.csv(SRC, stringsAsFactors = FALSE)
d <- d[d$mode %in% c("personalized","reference"), ]
cat("cohorts:", paste(unique(d$cohort), collapse=", "), "\n")
cat("models :", paste(unique(d$model),  collapse=", "), "\n\n")

COHORT_LEV <- c("LBP-Brain","LBP-Blood","ROSMAP-Brain","COVID-Blood")
COHORT_LEV <- COHORT_LEV[COHORT_LEV %in% unique(d$cohort)]
MODEL_LEV  <- c("Enformer","Basenji2","ExPecto","Xpresso")
MODEL_LEV  <- MODEL_LEV[MODEL_LEV %in% unique(d$model)]
fam_cols   <- c(Enformer="#0072B2", Basenji2="#009E73", ExPecto="#D55E00", Xpresso="#CC79A7")

best_track <- function(x, mdl) {
  z <- x[x$model == mdl & x$mode == "personalized", ]
  if (!nrow(z)) return(NA_character_)
  a <- aggregate(spearman ~ track, z, function(v) median(v, na.rm = TRUE))
  a$track[which.max(a$spearman)]
}
fmt_p <- function(p) ifelse(is.na(p), "p = NA",
                    ifelse(p < 1e-16, "p < 1e-16",
                    ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p))))

build <- function(metric) {
  diffs <- list(); stats <- list()
  for (co in COHORT_LEV) {
    x <- d[d$cohort == co, ]
    for (mdl in MODEL_LEV) {
      bt <- best_track(x, mdl)
      if (is.na(bt)) { cat("  skip (no personalized):", co, mdl, "\n"); next }
      p <- x[x$model==mdl & x$track==bt & x$mode=="personalized", c("sample", metric)]
      r <- x[x$model==mdl & x$track==bt & x$mode=="reference",    c("sample", metric)]
      if (!nrow(r)) { cat("  skip (no reference):", co, mdl, bt, "\n"); next }
      names(p)[2] <- "personal"; names(r)[2] <- "reference"
      j <- merge(p, r, by = "sample")
      j <- j[is.finite(j$personal) & is.finite(j$reference), ]
      if (nrow(j) < 5) { cat("  skip (n<5):", co, mdl, "\n"); next }
      dv <- j$personal - j$reference
      w  <- tryCatch(suppressWarnings(wilcox.test(j$personal, j$reference, paired=TRUE, conf.int=TRUE)),
                     error = function(e) NULL)
      ci <- if (is.null(w) || is.null(w$conf.int)) c(NA_real_, NA_real_) else as.numeric(w$conf.int)
      diffs[[paste(co,mdl)]] <- data.frame(cohort=co, model=mdl, diff=dv, stringsAsFactors=FALSE)
      stats[[paste(co,mdl)]] <- data.frame(cohort=co, model=mdl, track=bt, metric=metric,
        n=nrow(j), med_personal=median(j$personal), med_reference=median(j$reference),
        med_diff=median(dv), max_abs_diff=max(abs(dv)), pct_personal_higher=100*mean(dv > 0),
        ci_lo=ci[1], ci_hi=ci[2], p_value=if (is.null(w)) NA_real_ else w$p.value,
        stringsAsFactors=FALSE)
    }
  }
  D <- do.call(rbind, diffs); S <- do.call(rbind, stats)
  if (is.null(D)) { cat("no pairings for", metric, "\n"); return(NULL) }
  S$p_bonf <- pmin(1, S$p_value * nrow(S))
  S$txt <- sprintf("median %+.5f   %s", S$med_diff, fmt_p(S$p_bonf))

  D$model  <- factor(D$model,  levels = rev(MODEL_LEV))
  S$model  <- factor(S$model,  levels = rev(MODEL_LEV))
  D$cohort <- factor(D$cohort, levels = COHORT_LEV)
  S$cohort <- factor(S$cohort, levels = COHORT_LEV)
  D$fill   <- factor(as.character(D$model), levels = names(fam_cols))

  rng <- range(D$diff, na.rm = TRUE); span <- max(diff(rng), 1e-6); pad <- span * 1.15
  stat_word <- if (metric == "spearman") "Spearman" else "Pearson"
  ylab_expr <- if (metric == "spearman")
      expression("Difference in cross-gene Spearman " * rho * "  (personalized - reference)")
    else
      expression("Difference in cross-gene Pearson " * italic(r) * "  (personalized - reference)")

  p <- ggplot(D, aes(model, diff, fill = fill)) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_violin(scale = "width", width = 0.85, alpha = 0.55, colour = NA, trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.size = 0.2, alpha = 0.95,
                 colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
    geom_text(data = S, aes(model, y = 0.0098, label = txt),
              hjust = 1, size = 2.8, colour = "grey25", inherit.aes = FALSE) +
    scale_fill_manual(values = fam_cols, drop = TRUE, name = NULL) +
    scale_y_continuous(breaks = seq(-0.01, 0.01, 0.005)) +
    coord_flip(ylim = c(-0.01, 0.01)) +
    facet_wrap(~ cohort, ncol = 2) +
    labs(x = NULL, y = ylab_expr) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top",
          strip.text = element_text(face = "bold"))

  fn <- sprintf("paired_diff_personal_vs_reference_allcohorts_%s", metric)
  ggsave(file.path(PLOT_DIR, paste0(fn, ".pdf")), p, width = 15, height = 8)
  ggsave(file.path(PLOT_DIR, paste0(fn, ".png")), p, width = 15, height = 8, dpi = 300)
  cat("\n=== ", stat_word, " ===\n", sep = "")
  print(as.data.frame(S[, c("cohort","model","n","med_diff","max_abs_diff",
                            "pct_personal_higher","p_value","p_bonf")]),
        row.names = FALSE, digits = 4)
  cat("saved: ", fn, ".pdf/.png\n", sep = "")
  S
}

allS <- do.call(rbind, Filter(Negate(is.null), lapply(c("spearman","pearson"), build)))
if (!is.null(allS)) {
  f <- file.path(BASE, "evaluation", "wilcoxon_personal_vs_reference_allcohorts.csv")
  write.csv(allS[, setdiff(names(allS), "txt")], f, row.names = FALSE)
  cat("\nwrote", f, "\n")
}
