#!/usr/bin/env Rscript
## Fig. 1 and Supp. Figs 1-2: cross-gene correlation for every arm in all four
## cohorts (clean = best track per model; alltracks = every CAGE track).
## Reads cross_gene_correlation/01_score_all_cohorts.R output; writes cross_gene_allcohorts_expanded_*
suppressPackageStartupMessages(library(ggplot2))

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
SRC <- file.path(BASE, "evaluation", "allcohorts", "arm_scores_allcohorts_qcexact.csv")
stopifnot(file.exists(SRC))

COHORTS <- c("LBP-Brain","LBP-Blood","ROSMAP-Brain","COVID-Blood")
fam_cols <- c("Gene-mean" = "grey45", "Enformer" = "#0072B2",
              "Basenji2" = "#009E73", "ExPecto" = "#D55E00",
              "Xpresso" = "#CC79A7", "Borzoi" = "#A11D33")
family_of <- function(m) {
  m <- as.character(m)
  ifelse(grepl("Gene-mean", m), "Gene-mean",
  ifelse(grepl("Enformer", m), "Enformer",
  ifelse(grepl("Basenji", m),  "Basenji2",
  ifelse(grepl("ExPecto", m),  "ExPecto",
  ifelse(grepl("Xpresso", m),  "Xpresso",
  ifelse(grepl("Borzoi", m),   "Borzoi", "Other"))))))
}
short_arm <- function(a) {
  a <- gsub("CAGE[:_]\\s*", "", a); a <- gsub("_", " ", a)
  a <- gsub("(,\\s*)+", ", ", a)
  a <- gsub("\\s*,\\s*$", "", a); a <- gsub("\\s+", " ", a); trimws(a)
}
title_case <- function(s) gsub("(^|[ ])([a-z])", "\\1\\U\\2", s, perl = TRUE)
pretty_arm <- function(a) {
  out <- a
  out <- sub("^(ExPecto|Xpresso): .*$", "\\1", out)
  ismod <- !grepl("^Gene-mean", out) & !out %in% c("ExPecto","Xpresso")
  out[ismod] <- title_case(out[ismod]); out
}

d <- read.csv(SRC, stringsAsFactors = FALSE)
d$family <- factor(family_of(d$model), levels = names(fam_cols))
d$arm    <- pretty_arm(short_arm(d$arm))
d <- d[!is.na(d$family), ]
d$family <- droplevels(d$family)   # hide families with no data (e.g. Borzoi)
COHORTS <- COHORTS[COHORTS %in% unique(d$cohort)]
cat("cohorts:", paste(COHORTS, collapse=", "), "\n")
for (co in COHORTS)
  cat(sprintf("  %-14s N = %d samples\n", co, length(unique(d$sample[d$cohort==co]))))

axis_cfg <- list(
  spearman = list(lab = expression("Cross-gene Spearman " * rho * " (per sample)"),
                  lo = 0.0, brks = seq(0.0, 1.0, 0.2)),
  pearson  = list(lab = "Cross-gene Pearson r (per sample)",
                  lo = 0.0, brks = seq(0.0, 1.0, 0.2)))
LABEL_X <- 1.04; VIEW_HI <- 1.14

make_plot <- function(df, co, metric, best_track_only = FALSE) {
  cfg <- axis_cfg[[metric]]
  x <- df[df$cohort == co, ]; x$value <- x[[metric]]; x <- x[is.finite(x$value), ]
  if (!nrow(x)) return(NULL)
  med <- aggregate(value ~ arm + family, x, median)
  if (best_track_only) {
    gm   <- med[med$family == "Gene-mean", ]
    mods <- med[med$family != "Gene-mean", ]
    if (nrow(mods)) {
      best <- do.call(rbind, by(mods, droplevels(mods$family), function(z) z[which.max(z$value), ]))
      keep <- c(gm$arm, best$arm)
    } else keep <- gm$arm
    x <- x[x$arm %in% keep, ]; med <- med[med$arm %in% keep, ]
  }
  x$arm <- factor(x$arm, levels = med$arm[order(med$value)])
  ceil  <- med$value[grepl("\\(LOO\\)", med$arm)][1]      # this cohort's own LOO gene-mean
  n <- length(unique(x$sample))
  p <- ggplot(x, aes(arm, value, fill = family))
  if (is.finite(ceil)) p <- p + geom_hline(yintercept = ceil, linetype = "dashed",
                                            colour = "grey55", linewidth = 0.4)
  p + geom_violin(scale = "width", width = 0.9, alpha = 0.5, colour = NA, trim = TRUE) +
    geom_boxplot(width = 0.16, outlier.size = 0.25, alpha = 0.95,
                 colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
    geom_text(data = med, aes(arm, label = sprintf("%.2f", value)),
              y = LABEL_X, hjust = 0, size = 2.6, colour = "grey25", inherit.aes = FALSE) +
    scale_fill_manual(values = fam_cols, drop = FALSE, name = NULL) +
    scale_y_continuous(breaks = cfg$brks) +
    coord_flip(ylim = c(cfg$lo, VIEW_HI)) +
    labs(title = co, x = NULL, y = cfg$lab,
         subtitle = sprintf("%d samples; dashed = within-cohort gene-mean (LOO)", n)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top",
          plot.title = element_text(face = "bold"))
}

save_combined <- function(metric, best, tag) {
  ps <- Filter(Negate(is.null), lapply(COHORTS, function(co) make_plot(d, co, metric, best)))
  if (!length(ps)) { cat("nothing to plot:", metric, tag, "\n"); return(invisible()) }
  fn <- sprintf("cross_gene_allcohorts_expanded_qc_%s_%s", tag, metric)
  h  <- if (tag == "clean") 9 else 13
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    g <- Reduce(`+`, ps) + plot_layout(ncol = 2, guides = "collect") &
         theme(legend.position = "top")
    ggsave(file.path(PLOT_DIR, paste0(fn, ".pdf")), g, width = 17, height = h, limitsize = FALSE)
    ggsave(file.path(PLOT_DIR, paste0(fn, ".png")), g, width = 17, height = h, dpi = 300, limitsize = FALSE)
  } else {
    for (i in seq_along(ps))
      ggsave(file.path(PLOT_DIR, sprintf("%s_%s.pdf", fn, COHORTS[i])), ps[[i]], width = 13, height = h/2)
  }
  cat("saved:", fn, "\n")
}

for (m in c("spearman", "pearson"))
  for (bt in list(list(TRUE, "clean"), list(FALSE, "alltracks")))
    save_combined(m, bt[[1]], bt[[2]])
