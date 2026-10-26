#!/usr/bin/env Rscript
## Fig. 1 and Supp. Figs 1-2: cross-gene correlation of each model and each gene-mean
## baseline in all four cohorts (clean = best track per model; alltracks = every CAGE track).
## Both cohorts of a tissue are compared with the same gene-mean baselines:
##   brain (LBP-Brain, ROSMAP-Brain): LBP PM, MSBB, GTEx, FANTOM5 CAGE
##   blood (LBP-Blood, COVID-Blood):  GTEx Blood, FANTOM5 CAGE
## Panels: A LBP-Brain, B ROSMAP-Brain, C LBP-Blood, D COVID-Blood. The best Enformer and
## Basenji2 track is the one with the highest median; numbers right of each violin are means.
## Reads cross_gene_correlation/01_score_all_cohorts.R output; writes cross_gene_allcohorts_expanded_qc_*
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
SRC <- file.path(BASE, "evaluation", "allcohorts", "arm_scores_allcohorts_qcexact.csv")
stopifnot(file.exists(SRC))

COHORTS <- c("LBP-Brain", "ROSMAP-Brain", "LBP-Blood", "COVID-Blood")
REFS <- list(brain = c("Gene-mean: LBP PM", "Gene-mean: MSBB", "Gene-mean: GTEx", "Gene-mean: FANTOM5 CAGE"),
             blood = c("Gene-mean: GTEx Blood", "Gene-mean: FANTOM5 CAGE"))
tissue_of <- function(co) if (grepl("Brain", co)) "brain" else "blood"
FAMS     <- c("Gene-mean", "Enformer", "Basenji2", "ExPecto", "Xpresso")
fam_cols <- c("Gene-mean" = "grey45", "Enformer" = "#0072B2", "Basenji2" = "#009E73",
              "ExPecto" = "#D55E00", "Xpresso" = "#CC79A7")

family_of <- function(m) {
  m <- as.character(m)
  ifelse(grepl("Gene-mean", m), "Gene-mean",
  ifelse(grepl("Enformer", m), "Enformer",
  ifelse(grepl("Basenji", m),  "Basenji2",
  ifelse(grepl("ExPecto", m),  "ExPecto",
  ifelse(grepl("Xpresso", m),  "Xpresso", "Other")))))
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
  ismod <- !grepl("^Gene-mean", out) & !out %in% c("ExPecto", "Xpresso")
  out[ismod] <- title_case(out[ismod]); out
}

JOBS <- expand.grid(metric = c("spearman", "pearson"), tag = c("clean", "alltracks"), stringsAsFactors = FALSE)
JOBS$fn <- sprintf("cross_gene_allcohorts_expanded_qc_%s_%s", JOBS$tag, JOBS$metric)

## keep the four models and only the chosen gene means for each tissue
d <- read.csv(SRC, stringsAsFactors = FALSE)
d$family <- family_of(d$model)
d <- d[d$family %in% FAMS & d$cohort %in% COHORTS, ]
d$arm <- pretty_arm(short_arm(d$arm))
is_ref <- mapply(function(a, co) a %in% REFS[[tissue_of(co)]], d$arm, d$cohort)
d <- d[d$family != "Gene-mean" | is_ref, ]
for (co in COHORTS) {
  got  <- unique(d$arm[d$cohort == co & d$family == "Gene-mean"])
  miss <- setdiff(REFS[[tissue_of(co)]], got)
  if (length(miss)) stop(co, ": gene mean not found: ", paste(miss, collapse = ", "))
  cat(sprintf("  %-13s %4d samples | %s\n", co, length(unique(d$sample[d$cohort == co])),
              paste(sub("Gene-mean: ", "", REFS[[tissue_of(co)]]), collapse = ", ")))
}

axis_cfg <- list(
  spearman = list(lab = expression("Cross-gene Spearman " * rho * " (per sample)"), brks = seq(0, 1, 0.2)),
  pearson  = list(lab = "Cross-gene Pearson r (per sample)", brks = seq(0, 1, 0.2)))
LABEL_X <- 1.04; VIEW_HI <- 1.14

make_plot <- function(co, metric, best_track_only) {
  cfg <- axis_cfg[[metric]]
  x <- d[d$cohort == co, ]; x$value <- x[[metric]]; x <- x[is.finite(x$value), ]
  med  <- aggregate(value ~ arm + family, x, median)      # picks the best track
  mods <- med[med$family != "Gene-mean", ]
  if (best_track_only) {
    mods <- do.call(rbind, lapply(split(mods, mods$family), function(z) z[which.max(z$value), ]))
    x   <- x[x$family == "Gene-mean" | x$arm %in% mods$arm, ]
  }
  avg  <- aggregate(value ~ arm + family, x, mean)        # numbers shown
  mods <- mods[order(match(mods$family, FAMS), -mods$value), ]
  rows <- rev(c(REFS[[tissue_of(co)]], mods$arm))          # first row is drawn at the top
  x$arm <- factor(x$arm, levels = rows); avg$arm <- factor(avg$arm, levels = rows)
  x$family <- factor(x$family, levels = FAMS)
  ggplot(x, aes(arm, value, fill = family)) +
    geom_violin(scale = "width", width = 0.9, alpha = 0.5, colour = NA, trim = TRUE) +
    geom_boxplot(width = 0.16, outlier.size = 0.25, alpha = 0.95,
                 colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
    geom_text(data = avg, aes(arm, label = sprintf("%.2f", value)),
              y = LABEL_X, hjust = 0, size = 2.6, colour = "grey25", inherit.aes = FALSE) +
    scale_fill_manual(values = fam_cols, drop = FALSE, name = NULL) +
    scale_y_continuous(breaks = cfg$brks) +
    coord_flip(ylim = c(0, VIEW_HI)) +
    labs(title = co, x = NULL, y = cfg$lab) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top",
          plot.title = element_text(face = "bold"))
}

for (i in seq_len(nrow(JOBS))) {
  j  <- JOBS[i, ]
  ps <- lapply(seq_along(COHORTS), function(k)
          make_plot(COHORTS[k], j$metric, j$tag == "clean") + labs(tag = LETTERS[k]))
  ## row 1: legend across both columns; rows 2-3: A B / C D
  g  <- wrap_plots(c(list(guide_area()), ps), design = "AA\nBC\nDE",
                   heights = c(0.06, 1, 1), guides = "collect") &
        theme(legend.position = "top", legend.direction = "horizontal")
  h  <- if (j$tag == "clean") 9 else 13
  ggsave(file.path(PLOT_DIR, paste0(j$fn, ".pdf")), g, width = 17, height = h, limitsize = FALSE)
  ggsave(file.path(PLOT_DIR, paste0(j$fn, ".png")), g, width = 17, height = h, dpi = 300, limitsize = FALSE)
  cat("saved:", j$fn, "\n")
}

## Fig. 1 values (clean, Spearman): the mean shown in the figure, and the median
x <- d[is.finite(d$spearman), ]
tab <- merge(aggregate(cbind(mean = spearman) ~ cohort + arm + family, x, mean),
             aggregate(cbind(median = spearman) ~ cohort + arm + family, x, median))
tab <- tab[order(match(tab$cohort, COHORTS), match(tab$family, FAMS), -tab$median), ]
tab$mean <- round(tab$mean, 3); tab$median <- round(tab$median, 3)
print(tab[tab$family == "Gene-mean" | !duplicated(tab[, c("cohort", "family")]), c("cohort", "arm", "mean", "median")], row.names = FALSE)
