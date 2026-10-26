#!/usr/bin/env Rscript
## Supp. Fig. 4: gene-mean baseline from each dataset scored against every other
## dataset, on the scored panel (bed, all cohorts and FANTOM5 symbols intersected).
## Diagonal cells are leave-one-out; FANTOM5 has no individual samples, so its column is a
## gene-mean-to-gene-mean correlation and its diagonal is blank.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUT_DIR <- file.path(BASE, "Null")
FAN     <- read.csv(FANTOM5_MEANS, stringsAsFactors = FALSE)

strip <- function(x) sub("\\..*", "", x)
prep <- function(rds, ord, bed, fcol, tag) {
  M <- readRDS(file.path(BASE, rds))[ord]
  for (n in names(M)) {
    rownames(M[[n]]) <- strip(rownames(M[[n]]))
    M[[n]] <- M[[n]][!duplicated(rownames(M[[n]])), , drop = FALSE]
  }
  b   <- read.table(file.path(BASE, bed), sep = "\t", stringsAsFactors = FALSE,
                    col.names = c("chrom", "tss", "strand", "gene_id", "gene_name"))
  sym <- setNames(b$gene_name, strip(b$gene_id))
  p   <- Reduce(intersect, lapply(M, rownames))
  p   <- intersect(p, strip(b$gene_id))
  p   <- p[sym[p] %in% FAN$gene]
  cat(sprintf("%s: %d genes\n", tag, length(p)))
  list(M = lapply(M, function(m) m[p, , drop = FALSE]), panel = p,
       fan = log2(FAN[[fcol]][match(sym[p], FAN$gene)] + 1))
}
pairwise <- function(D, method) {
  M <- D$M; gm <- lapply(M, rowMeans); gm$fantom <- D$fan
  src <- names(gm); tgt <- c(names(M), "fantom")
  R <- matrix(NA_real_, length(src), length(tgt), dimnames = list(src, tgt))
  loo <- function(X) { N <- ncol(X); rs <- rowSums(X)
    vapply(seq_len(N), function(i) cor((rs - X[, i]) / (N - 1), X[, i], method = method), numeric(1)) }
  for (s in src) for (t in tgt) {
    R[s, t] <-
      if (t == "fantom") { if (s == "fantom") NA_real_ else cor(gm[[s]], D$fan, method = method, use = "complete.obs") }
      else if (s == t)   median(loo(M[[t]]))
      else               median(apply(M[[t]], 2, function(col) cor(gm[[s]], col, method = method, use = "complete.obs")))
  }
  R
}

BR <- prep("Null/mats.rds",       c("lbp_liv", "lbp_pm", "msbb", "rosmap", "gtex"),
           "tss_hg38_mane_expressed_brain.bed", "fantom_brain_mean", "BRAIN")
BL <- prep("Null/mats_blood.rds", c("lbp_blood", "covid_blood", "gtex_blood"),
           "tss_hg38_mane_expressed_blood.bed", "fantom_blood_mean", "BLOOD")
brain_sp <- pairwise(BR, "spearman"); brain_pe <- pairwise(BR, "pearson")
blood_sp <- pairwise(BL, "spearman"); blood_pe <- pairwise(BL, "pearson")
for (nm in c("brain_sp", "brain_pe", "blood_sp", "blood_pe"))
  write.csv(get(nm), file.path(OUT_DIR, sprintf("transfer_%s_scoredpanel.csv", nm)))

brain_ord <- c("lbp_liv", "lbp_pm", "msbb", "rosmap", "gtex", "fantom")
blood_ord <- c("lbp_blood", "covid_blood", "gtex_blood", "fantom")
labs_nice <- c(lbp_liv = "LBP living", lbp_pm = "LBP PM", msbb = "MSBB", rosmap = "ROSMAP", gtex = "GTEx",
               lbp_blood = "LBP living", covid_blood = "COVID-19", gtex_blood = "GTEx", fantom = "FANTOM5 CAGE")
lims <- c(floor(min(brain_sp, brain_pe, blood_sp, blood_pe, na.rm = TRUE) * 10) / 10, 1)
cols <- c("#FFF7EC", "#FDBB84", "#D7301F")

make_hm <- function(M, ord, title) {
  df <- as.data.frame(as.table(M)); names(df) <- c("source", "target", "val")
  df <- df[df$source %in% ord & df$target %in% ord, ]
  df$source <- factor(labs_nice[as.character(df$source)], levels = rev(unique(labs_nice[ord])))
  df$target <- factor(labs_nice[as.character(df$target)], levels = unique(labs_nice[ord]))
  ggplot(df, aes(target, source, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(data = df[is.finite(df$val), ],
              aes(label = sprintf("%.2f", val)), colour = "grey15", size = 3.2, show.legend = FALSE) +
    scale_fill_gradientn(colours = cols, limits = lims, name = "Median Correlation", na.value = "grey92") +
    coord_fixed() + labs(title = title, x = "target", y = "source (gene mean)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))
}
combined <- (make_hm(brain_sp, brain_ord, "Brain: Cross-gene Spearman Correlation") |
             make_hm(brain_pe, brain_ord, "Brain: Cross-gene Pearson Correlation")) /
            (make_hm(blood_sp, blood_ord, "Blood: Cross-gene Spearman Correlation") |
             make_hm(blood_pe, blood_ord, "Blood: Cross-gene Pearson Correlation")) +
            plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

ggsave(file.path(PLOT_DIR, "transfer_heatmaps_2x2_scoredpanel.pdf"), combined, width = 14, height = 11.5)
ggsave(file.path(PLOT_DIR, "transfer_heatmaps_2x2_scoredpanel.png"), combined, width = 14, height = 11.5, dpi = 300)
cat("\nsaved transfer_heatmaps_2x2_scoredpanel.pdf/.png\n")
