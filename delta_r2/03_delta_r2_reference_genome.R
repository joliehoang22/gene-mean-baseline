## ============================================================================
## Fig. 4 and Supp. Fig. 4 repeated with reference-genome predictions
## (Supp. Figs 5-6). Samples, gene panel, gene means, regression and track
## choice are the same as in 01/02. The only change is that every sample gets
## the reference-genome prediction instead of its personalized one.
##   null <- lm(observed ~ gene_mean);  full <- lm(observed ~ dl + gene_mean)
## Reads the per-sample tables from 01 and 02 for the exact sample lists.
## ============================================================================
suppressPackageStartupMessages({library(arrow); library(dplyr); library(ggplot2)})
arrow::set_cpu_count(1); arrow::set_io_thread_count(1)

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUT_DIR <- file.path(BASE, "evaluation", "complementarity")

z <- function(x) (x - mean(x)) / sd(x)
P <- function(...) file.path(BASE, "predictions", ...)
MODELS <- c("Enformer", "Basenji2", "ExPecto", "Xpresso")
LOG    <- c(Enformer = TRUE, Basenji2 = TRUE, ExPecto = FALSE, Xpresso = FALSE)

## personalized folders are only opened for their column names (candidate tracks)
PERS <- list(
  brain = list(Enformer = P("enformer_expanded_personal_brain"), Basenji2 = P("basenji2_expanded_personal_brain"),
               ExPecto  = c(P("expecto_expanded_personal_brain_gpu"), P("expecto_expanded_personal_brain")),
               Xpresso  = P("xpresso_expanded_personal_brain")),
  blood = list(Enformer = P("enformer_expanded_personal_blood"), Basenji2 = P("basenji2_expanded_personal_blood"),
               ExPecto  = c(P("expecto_expanded_personal_blood"), P("expecto_expanded_personal_blood_gpu")),
               Xpresso  = P("xpresso_expanded_personal_blood")))
REF <- function(model, tis) P(switch(model,
  Enformer = sprintf("enformer_reference_expanded_%s/reference_genome.parquet", tis),
  Basenji2 = sprintf("basenji2_reference_expanded_%s/reference_genome.parquet", tis),
  ExPecto  = sprintf("expecto_reference_expanded_%s/reference_genome.parquet",  tis),
  Xpresso  = sprintf("xpresso_reference_expanded_%s/seqonly.parquet",           tis)))
TARGETS   <- read.delim(file.path(BASE, "targets_human.txt"), stringsAsFactors = FALSE)
TRACK_IDX <- list(brain = scan(file.path(BASE, "brain_tracks.txt"), quiet = TRUE),
                  blood = scan(file.path(BASE, "blood_tracks.txt"), quiet = TRUE))
FANTOM    <- read.csv(FANTOM5_MEANS, stringsAsFactors = FALSE)

disc <- read.csv(file.path(OUT_DIR, "dR2_external_genemeans_shared.csv"),         stringsAsFactors = FALSE)
repl <- read.csv(file.path(OUT_DIR, "dR2_external_genemeans_replication_qc.csv"), stringsAsFactors = FALSE)

STUDY <- list(
  `LBP Brain`    = list(fig = "fig4",  tis = "brain", bed = "tss_hg38_mane_expressed_brain.bed",
                        mats = "Null/mats.rds", obs = "lbp_liv",
                        gms = c(GTEx = "gtex", ROSMAP = "rosmap", MSBB = "msbb"),
                        samples = unique(disc$sample[disc$tissue == "Brain"])),
  `LBP Blood`    = list(fig = "fig4",  tis = "blood", bed = "tss_hg38_mane_expressed_blood.bed",
                        mats = "Null/mats_blood.rds", obs = "lbp_blood",
                        gms = c(GTEx = "gtex_blood", `COVID-19` = "covid_blood"),
                        samples = unique(disc$sample[disc$tissue == "Blood"])),
  `ROSMAP Brain` = list(fig = "supp4", tis = "brain", bed = "tss_hg38_mane_expressed_brain.bed",
                        mats = "Null/mats.rds", obs = "rosmap",
                        gms = c(GTEx = "gtex", MSBB = "msbb", `LBP (PM)` = "lbp_pm"),
                        samples = unique(repl$sample[repl$cohort == "ROSMAP Brain"])),
  `COVID Blood`  = list(fig = "supp4", tis = "blood", bed = "tss_hg38_mane_expressed_blood.bed",
                        mats = "Null/mats_blood.rds", obs = "covid_blood",
                        gms = c(`GTEx Blood` = "gtex_blood", `LBP Blood` = "lbp_blood"),
                        samples = unique(repl$sample[repl$cohort == "COVID Blood"])))

candidate_tracks <- function(dirs) {
  for (d in dirs) {
    f <- list.files(d, pattern = "\\.parquet$", full.names = TRUE)
    if (length(f)) return(setdiff(names(read_parquet(f[1])), "gene_id"))
  }
  stop("no personalized parquet found in ", paste(dirs, collapse = ", "))
}

## reference prediction for one model: same candidate tracks and best-track rule
## as the personalized run (Spearman vs the mean observed profile)
load_ref <- function(model, tis, panel, obs_profile) {
  f <- REF(model, tis)
  tracks <- candidate_tracks(PERS[[tis]][[model]])
  cols   <- tracks
  if (model == "Basenji2") {                      # reference columns are track_<index>
    tg   <- TARGETS[TARGETS$index %in% TRACK_IDX[[tis]], ]
    cols <- paste0("track_", tg$index[match(tracks, tg$description)])
  }
  p <- as.data.frame(read_parquet(f, col_select = all_of(c("gene_id", cols))))
  V <- matrix(sapply(cols, function(cl) p[[cl]][match(panel, p$gene_id)]), nrow = length(panel),
              dimnames = list(panel, tracks))
  bt <- if (ncol(V) > 1) tracks[which.max(apply(V, 2, function(v)
          cor(v, obs_profile, method = "spearman", use = "complete.obs")))] else tracks[1]
  v <- V[, bt]
  if (LOG[[model]]) v <- log2(v + 1)
  list(v = v, track = bt)
}

run_all <- function(obs, dl, gm) {
  t(sapply(colnames(dl), function(s) {
    ok <- is.finite(obs[, s]) & is.finite(dl[, s]) & is.finite(gm)
    Y <- z(obs[ok, s]); D <- z(dl[ok, s]); Mm <- z(gm[ok])
    n <- lm(Y ~ Mm); f <- lm(Y ~ D + Mm)
    c(r2_null = summary(n)$r.squared, r2_full = summary(f)$r.squared,
      dR2 = summary(f)$r.squared - summary(n)$r.squared,
      beta_gm_null = unname(coef(n)["Mm"]), beta_dl = unname(coef(f)["D"]),
      beta_gm = unname(coef(f)["Mm"]), n_genes = sum(ok))
  }))
}

out <- list()
for (st in names(STUDY)) {
  cfg <- STUDY[[st]]
  bed <- read.table(file.path(BASE, cfg$bed), sep = "\t",
                    col.names = c("chrom","tss","strand","gene_id","gene_name"), stringsAsFactors = FALSE)
  matsL <- readRDS(file.path(BASE, cfg$mats)); obs_full <- matsL[[cfg$obs]]
  panel <- intersect(bed$gene_id, rownames(obs_full))
  for (k in cfg$gms) panel <- intersect(panel, rownames(matsL[[k]]))
  panel <- panel[bed$gene_name[match(panel, bed$gene_id)] %in% FANTOM$gene]
  samp <- intersect(cfg$samples, colnames(obs_full))
  obs  <- obs_full[panel, samp, drop = FALSE]
  cat(sprintf("\n=== %s: %d genes x %d samples ===\n", st, length(panel), length(samp)))
  GM <- lapply(cfg$gms, function(k) rowMeans(matsL[[k]][panel, , drop = FALSE], na.rm = TRUE))
  for (mn in MODELS) {
    L <- load_ref(mn, cfg$tis, panel, rowMeans(obs))
    cat(sprintf("  %-9s reference track: %s\n", mn, L$track))
    dl <- matrix(L$v, nrow = length(panel), ncol = length(samp), dimnames = list(panel, samp))
    for (gname in names(GM)) {
      r <- as.data.frame(run_all(obs, dl, GM[[gname]]))
      r$sample <- rownames(r); rownames(r) <- NULL
      r$model <- mn; r$gene_mean <- gname; r$study <- st; r$fig <- cfg$fig; r$track <- L$track
      out[[paste(st, mn, gname)]] <- r
    }
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL

summ <- res %>% group_by(fig, study, model, gene_mean) %>%
  summarise(r2_null = median(r2_null), r2_full = median(r2_full), dR2 = median(dR2), pct_added = 100 * median(dR2),
            beta_gm_null = median(beta_gm_null), beta_dl = median(beta_dl), beta_gm = median(beta_gm),
            genes = median(n_genes), n = n(), .groups = "drop")
f4  <- res[res$fig == "fig4", ];   f4$tissue  <- sub("^LBP ", "", f4$study)
f4s <- summ[summ$fig == "fig4", ]; f4s$tissue <- sub("^LBP ", "", f4s$study)
s4  <- res[res$fig == "supp4", ];   s4$cohort  <- s4$study
s4s <- summ[summ$fig == "supp4", ]; s4s$cohort <- s4s$study
write.csv(f4,  file.path(OUT_DIR, "dR2_external_genemeans_shared_reference.csv"),                row.names = FALSE)
write.csv(f4s, file.path(OUT_DIR, "dR2_external_genemeans_shared_reference_summary.csv"),        row.names = FALSE)
write.csv(s4,  file.path(OUT_DIR, "dR2_external_genemeans_replication_qc_reference.csv"),        row.names = FALSE)
write.csv(s4s, file.path(OUT_DIR, "dR2_external_genemeans_replication_qc_reference_summary.csv"), row.names = FALSE)

## personalized vs reference, per sample
pers <- rbind(transform(disc, study = paste("LBP", tissue))[, c("study","sample","model","gene_mean","dR2")],
              transform(repl, study = cohort)[, c("study","sample","model","gene_mean","dR2")])
cmp <- merge(pers, res[, c("study","sample","model","gene_mean","dR2")],
             by = c("study","sample","model","gene_mean"), suffixes = c("_personal", "_reference")) %>%
  group_by(study, model, gene_mean) %>%
  summarise(samples = n(), pct_personal = 100 * median(dR2_personal), pct_reference = 100 * median(dR2_reference),
            median_abs_diff_pp = 100 * median(abs(dR2_reference - dR2_personal)), .groups = "drop")
write.csv(cmp, file.path(OUT_DIR, "dR2_personal_vs_reference_summary.csv"), row.names = FALSE)
print(as.data.frame(cmp), row.names = FALSE, digits = 3)

## figures: same layout and colours as Fig. 4 / Supp. Fig. 4
plot_dR2 <- function(d, s, rowvar, row_lev, gm_lev, gm_cols, stem) {
  d$model <- factor(d$model, levels = MODELS); s$model <- factor(s$model, levels = MODELS)
  d$row <- factor(d[[rowvar]], levels = row_lev); s$row <- factor(s[[rowvar]], levels = row_lev)
  d$gene_mean <- factor(d$gene_mean, levels = gm_lev); s$gene_mean <- factor(s$gene_mean, levels = gm_lev)
  s <- s %>% group_by(row, model) %>% arrange(gene_mean, .by_group = TRUE) %>%
    mutate(vjust = 1.7 + 1.6 * (row_number() - 1)) %>% ungroup()
  p <- ggplot(d, aes(dR2, fill = gene_mean)) +
    geom_histogram(bins = 45, position = "identity", alpha = 0.55, colour = NA) +
    geom_vline(data = s, aes(xintercept = dR2, colour = gene_mean), linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    geom_text(data = s, aes(x = Inf, y = Inf, label = sprintf("%s  %.2f%%", as.character(gene_mean), pct_added),
                            colour = gene_mean), hjust = 1.08, vjust = s$vjust, size = 3, show.legend = FALSE) +
    scale_fill_manual(values = gm_cols, drop = TRUE, name = NULL) +
    scale_colour_manual(values = gm_cols, drop = TRUE) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    facet_grid(row ~ model, scales = "free_y") +
    labs(x = expression("Additional variance explained (" * Delta * italic(R)^2 * ")"), y = "samples") +
    theme_bw(base_size = 11) + theme(legend.position = "top", strip.text = element_text(face = "bold"))
  ggsave(file.path(PLOT_DIR, paste0(stem, ".pdf")), p, width = 13, height = 7)
  ggsave(file.path(PLOT_DIR, paste0(stem, ".png")), p, width = 13, height = 7, dpi = 300)
}
plot_dR2(f4, f4s, "tissue", c("Brain", "Blood"), c("GTEx", "ROSMAP", "MSBB", "COVID-19"),
         c(GTEx = "#7048E8", ROSMAP = "#0CA678", MSBB = "#F59F00", `COVID-19` = "#E03131"),
         "delta_rsquared_external_genemeans_shared_reference")                      # Supp. Fig. 5
plot_dR2(s4, s4s, "cohort", c("ROSMAP Brain", "COVID Blood"), c("GTEx", "MSBB", "LBP (PM)", "GTEx Blood", "LBP Blood"),
         c(GTEx = "#7048E8", MSBB = "#F59F00", `LBP (PM)` = "#0CA678", `GTEx Blood` = "#7048E8", `LBP Blood` = "#E03131"),
         "delta_rsquared_external_genemeans_replication_qc_reference")              # Supp. Fig. 6
