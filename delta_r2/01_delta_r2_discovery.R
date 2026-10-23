## ============================================================================
## How much does each model add to a gene mean from another cohort?
## Shared gene panel (13,061 brain / 10,562 blood).
##   null <- lm(observed ~ gene_mean);  full <- lm(observed ~ dl + gene_mean)
## Draws Fig. 4 and writes the Supp. Table 1 summary.
## All gene means come from other cohorts (an LBP LOO gene mean would be circular):
##   Brain: GTEx, ROSMAP, MSBB      Blood: GTEx, COVID-19
## ============================================================================
library(arrow); library(ggplot2); library(dplyr)
arrow::set_cpu_count(1); arrow::set_io_thread_count(1)

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUT_DIR  <- file.path(BASE, "evaluation", "complementarity")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)

z <- function(x) (x - mean(x)) / sd(x)
P <- function(...) file.path(BASE, "predictions", ...)

CFG <- list(
  brain = list(bed = "tss_hg38_mane_expressed_brain.bed",
               mats = "Null/mats.rds", obs = "lbp_liv",
               smap = "sample_mapping.tsv", hdr = TRUE,
               gms = c(GTEx = "gtex", ROSMAP = "rosmap", MSBB = "msbb"),
               models = list(
                 Enformer = list(dirs = P("enformer_expanded_personal_brain"), log = TRUE),
                 Basenji2 = list(dirs = P("basenji2_expanded_personal_brain"), log = TRUE),
                 ExPecto  = list(dirs = c(P("expecto_expanded_personal_brain_gpu"),
                                          P("expecto_expanded_personal_brain")), log = FALSE),
                 Xpresso  = list(dirs = P("xpresso_expanded_personal_brain"),  log = FALSE))),
  blood = list(bed = "tss_hg38_mane_expressed_blood.bed",
               mats = "Null/mats_blood.rds", obs = "lbp_blood",
               smap = "samples_blood.tsv", hdr = FALSE,
               gms = c(GTEx = "gtex_blood", `COVID-19` = "covid_blood"),
               models = list(
                 Enformer = list(dirs = P("enformer_expanded_personal_blood"), log = TRUE),
                 Basenji2 = list(dirs = P("basenji2_expanded_personal_blood"), log = TRUE),
                 ExPecto  = list(dirs = c(P("expecto_expanded_personal_blood"),
                                          P("expecto_expanded_personal_blood_gpu")), log = FALSE),
                 Xpresso  = list(dirs = P("xpresso_expanded_personal_blood"),  log = FALSE))))

pick <- function(dirs, s) {
  for (d in dirs) { f <- file.path(d, paste0(s, ".parquet")); if (file.exists(f)) return(f) }
  NA_character_
}

load_dl <- function(dirs, panel, samp, log_pred, obs_profile) {
  files <- vapply(samp, function(s) pick(dirs, s), character(1))
  keep  <- !is.na(files); files <- files[keep]; s_ok <- samp[keep]
  tracks <- setdiff(names(as.data.frame(read_parquet(files[1]))), "gene_id")
  M <- lapply(tracks, function(t)
        matrix(NA_real_, length(panel), length(files), dimnames = list(panel, s_ok)))
  names(M) <- tracks
  for (j in seq_along(files)) {
    d <- as.data.frame(read_parquet(files[j])); i <- match(panel, d$gene_id)
    for (t in tracks) M[[t]][, j] <- d[[t]][i]
  }
  bt <- if (length(tracks) > 1)
          names(which.max(sapply(M, function(m)
            cor(rowMeans(m, na.rm = TRUE), obs_profile, method = "spearman",
                use = "complete.obs")))) else tracks[1]
  mat <- M[[bt]]
  if (log_pred) mat <- log2(mat + 1)
  list(mat = mat, samples = s_ok, track = bt)
}

run_all <- function(obs, dl, gm) {
  t(sapply(colnames(dl), function(s) {
    ok <- is.finite(obs[, s]) & is.finite(dl[, s]) & is.finite(gm)
    Y <- z(obs[ok, s]); D <- z(dl[ok, s]); Mm <- z(gm[ok])
    n <- lm(Y ~ Mm); f <- lm(Y ~ D + Mm)
    c(r2_null = summary(n)$r.squared,
      r2_full = summary(f)$r.squared,
      dR2     = summary(f)$r.squared - summary(n)$r.squared,
      beta_gm_null = unname(coef(n)["Mm"]),
      beta_dl = unname(coef(f)["D"]),
      beta_gm = unname(coef(f)["Mm"]),
      n_genes = sum(ok))
  }))
}

out <- list()
for (tis in names(CFG)) {
  cfg  <- CFG[[tis]]
  bed  <- read.table(file.path(BASE, cfg$bed), sep = "\t",
                     col.names = c("chrom","tss","strand","gene_id","gene_name"),
                     stringsAsFactors = FALSE)
  matsL <- readRDS(file.path(BASE, cfg$mats))
  obs_full <- matsL[[cfg$obs]]

  ## gene panel: bed genes present in observed AND every gene-mean source, with a FANTOM5 symbol match
  panel <- intersect(bed$gene_id, rownames(obs_full))
  for (k in cfg$gms) panel <- intersect(panel, rownames(matsL[[k]]))
  panel <- panel[bed$gene_name[match(panel, bed$gene_id)] %in%
                 read.csv(FANTOM5_MEANS, stringsAsFactors = FALSE)$gene]   # FANTOM5 symbol match, as in cross_gene_correlation

  ## samples with a valid (non-all-NA) parquet in all four models
  sm   <- read.delim(file.path(BASE, cfg$smap), header = cfg$hdr, stringsAsFactors = FALSE)
  cand <- intersect(unique(sm[[2]]), colnames(obs_full))
  ok <- vapply(cand, function(s) all(vapply(names(cfg$models), function(m) {
    f <- pick(cfg$models[[m]]$dirs, s); if (is.na(f)) return(FALSE)
    v <- tryCatch(as.data.frame(read_parquet(f))[[1]], error = function(e) NA)
    !all(is.na(v))
  }, logical(1))), logical(1))
  samp <- cand[ok]

  obs <- obs_full[panel, samp, drop = FALSE]
  cat(sprintf("\n=== %s: %d genes x %d samples | gene means: %s ===\n",
              tis, length(panel), length(samp), paste(names(cfg$gms), collapse = ", ")))

  GM <- lapply(cfg$gms, function(k)
          rowMeans(matsL[[k]][panel, , drop = FALSE], na.rm = TRUE))
  names(GM) <- names(cfg$gms)

  for (mn in names(cfg$models)) {
    spec <- cfg$models[[mn]]
    L <- load_dl(spec$dirs, panel, samp, spec$log, rowMeans(obs))
    cat(sprintf("  %-9s track: %-45s n=%d\n", mn, L$track, length(L$samples)))
    o <- obs[, L$samples, drop = FALSE]
    for (gname in names(GM)) {
      r <- as.data.frame(run_all(o, L$mat, GM[[gname]]))
      r$sample <- rownames(r); rownames(r) <- NULL
      r$model <- mn; r$gene_mean <- gname
      r$tissue <- if (tis == "brain") "Brain" else "Blood"
      out[[paste(tis, mn, gname)]] <- r
    }
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL

summ <- res %>% group_by(tissue, model, gene_mean) %>%
  summarise(r2_null = median(r2_null), r2_full = median(r2_full),
            dR2 = median(dR2), pct_added = 100 * median(dR2),
            beta_gm_null = median(beta_gm_null),
            beta_dl = median(beta_dl), beta_gm = median(beta_gm),
            genes = median(n_genes), n = n(), .groups = "drop")
cat("\n=== medians across samples (shared panel) ===\n")
print(as.data.frame(summ), row.names = FALSE, digits = 3)
write.csv(res,  file.path(OUT_DIR, "dR2_external_genemeans_shared.csv"),         row.names = FALSE)
write.csv(summ, file.path(OUT_DIR, "dR2_external_genemeans_shared_summary.csv"), row.names = FALSE)

res$model  <- factor(res$model,  levels = c("Enformer","Basenji2","ExPecto","Xpresso"))
res$tissue <- factor(res$tissue, levels = c("Brain","Blood"))
GM_LEV  <- c("GTEx","ROSMAP","MSBB","COVID-19")
gm_cols <- c(GTEx = "#7048E8", ROSMAP = "#0CA678", MSBB = "#F59F00", `COVID-19` = "#E03131")
res$gene_mean <- factor(res$gene_mean, levels = GM_LEV)

med <- summ
med$model     <- factor(med$model,  levels = levels(res$model))
med$tissue    <- factor(med$tissue, levels = levels(res$tissue))
med$gene_mean <- factor(med$gene_mean, levels = GM_LEV)
med <- med %>% group_by(tissue, model) %>% arrange(gene_mean, .by_group = TRUE) %>%
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
  facet_grid(tissue ~ model, scales = "free_y") +
  labs(title = "Extra cross-gene variance explained by adding the model to an external gene mean",
       subtitle = expression("shared gene panel; one value per sample; dashed lines and corner labels are the median " * Delta * italic(R)^2),
       x = expression("Additional variance explained (" * Delta * italic(R)^2 * ")"),
       y = "samples") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOT_DIR, "delta_rsquared_external_genemeans_shared.pdf"), p, width = 13, height = 7)
cat("\nsaved delta_rsquared_external_genemeans_shared.pdf\n")
#        Brain   Xpresso      ROSMAP   0.509  0.517  0.829   0.714     0.688     0.094
