## ============================================================================
## How much does each model add to a gene mean from another cohort?
## Shared gene panel (13,061 brain / 10,562 blood).
##   null <- lm(observed ~ gene_mean);  full <- lm(observed ~ dl + gene_mean)
## Gene means never come from the study cohort itself:
##   Brain (LBP LIV): GTEx, ROSMAP, MSBB, LBP PM, FANTOM5 CAGE
##   Blood (LBP LIV): GTEx, COVID-19, FANTOM5 CAGE
## FANTOM5 gene mean = log2(mean CAGE TPM + 1), matched by gene symbol, as in the scoring.
## Writes the tables behind Fig. 2 (brain), Supp. Fig. 5 (blood rows) and Supp. Tables 2-3.
## ============================================================================
library(arrow); library(dplyr)
arrow::set_cpu_count(1); arrow::set_io_thread_count(1)

source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUT_DIR <- file.path(BASE, "evaluation", "complementarity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FANTOM  <- read.csv(FANTOM5_MEANS, stringsAsFactors = FALSE)

z <- function(x) (x - mean(x)) / sd(x)
P <- function(...) file.path(BASE, "predictions", ...)

CFG <- list(
  brain = list(bed = "tss_hg38_mane_expressed_brain.bed",
               mats = "Null/mats.rds", obs = "lbp_liv",
               smap = "sample_mapping.tsv", hdr = TRUE, fcol = "fantom_brain_mean",
               gms = c(GTEx = "gtex", ROSMAP = "rosmap", MSBB = "msbb", `LBP PM` = "lbp_pm"),
               models = list(
                 Enformer = list(dirs = P("enformer_expanded_personal_brain"), log = TRUE),
                 Basenji2 = list(dirs = P("basenji2_expanded_personal_brain"), log = TRUE),
                 ExPecto  = list(dirs = c(P("expecto_expanded_personal_brain_gpu"),
                                          P("expecto_expanded_personal_brain")), log = FALSE),
                 Xpresso  = list(dirs = P("xpresso_expanded_personal_brain"),  log = FALSE))),
  blood = list(bed = "tss_hg38_mane_expressed_blood.bed",
               mats = "Null/mats_blood.rds", obs = "lbp_blood",
               smap = "samples_blood.tsv", hdr = FALSE, fcol = "fantom_blood_mean",
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
  sym   <- bed$gene_name[match(panel, bed$gene_id)]
  panel <- panel[sym %in% FANTOM$gene]; sym <- bed$gene_name[match(panel, bed$gene_id)]

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
  GM <- c(lapply(cfg$gms, function(k) rowMeans(matsL[[k]][panel, , drop = FALSE], na.rm = TRUE)),
          list(`FANTOM5 CAGE` = log2(FANTOM[[cfg$fcol]][match(sym, FANTOM$gene)] + 1)))
  cat(sprintf("\n=== %s: %d genes x %d samples | gene means: %s ===\n",
              tis, length(panel), length(samp), paste(names(GM), collapse = ", ")))

  for (mn in names(cfg$models)) {
    spec <- cfg$models[[mn]]
    L <- load_dl(spec$dirs, panel, samp, spec$log, rowMeans(obs))
    cat(sprintf("  %-9s track: %-45s n=%d\n", mn, L$track, length(L$samples)))
    o <- obs[, L$samples, drop = FALSE]
    for (gname in names(GM)) {
      r <- as.data.frame(run_all(o, L$mat, GM[[gname]]))
      r$sample <- rownames(r); rownames(r) <- NULL
      r$model <- mn; r$gene_mean <- gname; r$track <- L$track
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
cat("\nFig. 2 and Supp. Fig. 5 are drawn by figures/Fig2.R and figures/SuppFig5.R\n")
