#!/usr/bin/env Rscript
## Cross-gene correlations for every arm and every sample. Writes
## arm_scores_allcohorts_qcexact.csv, which the figures and tests use.
## Arms: every track of the four models, and the gene-mean baselines in Fig. 1
## (brain: LBP PM, MSBB, GTEx, FANTOM5 CAGE; blood: GTEx, FANTOM5 CAGE).
## Within a cohort, every arm uses the same genes and samples:
##   genes   = expanded bed ∩ observed ∩ all gene-mean sources ∩ FANTOM5 symbols
##   samples = mapping ∩ observed ∩ has a non-NA parquet in all four models,
##             minus samples that fail QC
## QC: a sample is dropped if its expression correlates poorly with its own
## cohort's leave-one-out gene mean (Spearman < 0.2), which points to a problem
## with the RNA-seq data. This drops two COVID-19 samples. They are listed in
## flagged_low_genemean_samples.csv for 02_score_personal_vs_reference.R and
## delta_r2/02_delta_r2_replication.R.
suppressPackageStartupMessages(library(arrow))
arrow::set_cpu_count(1); arrow::set_io_thread_count(1)
source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUTDIR <- file.path(BASE, "evaluation", "allcohorts")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
FANTOM <- read.csv(FANTOM5_MEANS, stringsAsFactors = FALSE)
MODELS <- c("Enformer","Basenji2","ExPecto","Xpresso")
QC_MIN_RHO <- 0.2
gid_of <- function(d) if ("gene_id" %in% names(d)) d$gene_id else rownames(d)
P <- function(...) file.path(BASE, "predictions", ...)
R <- function(...) file.path(BASE, "Replication", ...)

BRAIN_GM <- list(`Gene-mean: LBP PM`="lbp_pm", `Gene-mean: MSBB`="msbb", `Gene-mean: GTEx`="gtex")
BLOOD_GM <- list(`Gene-mean: GTEx Blood`="gtex_blood")
CFG <- list(
  `LBP-Brain` = list(
    bed="tss_hg38_mane_expressed_brain.bed", mats="Null/mats.rds", obs="lbp_liv",
    smap="sample_mapping.tsv", hdr=TRUE, fcol="fantom_brain_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=BRAIN_GM,
    dirs=list(Enformer=P("enformer_expanded_personal_brain"), Basenji2=P("basenji2_expanded_personal_brain"),
              ExPecto =P("expecto_expanded_personal_brain"),  Xpresso =P("xpresso_expanded_personal_brain"))),
  `LBP-Blood` = list(
    bed="tss_hg38_mane_expressed_blood.bed", mats="Null/mats_blood.rds", obs="lbp_blood",
    smap="samples_blood.tsv", hdr=FALSE, fcol="fantom_blood_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=BLOOD_GM,
    dirs=list(Enformer=P("enformer_expanded_personal_blood"), Basenji2=P("basenji2_expanded_personal_blood"),
              ExPecto =P("expecto_expanded_personal_blood"),  Xpresso =P("xpresso_expanded_personal_blood"))),
  `ROSMAP-Brain` = list(
    bed="tss_hg38_mane_expressed_brain.bed", mats="Null/mats.rds", obs="rosmap",
    smap="sample_mapping_rosmap.tsv", hdr=TRUE, fcol="fantom_brain_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=BRAIN_GM,
    dirs=list(Enformer=R("ROSMAP/predictions/enformer_expanded_personal_brain"),
              Basenji2=R("ROSMAP/predictions/basenji2_expanded_personal_brain"),
              ExPecto =R("ROSMAP/predictions/expecto_expanded_personal_brain"),
              Xpresso =R("ROSMAP/predictions/xpresso_expanded_personal_brain"))),
  `COVID-Blood` = list(
    bed="tss_hg38_mane_expressed_blood.bed", mats="Null/mats_blood.rds", obs="covid_blood",
    smap="samples_covid.tsv", hdr=FALSE, fcol="fantom_blood_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=BLOOD_GM,
    dirs=list(Enformer=R("COVID/predictions/enformer_expanded_personal_blood"),
              Basenji2=R("COVID/predictions/basenji2_expanded_personal_blood"),
              ExPecto =R("COVID/predictions/expecto_expanded_personal_blood"),
              Xpresso =R("COVID/predictions/xpresso_expanded_personal_blood"))))

pick <- function(dir, s) { f <- file.path(dir, paste0(s, ".parquet")); if (file.exists(f)) f else NA_character_ }

all_rows <- list(); flagged <- list()
for (co in names(CFG)) {
  cfg  <- CFG[[co]]
  bed  <- read.table(file.path(BASE,cfg$bed), sep="\t",
                     col.names=c("chrom","tss","strand","gene_id","gene_name"), stringsAsFactors=FALSE)
  mats <- readRDS(file.path(BASE,cfg$mats)); obs_full <- mats[[cfg$obs]]

  panel <- intersect(bed$gene_id, rownames(obs_full))
  for (s in cfg$gm) panel <- intersect(panel, rownames(mats[[s]]))
  panel <- panel[bed$gene_name[match(panel, bed$gene_id)] %in% FANTOM$gene]
  sym   <- bed$gene_name[match(panel, bed$gene_id)]

  sm   <- read.delim(file.path(BASE,cfg$smap), header=cfg$hdr, stringsAsFactors=FALSE)
  cand <- intersect(unique(sm[[2]]), colnames(obs_full))

  ok <- vapply(cand, function(s) all(vapply(MODELS, function(m) {
    f <- pick(cfg$dirs[[m]], s); if (is.na(f)) return(FALSE)
    v <- tryCatch(as.data.frame(read_parquet(f))[[1]], error=function(e) NA)
    !all(is.na(v))
  }, logical(1))), logical(1))
  samp <- cand[ok]

  ## QC: drop samples that barely match the cohort's own leave-one-out gene mean (only 2 samples were dropped in ROSMAP)
  o0  <- obs_full[panel, samp, drop=FALSE]; rs0 <- rowSums(o0); n0 <- ncol(o0)
  r0  <- vapply(seq_len(n0), function(i)
           cor((rs0 - o0[,i])/(n0-1), o0[,i], method="spearman", use="complete.obs"), numeric(1))
  bad <- which(r0 < QC_MIN_RHO)
  if (length(bad)) {
    flagged[[co]] <- data.frame(cohort=co, sample=samp[bad], spearman=r0[bad], stringsAsFactors=FALSE)
    cat(sprintf("  [QC] dropping %d sample(s) from %s (leave-one-out gene-mean rho < %.1f)\n", length(bad), co, QC_MIN_RHO))
    samp <- samp[-bad]
  }
  cat(sprintf("\n### %s: %d genes x %d samples (of %d in mapping, %d with observed expr)\n",
              co, length(panel), length(samp), length(unique(sm[[2]])), length(cand)))
  obs <- obs_full[panel, samp, drop=FALSE]; N <- ncol(obs)

  ## gene-mean baselines: the same predictor for every sample
  add <- function(gm, nm) {
    sp <- vapply(seq_len(N), function(i) cor(gm, obs[,i], method="spearman", use="complete.obs"), numeric(1))
    pe <- vapply(seq_len(N), function(i) cor(gm, obs[,i], use="complete.obs"), numeric(1))
    all_rows[[paste(co,nm)]] <<- data.frame(cohort=co, arm=nm, model="Gene-mean", track=nm,
      sample=samp, spearman=sp, pearson=pe, stringsAsFactors=FALSE)
  }
  for (nm in names(cfg$gm)) add(rowMeans(mats[[cfg$gm[[nm]]]][panel, , drop=FALSE], na.rm=TRUE), nm)
  add(log2(FANTOM[[cfg$fcol]][match(sym, FANTOM$gene)] + 1), "Gene-mean: FANTOM5 CAGE")

  ## models: one arm per track
  for (mi in seq_along(MODELS)) {
    mn <- MODELS[mi]; lg <- cfg$log[mi]; tracks <- NULL; sp <- pe <- NULL
    for (i in seq_along(samp)) {
      f <- pick(cfg$dirs[[mn]], samp[i]); if (is.na(f)) next
      p <- as.data.frame(read_parquet(f)); j <- match(panel, gid_of(p))
      use <- setdiff(names(p), "gene_id")
      if (is.null(tracks)) { tracks <- use
        sp <- pe <- matrix(NA_real_, length(samp), length(use), dimnames=list(samp,use)) }
      for (k in seq_along(use)) {
        v <- p[[use[k]]][j]; vp <- if (lg) log2(v+1) else v
        sp[i,k] <- cor(v,  obs[,i], method="spearman", use="complete.obs")
        pe[i,k] <- cor(vp, obs[,i], use="complete.obs")
      }
    }
    for (k in seq_along(tracks))
      all_rows[[paste(co,mn,k)]] <- data.frame(cohort=co,
        arm=sprintf("%s: %s", mn, tracks[k]), model=mn, track=tracks[k],
        sample=samp, spearman=sp[,k], pearson=pe[,k], stringsAsFactors=FALSE)
    cat(sprintf("  %-9s %d track(s), %d/%d scored\n", mn, length(tracks), sum(!is.na(sp[,1])), N))
  }
}

fl <- if (length(flagged)) do.call(rbind, flagged) else data.frame(cohort=character(0), sample=character(0), spearman=numeric(0))
write.csv(fl, file.path(OUTDIR,"flagged_low_genemean_samples.csv"), row.names=FALSE)
cat("\nQC-flagged samples:\n"); print(fl, row.names=FALSE)

res <- do.call(rbind, all_rows); rownames(res) <- NULL
write.csv(res, file.path(OUTDIR,"arm_scores_allcohorts_qcexact.csv"), row.names=FALSE)
cat("\nwrote", file.path(OUTDIR,"arm_scores_allcohorts_qcexact.csv"), "\n")
