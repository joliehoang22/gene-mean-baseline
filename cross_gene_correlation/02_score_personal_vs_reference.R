#!/usr/bin/env Rscript
## Scores every model track twice, with personalized and with reference-genome
## predictions, on the same genes and samples as 01 (QC-flagged samples dropped).
## Writes arm_scores_allcohorts_withref.csv with a `mode` column
## (personalized / reference), used by figures/Fig3.R.
suppressPackageStartupMessages(library(arrow))
arrow::set_cpu_count(1); arrow::set_io_thread_count(1)
source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUTDIR  <- file.path(BASE, "evaluation", "allcohorts")
FANTOM  <- read.csv(FANTOM5_MEANS, stringsAsFactors = FALSE)
TARGETS <- read.delim(file.path(BASE, "targets_human.txt"), stringsAsFactors = FALSE)
FLAGGED <- read.csv(file.path(OUTDIR, "flagged_low_genemean_samples.csv"), stringsAsFactors = FALSE)
MODELS  <- c("Enformer","Basenji2","ExPecto","Xpresso")
gid_of  <- function(d) if ("gene_id" %in% names(d)) d$gene_id else rownames(d)
P <- function(...) file.path(BASE, "predictions", ...)
R <- function(...) file.path(BASE, "Replication", ...)

## reference-genome predictions are shared by both cohorts of a tissue.
## Basenji2 names its reference columns track_<index>; the others use the track names.
REF_FILE <- function(model, tis) P(switch(model,
  Enformer = sprintf("enformer_reference_expanded_%s/reference_genome.parquet", tis),
  Basenji2 = sprintf("basenji2_reference_expanded_%s/reference_genome.parquet", tis),
  ExPecto  = sprintf("expecto_reference_expanded_%s/reference_genome.parquet",  tis),
  Xpresso  = sprintf("xpresso_reference_expanded_%s/seqonly.parquet",           tis)))

## gm = gene-mean sources, used only to define the same gene panel as 01
BRAIN_GM <- c("lbp_pm", "msbb", "gtex"); BLOOD_GM <- "gtex_blood"
CFG <- list(
  `LBP-Brain` = list(
    bed="tss_hg38_mane_expressed_brain.bed", mats="Null/mats.rds", obs="lbp_liv", tis="brain",
    smap="sample_mapping.tsv", hdr=TRUE, log=c(TRUE,TRUE,FALSE,FALSE), gm=BRAIN_GM,
    dirs=list(Enformer=P("enformer_expanded_personal_brain"), Basenji2=P("basenji2_expanded_personal_brain"),
              ExPecto =P("expecto_expanded_personal_brain"),  Xpresso =P("xpresso_expanded_personal_brain"))),
  `LBP-Blood` = list(
    bed="tss_hg38_mane_expressed_blood.bed", mats="Null/mats_blood.rds", obs="lbp_blood", tis="blood",
    smap="samples_blood.tsv", hdr=FALSE, log=c(TRUE,TRUE,FALSE,FALSE), gm=BLOOD_GM,
    dirs=list(Enformer=P("enformer_expanded_personal_blood"), Basenji2=P("basenji2_expanded_personal_blood"),
              ExPecto =P("expecto_expanded_personal_blood"),  Xpresso =P("xpresso_expanded_personal_blood"))),
  `ROSMAP-Brain` = list(
    bed="tss_hg38_mane_expressed_brain.bed", mats="Null/mats.rds", obs="rosmap", tis="brain",
    smap="sample_mapping_rosmap.tsv", hdr=TRUE, log=c(TRUE,TRUE,FALSE,FALSE), gm=BRAIN_GM,
    dirs=list(Enformer=R("ROSMAP/predictions/enformer_expanded_personal_brain"),
              Basenji2=R("ROSMAP/predictions/basenji2_expanded_personal_brain"),
              ExPecto =R("ROSMAP/predictions/expecto_expanded_personal_brain"),
              Xpresso =R("ROSMAP/predictions/xpresso_expanded_personal_brain"))),
  `COVID-Blood` = list(
    bed="tss_hg38_mane_expressed_blood.bed", mats="Null/mats_blood.rds", obs="covid_blood", tis="blood",
    smap="samples_covid.tsv", hdr=FALSE, log=c(TRUE,TRUE,FALSE,FALSE), gm=BLOOD_GM,
    dirs=list(Enformer=R("COVID/predictions/enformer_expanded_personal_blood"),
              Basenji2=R("COVID/predictions/basenji2_expanded_personal_blood"),
              ExPecto =R("COVID/predictions/expecto_expanded_personal_blood"),
              Xpresso =R("COVID/predictions/xpresso_expanded_personal_blood"))))

pick <- function(dir, s) { f <- file.path(dir, paste0(s, ".parquet")); if (file.exists(f)) f else NA_character_ }

all_rows <- list()
for (co in names(CFG)) {
  cfg  <- CFG[[co]]
  bed  <- read.table(file.path(BASE,cfg$bed), sep="\t",
                     col.names=c("chrom","tss","strand","gene_id","gene_name"), stringsAsFactors=FALSE)
  mats <- readRDS(file.path(BASE,cfg$mats)); obs_full <- mats[[cfg$obs]]

  panel <- intersect(bed$gene_id, rownames(obs_full))
  for (s in cfg$gm) panel <- intersect(panel, rownames(mats[[s]]))
  panel <- panel[bed$gene_name[match(panel, bed$gene_id)] %in% FANTOM$gene]

  sm   <- read.delim(file.path(BASE,cfg$smap), header=cfg$hdr, stringsAsFactors=FALSE)
  cand <- intersect(unique(sm[[2]]), colnames(obs_full))
  ok <- vapply(cand, function(s) all(vapply(MODELS, function(m) {
    f <- pick(cfg$dirs[[m]], s); if (is.na(f)) return(FALSE)
    v <- tryCatch(as.data.frame(read_parquet(f))[[1]], error=function(e) NA)
    !all(is.na(v))
  }, logical(1))), logical(1))
  samp <- setdiff(cand[ok], FLAGGED$sample[FLAGGED$cohort == co])
  cat(sprintf("\n### %s: %d genes x %d samples\n", co, length(panel), length(samp)))
  obs <- obs_full[panel, samp, drop=FALSE]; N <- ncol(obs)

  row <- function(mn, track, mode, sp, pe) data.frame(cohort=co, arm=sprintf("%s: %s", mn, track), model=mn,
    track=track, mode=mode, sample=samp, spearman=sp, pearson=pe, stringsAsFactors=FALSE)

  for (mi in seq_along(MODELS)) {
    mn <- MODELS[mi]; lg <- cfg$log[mi]; sp <- pe <- NULL

    ## personalized: one parquet per sample, scored as it is read
    for (i in seq_len(N)) {
      p <- as.data.frame(read_parquet(pick(cfg$dirs[[mn]], samp[i]))); j <- match(panel, gid_of(p))
      if (is.null(sp)) { tracks <- setdiff(names(p), "gene_id")
        sp <- pe <- matrix(NA_real_, N, length(tracks), dimnames=list(samp, tracks)) }
      for (tr in tracks) {
        v <- p[[tr]][j]; vp <- if (lg) log2(v+1) else v
        sp[i,tr] <- cor(v,  obs[,i], method="spearman", use="complete.obs")
        pe[i,tr] <- cor(vp, obs[,i], use="complete.obs")
      }
    }
    for (tr in tracks) all_rows[[paste(co, mn, tr, "pers")]] <- row(mn, tr, "personalized", sp[,tr], pe[,tr])

    ## reference genome: one parquet per tissue, the same prediction for every sample
    cols <- if (mn == "Basenji2") paste0("track_", TARGETS$index[match(tracks, TARGETS$description)]) else tracks
    rp <- as.data.frame(read_parquet(REF_FILE(mn, cfg$tis)))
    stopifnot(all(cols %in% names(rp)))
    for (k in seq_along(tracks)) {
      v <- rp[[cols[k]]][match(panel, gid_of(rp))]; vp <- if (lg) log2(v+1) else v
      all_rows[[paste(co, mn, tracks[k], "ref")]] <- row(mn, tracks[k], "reference",
        vapply(seq_len(N), function(i) cor(v,  obs[,i], method="spearman", use="complete.obs"), numeric(1)),
        vapply(seq_len(N), function(i) cor(vp, obs[,i], use="complete.obs"), numeric(1)))
    }
    cat(sprintf("  %-9s %d track(s), personalized and reference\n", mn, length(tracks)))
  }
}

res <- do.call(rbind, all_rows); rownames(res) <- NULL
write.csv(res, file.path(OUTDIR,"arm_scores_allcohorts_withref.csv"), row.names=FALSE)
cat("\nwrote", file.path(OUTDIR,"arm_scores_allcohorts_withref.csv"), "\n")
