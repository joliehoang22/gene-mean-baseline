#!/usr/bin/env Rscript
## Same scoring as 01, plus the reference-genome arm for every model, with a
## `mode` column (personalized / reference). Writes arm_scores_allcohorts_withref.csv.
## Within a cohort, every arm uses the same genes and samples:
##   genes   = expanded bed ∩ observed ∩ all gene-mean sources ∩ FANTOM5 symbols
##   samples = mapping ∩ observed ∩ has a non-NA parquet in all four models
suppressPackageStartupMessages(library(arrow))
arrow::set_cpu_count(1); arrow::set_io_thread_count(1)
source("config/paths.R")   # BASE, PLOT_DIR, FANTOM5_MEANS
OUTDIR <- file.path(BASE, "evaluation", "allcohorts")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
FANTOM <- read.csv(FANTOM5_MEANS,
                   stringsAsFactors = FALSE)
MODELS <- c("Enformer","Basenji2","ExPecto","Xpresso")
gid_of <- function(d) if ("gene_id" %in% names(d)) d$gene_id else rownames(d)
P <- function(...) file.path(BASE, "predictions", ...)
R <- function(...) file.path(BASE, "Replication", ...)

CFG <- list(
  `LBP-Brain` = list(
    bed="tss_hg38_mane_expressed_brain.bed", mats="Null/mats.rds", obs="lbp_liv",
    smap="sample_mapping.tsv", hdr=TRUE, fcol="fantom_brain_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=list(`Gene-mean: LBP LIV (LOO)`="LOO", `Gene-mean: LBP PM`="lbp_pm",
            `Gene-mean: MSBB`="msbb", `Gene-mean: GTEx`="gtex", `Gene-mean: ROSMAP`="rosmap"),
    dirs=list(Enformer=P("enformer_expanded_personal_brain"),
              Basenji2=P("basenji2_expanded_personal_brain"),
              ExPecto =c(P("expecto_expanded_personal_brain_gpu"), P("expecto_expanded_personal_brain")),
              Xpresso =P("xpresso_expanded_personal_brain"))),
  `LBP-Blood` = list(
    bed="tss_hg38_mane_expressed_blood.bed", mats="Null/mats_blood.rds", obs="lbp_blood",
    smap="samples_blood.tsv", hdr=FALSE, fcol="fantom_blood_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=list(`Gene-mean: LBP Blood (LOO)`="LOO", `Gene-mean: GTEx Blood`="gtex_blood",
            `Gene-mean: COVID Blood`="covid_blood"),
    dirs=list(Enformer=P("enformer_expanded_personal_blood"),
              Basenji2=P("basenji2_expanded_personal_blood"),
              ExPecto =c(P("expecto_expanded_personal_blood"), P("expecto_expanded_personal_blood_gpu")),
              Xpresso =P("xpresso_expanded_personal_blood"))),
  `ROSMAP-Brain` = list(
    bed="tss_hg38_mane_expressed_brain.bed", mats="Null/mats.rds", obs="rosmap",
    smap="sample_mapping_rosmap.tsv", hdr=TRUE, fcol="fantom_brain_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=list(`Gene-mean: ROSMAP (LOO)`="LOO", `Gene-mean: LBP LIV`="lbp_liv",
            `Gene-mean: LBP PM`="lbp_pm", `Gene-mean: MSBB`="msbb", `Gene-mean: GTEx`="gtex"),
    dirs=list(Enformer=R("ROSMAP/predictions/enformer_expanded_personal_brain"),
              Basenji2=R("ROSMAP/predictions/basenji2_expanded_personal_brain"),
              ExPecto =c(R("ROSMAP/predictions/expecto_expanded_personal_brain"),
                         R("ROSMAP/predictions/expecto_expanded_personal_brain_gpu")),
              Xpresso =R("ROSMAP/predictions/xpresso_expanded_personal_brain"))),
  `COVID-Blood` = list(
    bed="tss_hg38_mane_expressed_blood.bed", mats="Null/mats_blood.rds", obs="covid_blood",
    smap="samples_covid.tsv", hdr=FALSE, fcol="fantom_blood_mean", log=c(TRUE,TRUE,FALSE,FALSE),
    gm=list(`Gene-mean: COVID Blood (LOO)`="LOO", `Gene-mean: LBP Blood`="lbp_blood",
            `Gene-mean: GTEx Blood`="gtex_blood"),
    dirs=list(Enformer=R("COVID/predictions/enformer_expanded_personal_blood"),
              Basenji2=R("COVID/predictions/basenji2_expanded_personal_blood"),
              ExPecto =c(R("COVID/predictions/expecto_expanded_personal_blood"),
                         R("COVID/predictions/expecto_expanded_personal_blood_gpu")),
              Xpresso =R("COVID/predictions/xpresso_expanded_personal_blood"))))

pick <- function(dirs, s) { for (d in dirs) { f <- file.path(d, paste0(s,".parquet")); if (file.exists(f)) return(f) }; NA_character_ }

all_rows <- list()
for (co in names(CFG)) {
  cfg  <- CFG[[co]]
  bed  <- read.table(file.path(BASE,cfg$bed), sep="\t",
                     col.names=c("chrom","tss","strand","gene_id","gene_name"), stringsAsFactors=FALSE)
  mats <- readRDS(file.path(BASE,cfg$mats)); obs_full <- mats[[cfg$obs]]
  if (is.null(obs_full)) { cat("\n###", co, ": observed matrix", cfg$obs, "not found, skipping\n"); next }

  panel <- intersect(bed$gene_id, rownames(obs_full))
  for (s in cfg$gm) if (s != "LOO") panel <- intersect(panel, rownames(mats[[s]]))
  sym0  <- bed$gene_name[match(panel, bed$gene_id)]
  panel <- panel[!is.na(match(sym0, FANTOM$gene))]
  sym   <- bed$gene_name[match(panel, bed$gene_id)]

  sm   <- read.delim(file.path(BASE,cfg$smap), header=cfg$hdr, stringsAsFactors=FALSE)
  cand <- intersect(unique(sm[[2]]), colnames(obs_full))

  ok <- vapply(cand, function(s) all(vapply(MODELS, function(m) {
    f <- pick(cfg$dirs[[m]], s); if (is.na(f)) return(FALSE)
    v <- tryCatch(as.data.frame(read_parquet(f))[[1]], error=function(e) NA)
    !all(is.na(v))
  }, logical(1))), logical(1))
  samp <- cand[ok]
  ## --- QC: drop the samples flagged by 01 before gene means are built ---
  .exf <- file.path(BASE, "evaluation/allcohorts/flagged_low_genemean_samples.csv")
  if (!exists(".EXCL")) .EXCL <- if (file.exists(.exf)) read.csv(.exf, stringsAsFactors=FALSE) else data.frame(cohort=character(0), sample=character(0))
  .exs <- .EXCL$sample[.EXCL$cohort == co]
  if (length(.exs) && any(samp %in% .exs)) {
    cat(sprintf("  [QC] dropping %d flagged sample(s) from %s\n", sum(samp %in% .exs), co))
    samp <- samp[!(samp %in% .exs)]
  }
  cat(sprintf("\n### %s: %d genes x %d samples (of %d in mapping, %d with observed expr)\n",
              co, length(panel), length(samp), length(unique(sm[[2]])), length(cand)))
  if (!length(samp)) next
  obs <- obs_full[panel, samp, drop=FALSE]; N <- ncol(obs)

  add <- function(gm, nm) {
    sp <- vapply(seq_len(N), function(i) cor(gm[,i], obs[,i], method="spearman", use="complete.obs"), numeric(1))
    pe <- vapply(seq_len(N), function(i) cor(gm[,i], obs[,i], use="complete.obs"), numeric(1))
    all_rows[[paste(co,nm)]] <<- data.frame(cohort=co, arm=nm, model="Gene-mean", track=nm, mode="genemean",
      sample=samp, spearman=sp, pearson=pe, stringsAsFactors=FALSE)
  }
  rs  <- rowSums(obs)
  loo <- vapply(seq_len(N), function(i) (rs-obs[,i])/(N-1), numeric(length(panel)))
  dimnames(loo) <- list(panel, samp)
  for (nm in names(cfg$gm)) {
    s <- cfg$gm[[nm]]
    add(if (s=="LOO") loo else matrix(rowMeans(mats[[s]][panel,,drop=FALSE], na.rm=TRUE),
        nrow=length(panel), ncol=N, dimnames=list(panel,samp)), nm)
  }
  add(matrix(log2(FANTOM[[cfg$fcol]][match(sym, FANTOM$gene)] + 1),
             nrow=length(panel), ncol=N, dimnames=list(panel,samp)), "Gene-mean: FANTOM5 CAGE")

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
    if (is.null(tracks)) { cat(sprintf("  %-9s SKIP\n", mn)); next }
    for (k in seq_along(tracks))
      all_rows[[paste(co,mn,k)]] <- data.frame(cohort=co,
        arm=sprintf("%s: %s", mn, tracks[k]), model=mn, track=tracks[k], mode="personalized",
        sample=samp, spearman=sp[,k], pearson=pe[,k], stringsAsFactors=FALSE)
    cat(sprintf("  %-9s %d track(s), %d/%d scored\n", mn, length(tracks), sum(!is.na(sp[,1])), N))
    ## ---- reference-genome arm: same panel/samples, tracks matched to the personalized arm ----
    .PRED <- file.path(BASE, "predictions")
    if (!exists("TARGETS"))
      TARGETS <- read.delim(file.path(BASE, "targets_human.txt"),
                            stringsAsFactors = FALSE)
    .tis  <- if (grepl("Brain", co)) "brain" else "blood"
    .spec <- switch(mn,
      Enformer = list(f = sprintf("enformer_reference_expanded_%s/reference_genome.parquet", .tis), kind = "named"),
      Basenji2 = list(f = sprintf("basenji2_reference_expanded_%s/reference_genome.parquet", .tis), kind = "indexed"),
      ExPecto  = list(f = sprintf("expecto_reference_expanded_%s/reference_genome.parquet",  .tis), kind = "named"),
      Xpresso  = list(f = sprintf("xpresso_reference_expanded_%s/seqonly.parquet",           .tis), kind = "named"),
      NULL)
    .rf <- if (is.null(.spec)) NA_character_ else file.path(.PRED, .spec$f)
    if (is.na(.rf) || !file.exists(.rf)) {
      cat(sprintf("  %-9s reference: MISSING (%s)\n", mn, .rf))
    } else {
      .use <- NULL
      if (.spec$kind == "indexed") {
        .idx <- TARGETS$index[match(tracks, TARGETS$description)]
        if (anyNA(.idx)) cat(sprintf("  %-9s reference: FAIL, %d track name(s) not in targets_human.txt\n",
                                     mn, sum(is.na(.idx))))
        else { .use <- paste0("track_", .idx); .nms <- tracks }
      } else { .use <- tracks; .nms <- tracks }
      if (!is.null(.use)) {
        .avail <- open_dataset(.rf)$schema$names
        if (!all(.use %in% .avail)) {
          cat(sprintf("  %-9s reference: FAIL, %d/%d columns absent in %s\n",
                      mn, sum(!(.use %in% .avail)), length(.use), basename(.rf)))
        } else {
          .p <- as.data.frame(read_parquet(.rf, col_select = all_of(c("gene_id", .use))))
          .j <- match(panel, gid_of(.p))
          for (k in seq_along(.use)) {
            .v  <- .p[[.use[k]]][.j]; .vp <- if (lg) log2(.v + 1) else .v
            .spv <- vapply(seq_len(N), function(i) cor(.v,  obs[,i], method="spearman", use="complete.obs"), numeric(1))
            .pev <- vapply(seq_len(N), function(i) cor(.vp, obs[,i], use="complete.obs"), numeric(1))
            all_rows[[paste(co,mn,"ref",k)]] <- data.frame(cohort = co,
              arm = sprintf("%s: %s", mn, .nms[k]), model = mn, track = .nms[k], mode = "reference",
              sample = samp, spearman = .spv, pearson = .pev, stringsAsFactors = FALSE)
          }
          cat(sprintf("  %-9s reference: %d track(s) from %s\n", mn, length(.use), basename(.rf)))
        }
      }
    }
  }
}

res <- do.call(rbind, all_rows); rownames(res) <- NULL
write.csv(res, file.path(OUTDIR,"arm_scores_allcohorts_withref.csv"), row.names=FALSE)
cat("\nwrote", file.path(OUTDIR,"arm_scores_allcohorts_withref.csv"), "\n")
