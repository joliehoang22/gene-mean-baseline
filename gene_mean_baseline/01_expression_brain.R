################################################################################
## 01_expression_brain.R
##
## Observed expression for the five brain datasets. The results are saved as
## .rds files, which later scripts load.
##
## All five datasets are put on the same scale: log10(TPM + 1), with
## version-stripped ENSG IDs.
##   GTEx          published TPM  -> log
##   MSBB, ROSMAP  raw counts + featureCounts Length column -> TPM
##   LBP LIV/PM    raw counts + v30 union-exon lengths -> TPM
##
## We use TPM because it divides out transcript length, so it is proportional
## to molar abundance. CAGE tags the 5' end of each molecule once, so it is also
## molar. CPM would keep a transcript-length factor that CAGE does not have.
## Xpresso and ExPecto were also trained on RPKM- or TPM-scaled values.
##
## Runtime is about 10-15 min, mostly reading the GTF and the GTEx table.
##
## OUTPUTS (all under Null/)
##   mats.rds         list of 5 matrices, restricted to `genes`
##   genes.rds        expression-filtered common gene universe
##   common.rds       genes present in all 5 (pre-filter)
##   gene_means.rds   per-dataset gene-mean vector
##   gtex_log10tpm_full.rds   unfiltered GTEx, for threshold sweeps
################################################################################

library(data.table)
library(rtracklayer)

source("config/paths.R")   # BASE, need()
OUT_DIR <- file.path(BASE, "Null")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GTF_V30        <- need("GENCODE_V30_GTF")
GTEX_TPM_FILE  <- need("GTEX_TPM")
GTEX_META_FILE <- need("GTEX_SAMPLE_ATTR")
GTEX_TISSUE    <- "Brain - Frontal Cortex (BA9)"
MSBB_FILE      <- need("MSBB_COUNTS")
ROSMAP_FILE    <- need("ROSMAP_COUNTS")
LBP_FILE       <- need("LBP_BRAIN_RDS")
TPM_MIN      <- 1      # gene must reach >= this TPM ...
FRAC_SAMPLES <- 0.25   # ... in >= this fraction of GTEx samples

## ── helpers ─────────────────────────────────────────────────────────────────
load_gene_lengths <- function(gtf_path) {
  gtf   <- import(gtf_path, "gtf")
  exons <- gtf[gtf$type == "exon"]
  len   <- sum(width(GenomicRanges::reduce(split(exons, exons$gene_id))))
  names(len) <- sub("\\..*", "", names(len))
  len[!duplicated(names(len))]
}

counts_to_log10tpm <- function(counts, gene_len) {
  counts <- as.matrix(counts)
  rownames(counts) <- sub("\\..*", "", rownames(counts))
  L  <- gene_len[rownames(counts)]
  ok <- !is.na(L)
  counts <- counts[ok, , drop = FALSE]; L <- L[ok]
  rpk <- counts / (L / 1000)
  tpm <- t(t(rpk) / colSums(rpk)) * 1e6
  log10(tpm + 1)
}

collapse_by_mean <- function(mat) {
  base <- sub("\\..*", "", rownames(mat))
  s <- rowsum(mat, base)
  s / as.numeric(table(base)[rownames(s)])
}

## counts matrix straight out of a featureCounts table that carries `Length`
fc_to_log10tpm <- function(path) {
  d <- fread(path, data.table = FALSE)
  rownames(d) <- d$Geneid
  len <- setNames(d$Length, sub("\\..*", "", d$Geneid))
  sample_cols <- setdiff(colnames(d),
                         c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
  cnt <- as.matrix(d[, sample_cols])
  rownames(cnt) <- sub("\\..*", "", rownames(d))
  rpk <- cnt / (len / 1000)
  log10(t(t(rpk) / colSums(rpk)) * 1e6 + 1)
}

## ── GTEx ────────────────────────────────────────────────────────────────────
cat("GTEx...\n")
gtex_tpm  <- fread(GTEX_TPM_FILE, data.table = FALSE)
gtex_meta <- fread(GTEX_META_FILE, data.table = FALSE)
fc_cols   <- intersect(gtex_meta$SAMPID[gtex_meta$SMTSD == GTEX_TISSUE],
                       colnames(gtex_tpm))
gtex_mat  <- as.matrix(gtex_tpm[, fc_cols, drop = FALSE])
rownames(gtex_mat) <- gtex_tpm$Name
gtex_log10tpm <- log10(collapse_by_mean(gtex_mat) + 1)
cat("  ", nrow(gtex_log10tpm), "genes x", ncol(gtex_log10tpm), "samples\n")
rm(gtex_tpm, gtex_mat); gc()

## ── MSBB / ROSMAP ───────────────────────────────────────────────────────────
cat("MSBB...\n");   msbb_log10tpm   <- fc_to_log10tpm(MSBB_FILE)
cat("  ", dim(msbb_log10tpm), "\n")
cat("ROSMAP...\n"); rosmap_log10tpm <- fc_to_log10tpm(ROSMAP_FILE)
cat("  ", dim(rosmap_log10tpm), "\n")

## ── LBP living + postmortem ─────────────────────────────────────────────────
cat("LBP (GTF import is the slow part)...\n")
len30   <- load_gene_lengths(GTF_V30)
lbp_raw <- readRDS(LBP_FILE)

metadata <- lbp_raw$covariates
pt_idx   <- grepl("^PT-", metadata$IID_ISMMS)     # PT- prefix marks living donors

liv_counts <- lbp_raw$featureCounts[,  pt_idx]
pm_counts  <- lbp_raw$featureCounts[, !pt_idx]
cat("  living:", ncol(liv_counts), " postmortem:", ncol(pm_counts), "\n")

lbp_liv_log10tpm <- counts_to_log10tpm(liv_counts, len30)
lbp_pm_log10tpm  <- counts_to_log10tpm(pm_counts,  len30)

write.table(lbp_liv_log10tpm, file.path(BASE, "lbp_liv_log10tpm.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
rm(lbp_raw); gc()

## ── common universe + expression filter ─────────────────────────────────────
mats <- list(gtex    = gtex_log10tpm,
             lbp_liv = lbp_liv_log10tpm,
             lbp_pm  = lbp_pm_log10tpm,
             rosmap  = rosmap_log10tpm,
             msbb    = msbb_log10tpm)

common    <- Reduce(intersect, lapply(mats, rownames))
## The inverse has to use the same log base as above. With log10 data and a 2^
## inverse, the TPM threshold would be wrong and `genes` would change size.
gtex_lin  <- 10^gtex_log10tpm - 1                      # back to linear TPM
expressed <- rownames(gtex_lin)[rowMeans(gtex_lin >= TPM_MIN) >= FRAC_SAMPLES]
genes     <- intersect(common, expressed)
cat("\ncommon genes:", length(common),
    "| after expr filter:", length(genes), "\n")

mats       <- lapply(mats, function(m) m[genes, , drop = FALSE])
gene_means <- lapply(mats, rowMeans)
print(sapply(mats, dim))

## ── save ────────────────────────────────────────────────────────────────────
saveRDS(mats,         file.path(OUT_DIR, "mats.rds"))
saveRDS(genes,        file.path(OUT_DIR, "genes.rds"))
saveRDS(common,       file.path(OUT_DIR, "common.rds"))
saveRDS(gene_means,   file.path(OUT_DIR, "gene_means.rds"))
saveRDS(gtex_log10tpm, file.path(OUT_DIR, "gtex_log10tpm_full.rds"))  # for sweeps

cat("\nSaved to", OUT_DIR, "\n")
cat("Next: 02_expression_blood.R\n")