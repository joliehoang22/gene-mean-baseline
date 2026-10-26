################################################################################
## 02_expression_blood.R
##
## Blood version of 01_expression_brain.R. Later scripts load the .rds files
## it writes.
##
## Three cohorts: GTEx whole blood, LBP LIV blood and the Mount Sinai COVID-19 Biobank.
##
## The counts from the COVID-19 Biobank were generated with GENCODE v30, the same
## annotation as LBP, so they go through the same counts -> TPM steps.
##
## Same scale as brain: log10(TPM + 1), version-stripped ENSG IDs. TPM is used
## for the reason given in 01_expression_brain.R.
##
## The expression filter is defined on GTEx whole blood.
##
## OUTPUTS (Null/)
##   mats_blood.rds        list(gtex_blood, lbp_blood, covid_blood), restricted to genes_blood
##   genes_blood.rds       expression-filtered common gene universe
################################################################################

library(data.table)
library(rtracklayer)

source("config/paths.R")   # BASE, need()
OUT_DIR <- file.path(BASE, "Null")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GTF_V30        <- need("GENCODE_V30_GTF")
GTEX_TPM_FILE  <- need("GTEX_TPM")
GTEX_META_FILE <- need("GTEX_SAMPLE_ATTR")
GTEX_TISSUE    <- "Whole Blood"
LBP_FILE       <- need("LBP_BLOOD_RDS")
COVID_FILE     <- need("COVID_COUNTS")
TPM_MIN      <- 1
FRAC_SAMPLES <- 0.25

## ── helpers (identical to the brain script) ─────────────────────────────────
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

## ── GTEx Whole Blood ────────────────────────────────────────────────────────
cat("GTEx Whole Blood...\n")
gtex_tpm  <- fread(GTEX_TPM_FILE, data.table = FALSE)
gtex_meta <- fread(GTEX_META_FILE, data.table = FALSE)
stopifnot(GTEX_TISSUE %in% gtex_meta$SMTSD)

wb_cols  <- intersect(gtex_meta$SAMPID[gtex_meta$SMTSD == GTEX_TISSUE],
                      colnames(gtex_tpm))
gtex_mat <- as.matrix(gtex_tpm[, wb_cols, drop = FALSE])
rownames(gtex_mat) <- gtex_tpm$Name
gtex_blood <- log10(collapse_by_mean(gtex_mat) + 1)
cat("  ", nrow(gtex_blood), "genes x", ncol(gtex_blood), "samples\n")
rm(gtex_tpm, gtex_mat); gc()

## ── LBP blood ───────────────────────────────────────────────────────────────
cat("LBP blood (GTF import is the slow part)...\n")
len30   <- load_gene_lengths(GTF_V30)
lbp_raw <- readRDS(LBP_FILE)

blood_cols <- grep("BLOOD", colnames(lbp_raw), value = TRUE)
cat("  blood samples:", length(blood_cols), "\n")
lbp_blood <- counts_to_log10tpm(lbp_raw[, blood_cols, drop = FALSE], len30)
cat("  ", dim(lbp_blood), "\n")
rm(lbp_raw); gc()

## ── Mount Sinai COVID-19 Biobank blood ──────────────────────────────────────
## Same GENCODE v30 annotation as LBP, so the same counts -> TPM path applies.
cat("COVID biobank blood...\n")
cov_raw <- fread(COVID_FILE, data.table = FALSE)
rownames(cov_raw) <- cov_raw[[1]]
cov_cnt <- as.matrix(cov_raw[, -1, drop = FALSE])

## stop if the matrix is not integer counts: a pre-normalized matrix cannot be
## converted to TPM, and nothing later would flag it
if (!all(cov_cnt == floor(cov_cnt), na.rm = TRUE))
  stop("COVID matrix is not integer counts; check what normalization was applied")
cat("  ", nrow(cov_cnt), "genes x", ncol(cov_cnt), "samples\n")
cat("   first gene id:", rownames(cov_cnt)[1], "\n")

covid_blood <- counts_to_log10tpm(cov_cnt, len30)
cat("  after length matching:", nrow(covid_blood), "genes\n")
if (nrow(covid_blood) < 0.5 * nrow(cov_cnt))
  warning("fewer than half of COVID genes matched v30 lengths; check the ID format")
rm(cov_raw, cov_cnt); gc()

## ── common universe + expression filter (defined on GTEx Whole Blood) ───────
mats <- list(gtex_blood  = gtex_blood,
             lbp_blood   = lbp_blood,
             covid_blood = covid_blood)

common      <- Reduce(intersect, lapply(mats, rownames))
gtex_lin    <- 10^gtex_blood - 1                      # same log base as above
expressed   <- rownames(gtex_lin)[rowMeans(gtex_lin >= TPM_MIN) >= FRAC_SAMPLES]
genes_blood <- intersect(common, expressed)
cat("\ncommon genes:", length(common),
    "| after expr filter:", length(genes_blood), "\n")

mats <- lapply(mats, function(m) m[genes_blood, , drop = FALSE])
print(sapply(mats, dim))

## ── save ────────────────────────────────────────────────────────────────────
saveRDS(mats,        file.path(OUT_DIR, "mats_blood.rds"))
saveRDS(genes_blood, file.path(OUT_DIR, "genes_blood.rds"))

cat("\nSaved to", OUT_DIR, "\n")
cat("Next: 03_fantom5_gene_means.R\n")
