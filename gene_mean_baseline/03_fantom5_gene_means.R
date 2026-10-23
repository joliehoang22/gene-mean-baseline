############################################################
# 03_fantom5_gene_means.R
#
# FANTOM5 CAGE gene-mean baseline. For each gene, the mean CAGE expression
# across the FANTOM5 libraries behind the Enformer/Basenji2 tracks we score
# (7 brain libraries in brain_tracks.txt, 4 whole-blood libraries in
# blood_tracks.txt). These libraries are part of the models' training data.
#
# Input:  FANTOM5 gene-level TPM table (hg38_fair+new_genes_phase1and2_tpm.osc.txt.gz)
# Output: $FANTOM5_MEANS, one row per gene symbol:
#         gene, fantom_brain_mean, fantom_blood_mean  (TPM; the scorer applies log2(x+1))
############################################################

suppressPackageStartupMessages(library(data.table))
source("config/paths.R")   # BASE, FANTOM5_MEANS, need()

f <- need("FANTOM5_CAGE")

# FANTOM5 library IDs (CNhs...) for the tracks in brain_tracks.txt / blood_tracks.txt
brain_ids <- c("CNhs10617","CNhs10647","CNhs12005","CNhs11796","CNhs12310","CNhs12338","CNhs13796")
blood_ids <- c("CNhs11675","CNhs11671","CNhs11761","CNhs11949")

# The .osc file opens with "##" metadata lines; count them so fread can skip.
con <- gzfile(f, "rt")
skip <- 0L
repeat {
  x <- readLines(con, n = 1)
  if (!length(x) || !startsWith(x, "##")) break
  skip <- skip + 1L
}
close(con)

d <- fread(cmd = paste("zcat", shQuote(f)), skip = skip)

# Column names embed the CNhs ID; each ID has to hit exactly one column.
get_cols <- function(x) sapply(x, function(id) {
  hit <- grep(id, names(d), value = TRUE)
  if (length(hit) != 1) stop(id, " matched ", length(hit), " columns")
  hit
})
brain_cols <- get_cols(brain_ids)
blood_cols <- get_cols(blood_ids)

out <- data.table(
  gene              = d[[1]],
  fantom_brain_mean = rowMeans(d[, ..brain_cols], na.rm = TRUE),
  fantom_blood_mean = rowMeans(d[, ..blood_cols], na.rm = TRUE)
)
# Expect 23,274 genes: 20,057 with a nonzero brain mean, 15,496 nonzero in blood.

fwrite(out, FANTOM5_MEANS)
cat("wrote", FANTOM5_MEANS, "--", nrow(out), "genes\n")
