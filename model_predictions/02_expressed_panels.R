## 02_expressed_panels.R
## Gene panel for each tissue: the MANE Select TSSs of expressed genes
## Uses the expressed-gene sets saved by gene_mean_baseline/01-02_expression_*.R (Null/genes.rds,
## Null/genes_blood.rds), so the expression filter is defined in one place.
source("config/paths.R")   # BASE, need()
mane <- read.table(file.path(BASE, "tss_hg38_mane_all_autosomal.bed"), sep = "\t",
                   col.names = c("chrom","tss","strand","gene_id","gene_name"),
                   stringsAsFactors = FALSE)
cat("MANE Select autosomal genes:", nrow(mane), "\n\n")

for (tis in c("brain","blood")) {
  gfile <- if (tis == "brain") "Null/genes.rds" else "Null/genes_blood.rds"
  genes <- readRDS(file.path(BASE, gfile))
  keep  <- mane[mane$gene_id %in% genes, ]
  out   <- file.path(BASE, sprintf("tss_hg38_mane_expressed_%s.bed", tis))
  write.table(keep, out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("%-5s expressed universe %5d | MANE and expressed %5d -> %s\n",
              tis, length(genes), nrow(keep), basename(out)))
}
