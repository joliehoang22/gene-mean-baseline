## ============================================================================
## Paired tests quoted in the text: a model vs a gene-mean baseline, using the
## per-sample cross-gene correlations behind Fig. 1.
##   CMP  : the comparisons quoted in the text (add lines as needed)
##   grid : every model x every outside gene-mean baseline, per cohort
## A model name means its top-performing track in that cohort, as in Fig. 1.
## Main test: paired t-test at participant level (a participant's samples are
## averaged first; participant = column 1 of the cohort's sample table).
## ============================================================================
suppressMessages(library(data.table))
source("config/paths.R")   # BASE
OUTDIR <- file.path(BASE, "evaluation", "allcohorts")
METRIC <- "spearman"

CMP <- list(
  list(cohort = "LBP-Brain",    a = "Enformer", b = "Gene-mean: LBP PM"),
  list(cohort = "ROSMAP-Brain", a = "Enformer", b = "Gene-mean: LBP PM"),
  list(cohort = "LBP-Blood",    a = "Enformer", b = "Gene-mean: GTEx Blood"),
  list(cohort = "LBP-Blood",    a = "Basenji2", b = "Gene-mean: GTEx Blood"),
  list(cohort = "LBP-Blood",    a = "ExPecto",  b = "Gene-mean: GTEx Blood"),
  list(cohort = "LBP-Blood",    a = "Xpresso",  b = "Gene-mean: GTEx Blood")
)

d <- fread(file.path(OUTDIR, "arm_scores_allcohorts_qcexact.csv"))
d[, r := get(METRIC)]

MAPS <- list(`LBP-Brain`    = c("sample_mapping.tsv", TRUE),        `LBP-Blood`   = c("samples_blood.tsv", FALSE),
             `ROSMAP-Brain` = c("sample_mapping_rosmap.tsv", TRUE), `COVID-Blood` = c("samples_covid.tsv", FALSE))
PID <- lapply(MAPS, function(m) { x <- read.delim(file.path(BASE, m[1]), header = as.logical(m[2]), stringsAsFactors = FALSE)
                                  setNames(as.character(x[[1]]), x[[2]]) })

arm_values <- function(co, a) {
  x <- d[cohort == co]
  if (a %in% x$arm) return(list(v = x[arm == a, .(sample, r)], label = a))
  m <- x[model == a]
  best <- m[, .(med = median(r, na.rm = TRUE)), by = track][which.max(med), track]
  list(v = m[track == best, .(sample, r)], label = paste0(a, ": ", best))
}

run <- function(co, a, b) {
  A <- arm_values(co, a); B <- arm_values(co, b)
  m <- merge(A$v, B$v, by = "sample", suffixes = c("_a", "_b"))[is.finite(r_a) & is.finite(r_b)]
  pid <- unname(PID[[co]][m$sample])
  pa <- tapply(m$r_a, pid, mean); pb <- tapply(m$r_b, pid, mean)
  tp <- t.test(pb, pa, paired = TRUE)
  ts <- t.test(m$r_b, m$r_a, paired = TRUE)
  data.table(cohort = co, A = A$label, B = B$label,
             n_samples = nrow(m), n_participants = length(pa), participants_2plus = sum(table(pid) > 1),
             median_A = median(m$r_a), median_B = median(m$r_b), mean_A = mean(m$r_a), mean_B = mean(m$r_b),
             pct_samples_B_higher = 100 * mean(m$r_b > m$r_a),
             t_part = unname(tp$statistic), df_part = unname(tp$parameter), p_part = tp$p.value,
             t_samp = unname(ts$statistic), df_samp = unname(ts$parameter), p_samp = ts$p.value,
             p_wilcox_samp = suppressWarnings(wilcox.test(m$r_b, m$r_a, paired = TRUE, exact = FALSE)$p.value))
}

res <- rbindlist(lapply(CMP, function(x) run(x$cohort, x$a, x$b)))
fwrite(res, file.path(OUTDIR, sprintf("paired_stats_%s.csv", METRIC)))
print(res[, .(cohort, A, B, n_participants, median_A, median_B, t_part, df_part, p_part)], digits = 3)

grid <- rbindlist(lapply(unique(d$cohort), function(co) {
  gms <- grep("LOO", unique(d[cohort == co & model == "Gene-mean", arm]), value = TRUE, invert = TRUE)
  rbindlist(lapply(c("Enformer", "Basenji2", "ExPecto", "Xpresso"), function(mo)
    rbindlist(lapply(gms, function(g) run(co, mo, g)))))
}))
fwrite(grid, file.path(OUTDIR, sprintf("paired_stats_grid_%s.csv", METRIC)))
