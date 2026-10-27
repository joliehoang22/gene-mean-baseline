## ============================================================================
## Paired t-tests behind Fig. 1 and Supp. Table 1: every model against every
## gene-mean baseline in each cohort (4 models x 12 cohort-baseline pairs = 48 tests).
## A model means its top-performing track in that cohort (highest median), as in Fig. 1.
## Each participant's samples are averaged first, so tests are at participant level
## (participant = column 1 of the cohort's sample table).
## ============================================================================
suppressMessages(library(data.table))
source("config/paths.R")   # BASE
OUTDIR <- file.path(BASE, "evaluation", "allcohorts")
METRIC <- "spearman"

d <- fread(file.path(OUTDIR, "arm_scores_allcohorts_qcexact.csv"))
d[, r := get(METRIC)]

MAPS <- list(`LBP-Brain`    = c("sample_mapping.tsv", TRUE),        `LBP-Blood`   = c("samples_blood.tsv", FALSE),
             `ROSMAP-Brain` = c("sample_mapping_rosmap.tsv", TRUE), `COVID-Blood` = c("samples_covid.tsv", FALSE))
PID <- lapply(MAPS, function(m) { x <- read.delim(file.path(BASE, m[1]), header = as.logical(m[2]), stringsAsFactors = FALSE)
                                  setNames(as.character(x[[1]]), x[[2]]) })

## two-sided p from t and df on the log scale, so p-values below 1e-308 are not printed as 0
log10_p <- function(t, df) (log(2) + pt(-abs(t), df, log.p = TRUE)) / log(10)
p_text  <- function(l10) { e <- floor(l10); m <- round(10^(l10 - e), 1)
                           e[m == 10] <- e[m == 10] + 1; m[m == 10] <- 1; sprintf("%.1fe%d", m, e) }

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
  tp <- t.test(pb, pa, paired = TRUE); l10 <- unname(log10_p(tp$statistic, tp$parameter))
  data.table(cohort = co, model = A$label, gene_mean = B$label,
             n_samples = nrow(m), n_participants = length(pa),
             mean_model = mean(m$r_a), mean_genemean = mean(m$r_b),
             median_model = median(m$r_a), median_genemean = median(m$r_b),
             t = unname(tp$statistic), df = unname(tp$parameter), log10_p = l10, p = p_text(l10))
}

res <- rbindlist(lapply(unique(d$cohort), function(co) {
  gms <- unique(d[cohort == co & model == "Gene-mean", arm])
  rbindlist(lapply(c("Enformer", "Basenji2", "ExPecto", "Xpresso"), function(mo)
    rbindlist(lapply(gms, function(g) run(co, mo, g)))))
}))
fwrite(res, file.path(OUTDIR, sprintf("paired_stats_%s.csv", METRIC)))
print(res[, .(cohort, model, gene_mean, n_participants, mean_model, mean_genemean, t, df, p)], digits = 3)
cat(sprintf("\n%d tests; largest p = %s\n", nrow(res), res[which.max(log10_p), p]))
