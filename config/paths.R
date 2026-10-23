## Reads back the paths exported by config/paths.sh.
## Run R scripts from the repo root, after `source config/paths.sh`.

need <- function(v) {
  x <- Sys.getenv(v)
  if (!nzchar(x)) stop(v, " is not set; run `source config/paths.sh` first", call. = FALSE)
  x
}

BASE          <- need("PROJECT_DIR")
PLOT_DIR      <- Sys.getenv("PLOT_DIR", file.path(BASE, "figures"))
FANTOM5_MEANS <- Sys.getenv("FANTOM5_MEANS", file.path(BASE, "fantom5_cage_gene_means.csv"))
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
