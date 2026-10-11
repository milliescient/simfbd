#!/usr/bin/env Rscript
# How much rate variation through time does a given GMRF_SD actually produce?
a <- commandArgs(TRUE); b <- a[1]
tv <- read.delim(file.path(b, "true_vals.tsv"))
fold <- function(pre) {
  m <- as.matrix(tv[, grep(paste0("^", pre), names(tv))])
  apply(m, 1, function(r) max(r) / min(r))
}
l <- fold("lambda"); p <- fold("psi")
cat(sprintf("  %-14s reps %3d   lambda max/min across intervals: median %5.1fx  p90 %6.1fx   psi: median %5.1fx\n",
            basename(b), nrow(tv), median(l), quantile(l, .9), median(p)))
