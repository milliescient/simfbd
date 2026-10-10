#!/usr/bin/env Rscript
# ESS of the rate parameters against record size, for an existing run.
#   ess_vs_taxa.R <run-dir> <batch-dir>
suppressMessages(library(posterior))
a <- commandArgs(TRUE); run <- a[1]; batch <- a[2]
fs <- list.files(file.path(run, "output"), pattern = "^rep_.*\\.log$", full.names = TRUE)
rows <- list()
for (f in fs) {
  n <- as.integer(sub(".*rep_(\\d+)\\.log$", "\\1", f))
  d <- try(read.delim(f, check.names = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || nrow(d) < 50) next
  keep <- grep("^(lambda|mu|psi)\\[", names(d), value = TRUE)
  if (!length(keep)) next
  burn <- seq_len(floor(nrow(d) / 4))
  e <- suppressWarnings(min(sapply(keep, function(k) ess_basic(d[[k]][-burn]))))
  tf <- file.path(batch, "specimens_complete", sprintf("taxa_%d.tsv", n))
  if (!file.exists(tf)) next
  tt <- read.delim(tf)
  rows[[length(rows) + 1]] <- data.frame(rep = n, taxa = length(unique(tt$taxon)),
                                         occ = nrow(tt), ess = e, gens = max(d$Iteration))
}
r <- do.call(rbind, rows)
cat(sprintf("\n  %s\n  %d chains, %.0f generations\n", basename(run), nrow(r), median(r$gens)))
cat(sprintf("  min ESS over rate parameters: median %.0f, p10 %.0f, %% below 200: %.0f%%\n",
            median(r$ess), quantile(r$ess, .1), 100 * mean(r$ess < 200)))
r$bin <- cut(r$taxa, breaks = unique(quantile(r$taxa, c(0, .33, .66, 1))), include.lowest = TRUE)
agg <- aggregate(cbind(ess, occ) ~ bin, r, median)
agg$n <- as.vector(table(r$bin))
cat("\n  taxa bin        n   median occ   median min-ESS\n")
for (i in seq_len(nrow(agg)))
  cat(sprintf("  %-14s %3d %12.0f %15.0f\n", agg$bin[i], agg$n[i], agg$occ[i], agg$ess[i]))
cat(sprintf("\n  correlation of min-ESS with taxa: %.2f;  with occurrences: %.2f\n",
            cor(r$ess, r$taxa, method = "spearman"), cor(r$ess, r$occ, method = "spearman")))
