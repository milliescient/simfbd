#!/usr/bin/env Rscript
# Does one batch serve both comparisons? Per replicate, the fraction of lineages that
# leave no fossil, which is what the PyRate arm needs to bite, against psi, which is what
# the classical estimators need to be defined at all.
a <- commandArgs(TRUE); b <- a[1]
e <- new.env(); load(file.path(b, "sim_list.RData"), envir = e)
tv <- read.delim(file.path(b, "true_vals.tsv"))
pc <- grep("^psi", names(tv), value = TRUE)
rows <- list()
for (r in seq_along(e$sims)) {
  f <- file.path(b, "specimens_complete", sprintf("taxa_%d.tsv", r))
  if (!file.exists(f)) next
  s <- e$sims[[r]]; tx <- read.delim(f)
  seen <- rep(FALSE, length(s$TS))
  seen[as.integer(sub("^t", "", unique(as.character(tx$taxon))))] <- TRUE
  rows[[length(rows) + 1]] <- data.frame(
    rep = r, psi = mean(as.numeric(tv[r, pc])),
    lineages = length(s$TS), unsamp = 100 * mean(!seen))
}
x <- do.call(rbind, rows)
x <- x[order(x$psi), ]
k <- floor(nrow(x) / 4)
cat(sprintf("\n  %s: %d reps, %d lineages total\n", basename(b), nrow(x), sum(x$lineages)))
cat("  psi quartile    median psi   median lineages   %% lineages leaving no fossil\n")
for (i in 1:4) {
  s <- x[((i - 1) * k + 1):(if (i == 4) nrow(x) else i * k), ]
  cat(sprintf("  Q%d %-12s %10.3f %17.0f %28.1f\n", i, "",
              median(s$psi), median(s$lineages), median(s$unsamp)))
}
cat(sprintf("\n  overall: %.1f%% of all lineages leave no fossil\n",
            100 * sum(x$lineages * x$unsamp / 100) / sum(x$lineages)))
