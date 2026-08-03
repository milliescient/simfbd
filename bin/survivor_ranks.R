#!/usr/bin/env Rscript
# What forbidding the point mass does to the RATES, in SBC ranks.
#
# Ranks, not point estimates: (est-true)/true divides by a lognormal draw that can be near
# zero, which reads as a large bias even when the posterior is centred. A calibrated rate has
# mean rank 0.5, with SE 1/sqrt(12n).
suppressMessages(library(posterior))
setwd("/research/phyloworks/simfbd")
BATCH <- "sims/survivors"
ESS_MIN <- 100
COLS <- c(sprintf("lambda[%d]", 1:3), sprintf("mu[%d]", 1:3), sprintf("psi[%d]", 1:3))

# A batch without true_vals_order predates the fix and holds the rates oldest-first, which is
# paleobuddy's order. Reading one as the other silently compares interval 1 against interval 3.
OLDEST_FIRST <- !any(grepl("true_vals_order",
                           readLines(file.path(BATCH, "manifest.tsv"), warn = FALSE)))
idx <- if (OLDEST_FIRST) 3:1 else 1:3
PARS <- c(paste0("lambda", idx), paste0("mu", idx), paste0("psi", idx))

truth <- read.table(file.path(BATCH, "true_vals.tsv"), header = TRUE, sep = "\t")

ranks_for <- function(arm) {
  out <- file.path(BATCH, "runs", arm, "output")
  res <- list()
  for (f in list.files(out, pattern = "^skyfbdr_[0-9]+\\.log$", full.names = TRUE)) {
    r <- as.integer(sub(".*_([0-9]+)\\.log$", "\\1", f))
    x <- try(read.table(f, header = TRUE, sep = "\t", check.names = FALSE), silent = TRUE)
    if (inherits(x, "try-error") || nrow(x) < 500) next
    x <- x[floor(nrow(x)*0.25):nrow(x), ]
    if (min(sapply(COLS, function(nm) ess_basic(x[[nm]]))) < ESS_MIN) next
    res[[as.character(r)]] <- vapply(seq_along(COLS),
      function(j) mean(x[[COLS[j]]] < truth[r, PARS[j]]), numeric(1))
  }
  res
}

a <- ranks_for("survivors_true"); b <- ranks_for("survivors_false")
both <- intersect(names(a), names(b))
cat(sprintf("\npaired replicates: %d (true %d, false %d)\n", length(both), length(a), length(b)))
if (length(both) < 5) { cat("not enough yet\n"); quit() }

A <- do.call(rbind, a[both]); B <- do.call(rbind, b[both])
n <- length(both); se <- 1/sqrt(12*n)
cat(sprintf("mean rank, 0.5 is calibrated, SE %.3f\n\n", se))
cat(sprintf("%-10s %18s %18s %12s\n", "", "survivors=TRUE", "survivors=FALSE", "paired diff"))
for (j in seq_along(PARS)) {
  ra <- mean(A[, j]); rb <- mean(B[, j])
  d <- A[, j] - B[, j]
  sd_ <- sd(d)/sqrt(n)
  cat(sprintf("%-10s %10.3f (z%+5.1f) %10.3f (z%+5.1f) %+7.3f (z%+5.1f)\n",
              COLS[j], ra, (ra-0.5)/se, rb, (rb-0.5)/se, mean(d),
              if (sd_ > 0) mean(d)/sd_ else 0))
}
cat("\nrank below 0.5 means the posterior sits above the truth\n")
