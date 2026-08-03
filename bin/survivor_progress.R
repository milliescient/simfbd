#!/usr/bin/env Rscript
# The survivors comparison, over whatever replicates both arms have finished.
#
# Paired on the replicate: each arm fits the same record, so the difference between them is
# the point mass and nothing else.
suppressMessages(library(posterior))
setwd("/research/phyloworks/simfbd")
BATCH <- "sims/survivors"
ESS_MIN <- 100
load(file.path(BATCH, "sim_list.RData"))

score <- function(arm) {
  out <- file.path(BATCH, "runs", arm, "output")
  res <- list()
  for (f in list.files(out, pattern = "^skyfbdr_[0-9]+\\.log$", full.names = TRUE)) {
    r <- as.integer(sub(".*_([0-9]+)\\.log$", "\\1", f))
    x <- try(read.table(f, header = TRUE, sep = "\t", check.names = FALSE), silent = TRUE)
    if (inherits(x, "try-error") || nrow(x) < 500) next
    x <- x[floor(nrow(x)*0.25):nrow(x), ]
    rc <- grep("^(lambda|mu|psi)", colnames(x), value = TRUE)
    if (min(sapply(rc, function(nm) ess_basic(x[[nm]]))) < ESS_MIN) next
    dd <- x[, grep("^dd\\[", colnames(x)), drop = FALSE]
    sp <- read.table(file.path(BATCH, "specimens_complete", sprintf("taxa_%d.tsv", r)),
                     header = TRUE, sep = "\t")
    ord <- sort(unique(as.character(sp$taxon)), method = "radix")
    if (length(ord) != ncol(dd)) next
    k <- as.integer(sub("^t", "", ord))
    res[[as.character(r)]] <- c(p = mean(colMeans(dd == 0)),
                                obs = mean(is.na(sims[[r]]$TE[k])))
  }
  res
}

a <- score("survivors_true"); b <- score("survivors_false")
both <- intersect(names(a), names(b))
cat(sprintf("\nfinished: survivors_true %d, survivors_false %d, paired %d\n",
            length(a), length(b), length(both)))
if (length(both) < 2) { cat("not enough paired replicates yet\n"); quit() }

A <- do.call(rbind, a[both]); B <- do.call(rbind, b[both])
cat(sprintf("%-16s %9s %10s %11s %8s\n", "arm", "mean p", "observed", "err", "z"))
for (nm in c("survivors_true", "survivors_false")) {
  M <- if (nm == "survivors_true") A else B
  e <- M[, "obs"] - M[, "p"]
  se <- sd(e)/sqrt(length(e))
  cat(sprintf("%-16s %9.3f %10.3f %+11.4f %8.2f\n",
              nm, mean(M[, "p"]), mean(M[, "obs"]), mean(e), mean(e)/se))
}
d <- (A[, "obs"] - A[, "p"]) - (B[, "obs"] - B[, "p"])
cat(sprintf("\npaired difference (true - false) %+.4f +/- %.4f  z %+.2f\n",
            mean(d), sd(d)/sqrt(length(d)), mean(d)/(sd(d)/sqrt(length(d)))))
cat("positive err = under-credits survival\n")
