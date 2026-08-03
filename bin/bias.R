#!/usr/bin/env Rscript
# Bias of dnBDS against dnFBDRP on identical data, binned by the fraction of simulated
# lineages that were sampled. dnBDS assumes complete lineage sampling, so the unsampled
# fraction is the misspecification under measurement. Batches differ only in the origin
# prior, which is the lever on tree size and hence on how many lineages go unsampled.
#
#   Rscript bin/bias.R sims/pilot sims/age20 sims/age26 sims/age32

suppressMessages(library(posterior))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("usage: bias.R <batch dir> [<batch dir> ...]")

# which two arms to pair. FBD_RUN swaps in the PyRate-prior arm without touching the rest,
# since dnBDS has no origin and its arm is shared.
RUNS <- c(bds = "bds_complete", fbd = Sys.getenv("FBD_RUN", "fbdrp_complete"))
RATES <- c("lambda", "mu", "psi")
BURNIN <- 0.25
# Chains on the larger trees occasionally freeze, and a frozen chain reports whatever
# state it stalled in as a posterior mean. A replicate counts only if BOTH of its chains
# cleared this, so the pair is always the same data fitted two ways.
ESS_MIN <- 100

manifest <- function(dir, key) {
  m <- read.table(file.path(dir, "manifest.tsv"), sep = "\t", header = FALSE,
                  col.names = c("k", "v"), stringsAsFactors = FALSE)
  m$v[m$k == key]
}

# posterior mean of each rate in each interval, and the worst effective sample size over
# the nine of them. RevBayes numbers intervals opposite to the simulator, so the columns
# come back in simulator order.
post_means <- function(run_dir, rep, ni) {
  f <- file.path(run_dir, "output", sprintf("rep_%d.log", rep))
  if (!file.exists(f)) return(NULL)
  d <- try(read.table(f, header = TRUE, sep = "\t", check.names = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || nrow(d) < 50) return(NULL)
  d <- d[floor(nrow(d) * BURNIN):nrow(d), ]
  cols <- as.vector(outer(RATES, seq_len(ni), function(p, i) sprintf("%s[%d]", p, i)))
  if (!all(cols %in% colnames(d))) return(NULL)
  m <- sapply(RATES, function(p) rev(colMeans(d[, sprintf("%s[%d]", p, seq_len(ni)), drop = FALSE])))
  if (any(is.na(m))) return(NULL)
  list(means = m,                              # rows = simulator interval, cols = rate
       ess = min(apply(d[, cols], 2, posterior::ess_basic)))
}

# A replicate that never started is not missing at random: it failed because of what its
# own data looks like. Keep those separate from the ones a chain merely mixed badly on.
init_failed <- function(run_dir) {
  f <- file.path(run_dir, "failures.log")
  if (!file.exists(f) || file.size(f) == 0) return(integer(0))
  as.integer(sub(" .*", "", readLines(f)))
}

rows <- list()
excluded <- list()
for (batch in args) {
  batch <- sub("/$", "", batch)
  ni <- as.integer(manifest(batch, "NINTERVALS"))
  nreps <- as.integer(manifest(batch, "NREPS"))
  age_min <- as.numeric(manifest(batch, "AGE_MIN"))
  age_max <- as.numeric(manifest(batch, "AGE_MAX"))
  nums <- read.table(file.path(batch, "nums.tsv"), header = TRUE)
  tv <- read.table(file.path(batch, "true_vals.tsv"), header = TRUE)
  failed <- lapply(RUNS, function(r) init_failed(file.path(batch, "runs", r)))

  for (rep in seq_len(nreps)) {
    fits <- lapply(RUNS, function(r) post_means(file.path(batch, "runs", r), rep, ni))
    if (any(sapply(fits, is.null))) {
      why <- if (any(sapply(failed, function(v) rep %in% v))) "init" else "no log"
      excluded[[length(excluded) + 1]] <- data.frame(
        batch = basename(batch), rep = rep, why = why,
        models = paste(names(RUNS)[sapply(fits, is.null)], collapse = "+"),
        n_sp = nums$n_sp[rep], n_unsampled = nums$n_sp[rep] - nums$n_sampled[rep],
        stringsAsFactors = FALSE)
      next
    }
    truth <- sapply(RATES, function(p) as.numeric(tv[rep, paste0(p, seq_len(ni))]))

    row <- data.frame(batch = basename(batch), rep = rep,
                      origin = mean(c(age_min, age_max)), age = tv$age[rep],
                      n_sp = nums$n_sp[rep], n_sampled = nums$n_sampled[rep],
                      n_unsampled = nums$n_sp[rep] - nums$n_sampled[rep],
                      perc_sampled = nums$perc_sampled[rep],
                      ess = min(sapply(fits, function(f) f$ess)),
                      stringsAsFactors = FALSE)
    for (m in names(RUNS)) for (p in RATES) {
      # interval-averaged relative bias, plus each interval on its own
      fm <- fits[[m]]$means
      row[[paste0(m, "_", p)]] <- mean((fm[, p] - truth[, p]) / truth[, p])
      for (i in seq_len(ni))
        row[[sprintf("%s_%s%d", m, p, i)]] <- (fm[i, p] - truth[i, p]) / truth[i, p]
    }
    rows[[length(rows) + 1]] <- row
  }
}

D <- do.call(rbind, rows)
if (is.null(D)) stop("no replicate had both fits")
saveRDS(D, "bias.rds")

dropped <- D[D$ess < ESS_MIN, ]
D <- D[D$ess >= ESS_MIN, ]
cat(sprintf("%d replicates with both fits, from %d batches\n", nrow(D), length(unique(D$batch))))
cat(sprintf("%d dropped for min ESS below %d (median tree %.0f lineages, vs %.0f kept)\n",
            nrow(dropped), ESS_MIN,
            if (nrow(dropped)) median(dropped$n_sp) else 0, median(D$n_sp)))

X <- do.call(rbind, excluded)
if (!is.null(X)) {
  ini <- X[X$why == "init", ]
  cat(sprintf("%d never produced a pair: %d failed to initialise, %d had no log yet\n",
              nrow(X), nrow(ini), sum(X$why == "no log")))
  if (nrow(ini)) {
    # the comparison that matters: init failure concentrating in the informative stratum
    cat(sprintf("  init failures: %s\n",
                paste(sprintf("%s rep %d (%s, %d lineages, %d unsampled)", ini$batch, ini$rep,
                              ini$models, ini$n_sp, ini$n_unsampled), collapse = "; ")))
    cat(sprintf("  %.0f%% of them carry an unsampled lineage, against %.0f%% of those kept\n",
                100 * mean(ini$n_unsampled > 0), 100 * mean(D$n_unsampled > 0)))
  }
}
cat("lineages simulated: ", paste(range(D$n_sp), collapse = "-"),
    "   sampled: ", paste(range(D$n_sampled), collapse = "-"),
    "   fully sampled: ", sprintf("%.0f%%", 100 * mean(D$perc_sampled == 1)), "\n\n", sep = "")

by_batch <- function() {
  cat("by batch (origin prior is the only difference between them)\n")
  cat(sprintf("%-8s %-12s %5s %7s %7s %8s\n", "batch", "origin", "reps", "med n_sp", "max n_sp", "%full"))
  for (b in unique(D$batch)) {
    s <- D[D$batch == b, ]
    cat(sprintf("%-8s U(%4.0f,%4.0f) %5d %7.0f %8d %7.0f%%\n", b,
                s$origin[1] - 3, s$origin[1] + 3, nrow(s), median(s$n_sp),
                max(s$n_sp), 100 * mean(s$perc_sampled == 1)))
  }
}

# Both models share the data and the prior, so on data-poor replicates both are pulled to
# the same prior mean and the raw bias measures that pull, not the model. The paired
# difference cancels it and isolates the complete-lineage-sampling assumption.
bias_table <- function(bins, labels, var, title) {
  cat("\n", title, "\n", sep = "")
  cat(sprintf("%-10s %5s %6s %6s | %-23s | %-23s\n", "", "reps", "med n", "minESS",
              "BDS-FBDRP paired diff", "raw bias  BDS / FBDRP"))
  g <- cut(D[[var]], breaks = bins, labels = labels, include.lowest = TRUE)
  for (lv in levels(g)) {
    s <- D[which(g == lv), ]
    if (nrow(s) == 0) next
    d <- sapply(RATES, function(p) s[[paste0("bds_", p)]] - s[[paste0("fbd_", p)]])
    d <- matrix(d, nrow = nrow(s))
    se <- apply(d, 2, function(x) sd(x) / sqrt(length(x)))
    cat(sprintf("%-10s %5d %6.0f %6.0f | %s | %s\n", lv, nrow(s), median(s$n_sp), min(s$ess),
                paste(sprintf("%+.2f%s", colMeans(d),
                              ifelse(abs(colMeans(d)) > 2 * se, "*", " ")), collapse = " "),
                paste(sprintf("%+.2f/%+.2f", c(mean(s$bds_lambda), mean(s$bds_mu), mean(s$bds_psi)),
                              c(mean(s$fbd_lambda), mean(s$fbd_mu), mean(s$fbd_psi))),
                      collapse = " ")))
  }
  cat("  * = mean paired difference exceeds twice its standard error\n")
}

by_batch()
bias_table(c(0, 0.5, 0.75, 0.9, 0.999, 1), c("<0.50", "0.50-0.75", "0.75-0.90", "0.90-1.00", "1.00 (all)"),
           "perc_sampled", "binned by fraction of lineages sampled")
bias_table(c(-1, 0, 1, 3, 10, 1e6), c("0", "1", "2-3", "4-10", "11+"),
           "n_unsampled", "binned by count of unsampled lineages")

# The two densities differ even on a fully sampled record, since dnBDS has no survival term,
# so the intercept is that standing gap and the slope is what the unsampled lineages add.
cat("\npaired difference regressed on the unsampled count\n")
cat(sprintf("%-7s %18s %18s\n", "", "intercept (0 unsamp)", "slope per lineage"))
for (p in RATES) {
  y <- D[[paste0("bds_", p)]] - D[[paste0("fbd_", p)]]
  co <- summary(lm(y ~ D$n_unsampled))$coefficients
  cat(sprintf("%-7s %+8.3f +/- %-5.3f %+8.4f +/- %-6.4f%s\n", p,
              co[1, 1], co[1, 2], co[2, 1], co[2, 2],
              if (abs(co[2, 3]) > 2) "  *" else ""))
}
cat("  * = slope exceeds twice its standard error\n")

cat("\nSpearman rho against the unsampled count (paired, same data both models)\n")
for (p in RATES)
  cat(sprintf("  %-7s BDS %+0.3f   FBDRP %+0.3f   BDS-FBDRP %+0.3f\n", p,
              cor(D$n_unsampled, D[[paste0("bds_", p)]], method = "spearman"),
              cor(D$n_unsampled, D[[paste0("fbd_", p)]], method = "spearman"),
              cor(D$n_unsampled, D[[paste0("bds_", p)]] - D[[paste0("fbd_", p)]],
                  method = "spearman")))

cat("\nwrote bias.rds\n")
