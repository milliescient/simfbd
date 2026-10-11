# David Cerny's SBC diagnostics (rank hist / ECDF / ECDF-diff / coverage /
# true-vs-est) via the `SBC` package, adapted from skyfbdr_SBC.Rmd and generalised
# to N skyline intervals, with the interval-index convention as an explicit flag.
#
#   MAPPING = "direct"   : true index k == rb index k   (simfbd: youngest-first)
#   MAPPING = "reversed" : true index k == rb index (n+1-k)  (paleobuddy: oldest-first)
#
# This is the fix for the bug in the original skyfbdr_SBC.Rmd, which used the
# DIRECT mapping on paleobuddy data that needs REVERSED (see MISSPECIFICATIONS.md #1).
#
# Usage:  Rscript sbc_diagnostics.R <datadir> <outdir> <direct|reversed>
suppressMessages({library(SBC); library(posterior); library(gtools); library(ggplot2)})
a <- commandArgs(trailingOnly=TRUE)
datadir <- ifelse(length(a)>=1, a[1], "data")
outdir  <- ifelse(length(a)>=2, a[2], "output")
MAPPING <- ifelse(length(a)>=3, a[3], "direct")

logs <- mixedsort(list.files(outdir, pattern="fbd_.*\\.log", full.names=TRUE))
draws <- lapply(logs, function(x) {
  d <- read.table(x, header=TRUE)
  # rb columns lambda.1. -> lambda[1]
  colnames(d) <- gsub("\\.$", "]", gsub("\\.(?=[0-9])", "[", colnames(d), perl=TRUE))
  cols <- grep("^(lambda|mu|psi)\\[", colnames(d), value=TRUE)
  posterior::as_draws_matrix(d[, cols])
})

tv <- read.table(file.path(datadir, "true_vals.tsv"), header=TRUE)
colnames(tv) <- gsub("(\\D+)(\\d+)$", "\\1[\\2]", colnames(tv))
rates <- c("lambda","mu","psi")
ni <- sum(grepl("^lambda\\[", colnames(tv)))

# apply mapping: reorder the TRUE columns so they line up with rb's draw columns
if (MAPPING == "reversed") {
  for (r in rates) {
    cols <- paste0(r,"[",1:ni,"]")
    if (all(cols %in% colnames(tv))) tv[, cols] <- tv[, rev(cols)]
  }
}
sel <- unlist(lapply(rates, function(r) paste0(r,"[",1:ni,"]")))
true_vals <- posterior::as_draws_matrix(as.matrix(tv[, sel]))

res_raw <- SBC_results(stats=data.frame(sim_id=seq_along(draws)), fits=draws,
                       backend_diagnostics=NULL, default_diagnostics=data.frame(),
                       outputs=NULL, messages=NULL, warnings=NULL,
                       errors=vector("list", length(draws)))
ds  <- SBC_datasets(true_vals, generated=vector("list", length(draws)))
res <- recompute_SBC_statistics(res_raw, ds, backend=NULL)

pr <- file.path(outdir, paste0("sbc_", MAPPING, "_"))
ggsave(paste0(pr,"rank_hist.png"), plot_rank_hist(res),       width=12, height=8, dpi=120)
ggsave(paste0(pr,"ecdf_diff.png"), plot_ecdf_diff(res)+theme(legend.position="bottom"), width=12, height=8, dpi=120)
ggsave(paste0(pr,"coverage.png"),  plot_coverage(res),        width=12, height=8, dpi=120)
ggsave(paste0(pr,"true_vs_est.png"),plot_sim_estimated(res, alpha=0.5), width=12, height=8, dpi=120)
cat("SBC diagnostics written with prefix", pr, "(mapping:", MAPPING, ", n=", length(draws), ")\n")
