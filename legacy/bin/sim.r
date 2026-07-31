# Parameterized SBC simulator for the dnFBDRMatrix model set (FossilSim/TreeSim).
# Reads config from environment variables (exported by run_sbc.sh). Draws rates
# from the SAME priors infer.Rev uses (valid SBC), simulates with the validated
# FossilSim tools, and writes data in the form the chosen SIM expects.
#
# Usage:  Rscript sim_sbc.r <rep> <seed>     (config via env: SIM,NINTERVALS,...)
#
# NOTE (limitation, see MISSPECIFICATIONS.md #4): TreeSim has no age-conditioned
# skyline simulator, so we use sim.rateshift.taxa (taxa-conditioned) + a time
# guard to skip mu>lambda draws that cannot reach the extant-tip target. This
# mildly conditions the deepest-interval rates (lambda in the oldest interval) on
# survival -> a known residual SBC wrinkle we cannot remove without an
# age-conditioned skyline simulator.
suppressWarnings(suppressMessages({library(TreeSim); library(FossilSim); library(ape); library(dplyr)}))

a       <- commandArgs(trailingOnly = TRUE)
rep     <- as.integer(a[1]);  seed <- as.integer(a[2])
geti    <- function(k, d) { v <- Sys.getenv(k); if (v=="") d else v }
SIM     <- geti("SIM","complete")
NINTERVALS <- as.integer(geti("NINTERVALS","3"))
IW     <- as.numeric(geti("INTERVAL_WIDTH","4"))
RHO     <- as.numeric(geti("RHO","1"))
LMEAN<-as.numeric(geti("LMEAN","-2")); LSD<-as.numeric(geti("LSD","0.5"))
PMEAN<-as.numeric(geti("PMEAN","-0.125")); PSD<-as.numeric(geti("PSD","0.5"))
AGE_MIN<-as.numeric(geti("AGE_MIN","5")); AGE_MAX<-as.numeric(geti("AGE_MAX","15"))
NMIN    <- as.integer(geti("NMIN","5"))
OUT     <- geti("DATADIR", file.path(geti("ROOT","."), "data"))
set.seed(seed)

ni  <- NINTERVALS
bps <- if (ni>1) IW*(1:(ni-1)) else numeric(0)   # breakpoints, youngest-first
dir.create(file.path(OUT,"specimens"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"times"),     recursive=TRUE, showWarnings=FALSE)

draw_rates <- function() list(lambda=exp(rnorm(ni,LMEAN,LSD)),
                              mu    =exp(rnorm(ni,LMEAN,LSD)),
                              psi   =exp(rnorm(ni,PMEAN,PSD)))

sim_one <- function() {
  r <- draw_rates()
  if (ni > 1) {
    tr <- tryCatch({ setTimeLimit(elapsed=3, transient=TRUE)
      sim.rateshift.taxa(n=NMIN, numbsim=1, lambda=r$lambda, mu=r$mu,
                         frac=rep(1,ni), times=c(0,bps), complete=TRUE)[[1]] },
      error=function(e) NULL)
    setTimeLimit()
  } else {
    age <- runif(1, AGE_MIN, AGE_MAX)
    tr  <- tryCatch(sim.bd.age(age, 1, r$lambda, r$mu, complete=TRUE)[[1]], error=function(e) NULL)
  }
  if (is.null(tr) || !inherits(tr,"phylo")) return(NULL)
  origin <- tr$root.edge + max(node.depth.edgelength(tr))
  tx <- tryCatch(sim.taxonomy(tr, beta=0), error=function(e) NULL)   # budding
  if (is.null(tx)) return(NULL)
  ia  <- c(0, bps, ceiling(origin)+1)                # interval boundaries for fossils
  exact <- (SIM=="complete")
  f <- tryCatch(sim.fossils.intervals(taxonomy=tx, interval.ages=ia,
                  rates=r$psi, use.exact.times=exact), error=function(e) NULL)
  if (is.null(f) || nrow(f)==0) return(NULL)
  st <- tx %>% group_by(sp) %>% summarise(ext=as.integer(any(round(end,8)==0)), .groups="drop")
  fsp <- sort(unique(f$sp))
  if (length(fsp) < NMIN) return(NULL)
  list(r=r, origin=origin, f=f, st=st, fsp=fsp)
}

s <- NULL
for (try in 1:400) { s <- sim_one(); if (!is.null(s)) break }
if (is.null(s)) stop("sim failed to produce data after 400 tries")

# ---- write data in the form the SIM expects ------------------------------
recs <- data.frame(taxon=character(), min_age=numeric(), max_age=numeric(), status=character())
status_of <- function(sp) ifelse(s$st$ext[s$st$sp==sp]==1, "extant", "extinct")
for (sp in s$fsp) {
  sub <- s$f[s$f$sp==sp, ]
  stat <- status_of(sp)
  if (SIM=="complete") {
    for (h in sub$hmin) recs <- rbind(recs, data.frame(taxon=paste0("t",sp), min_age=h, max_age=h, status=stat))
  } else if (SIM=="incomplete") {
    # oldest + youngest occurrence, as their binned intervals (drop intermediate count)
    oi <- which.max(sub$hmin); yi <- which.min(sub$hmin)
    rng <- unique(data.frame(min_age=c(sub$hmin[oi],sub$hmin[yi]), max_age=c(sub$hmax[oi],sub$hmax[yi])))
    for (j in 1:nrow(rng)) recs <- rbind(recs, data.frame(taxon=paste0("t",sp),
                                  min_age=rng$min_age[j], max_age=rng$max_age[j], status=stat))
  } else stop(paste("SIM not implemented in sim_sbc.r:", SIM))
}
write.table(recs, file.path(OUT,"specimens",paste0("taxa_",rep,".tsv")),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
if (ni>1) write.table(t(bps), file.path(OUT,"times",paste0("times_",rep,".tsv")),
                      row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# ---- append true values (youngest-first; matches rb index -> DIRECT mapping)
tvrow <- c(s$r$lambda, s$r$mu, s$r$psi, s$origin, length(s$fsp))
cat(paste(c(rep, sprintf("%.10g", tvrow)), collapse="\t"), "\n",
    file=file.path(OUT, paste0("truevals_",rep,".tsv")), sep="")
