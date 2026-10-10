#!/usr/bin/env Rscript
# Why a random-walk prior on skyline rates cannot be simulated from.
#
# Expected clade size is exponential in the log rates, so its upper tail is governed by the
# prior's marginal variance on them. A random walk's grows with the interval index; a
# stationary process's does not. Everything below measures that difference, and what it
# costs the simulator, which must reject draws it cannot build.
#
#   rate_prior_tails.R [out.tsv]        key<TAB>value, for the notebook page

set.seed(20261010)
N      <- as.integer(Sys.getenv("NDRAW", "20000"))
AMIN   <- 15; AMAX <- 21
LM     <- -1.5          # log lambda anchor mean
ANCH   <- 0.4           # anchor sd, shared by both priors
SIGMA  <- 0.4           # random-walk step sd
TAU    <- 0.5           # stationary sd
PHI    <- 0.5           # lag-1 autocorrelation
TANCH  <- 0.2; TSIG <- 0.2; TTAU <- 0.3   # the same, for log turnover

out <- new.env(); put <- function(k, v) assign(k, v, envir = out)

expected_lineages <- function(lambda, mu, age, shifts) {
  if (any(!is.finite(lambda)) || any(!is.finite(mu))) return(Inf)
  edges <- c(shifts, age); total <- 1; n <- 1
  for (j in seq_along(lambda)) {
    dt <- edges[j+1] - edges[j]; r <- lambda[j] - mu[j]
    if (!is.finite(r)) return(Inf)
    total <- total + lambda[j]*n*(if (abs(r) < 1e-12) dt else expm1(r*dt)/r)
    n <- n*exp(r*dt); if (!is.finite(total)) return(Inf)
  }
  total
}
rw <- function(m, anch, sd, n) cumsum(c(rnorm(1, m, anch), rnorm(n-1, 0, sd)))
ar <- function(m0, anch, tau, phi, n) {
  m <- rnorm(1, m0, anch); x <- numeric(n); x[1] <- rnorm(1, m, tau)
  s <- tau*sqrt(1-phi^2)
  if (n > 1) for (i in 2:n) x[i] <- m + phi*(x[i-1]-m) + rnorm(1, 0, s)
  x
}
hs <- function(m, anch, zeta, n)                       # Magee et al. 2020 horseshoe field
  cumsum(c(rnorm(1, m, anch), rnorm(n-1, 0, abs(rcauchy(1,0,zeta))*abs(rcauchy(n-1,0,1)))))

# One batch of draws under a named prior, returning expected clade size and lambda spread.
draws <- function(kind, nint, span = 12) {
  w <- span/nint; br <- w*seq_len(nint-1)
  en <- numeric(N); fold <- numeric(N)
  for (k in seq_len(N)) {
    if (kind == "gmrf_lm")      { a <- rw(LM,ANCH,SIGMA,nint); b <- rw(LM,ANCH,SIGMA,nint); lam <- rev(exp(a)); mu <- rev(exp(b)) }
    else if (kind == "gmrf")    { a <- rw(LM,ANCH,SIGMA,nint); t <- rw(0,TANCH,TSIG,nint);  lam <- rev(exp(a)); mu <- rev(exp(a+t)) }
    else if (kind == "hsmrf")   { a <- hs(LM,ANCH,0.2,nint);   t <- hs(0,TANCH,0.1,nint);   lam <- rev(exp(a)); mu <- rev(exp(a+t)) }
    else                        { a <- ar(LM,ANCH,TAU,PHI,nint); t <- ar(0,TANCH,TTAU,PHI,nint); lam <- rev(exp(a)); mu <- rev(exp(a+t)) }
    age <- runif(1, AMIN, AMAX)
    en[k] <- expected_lineages(lam, mu, age, c(0, age - rev(br)))
    fold[k] <- if (all(is.finite(lam)) && all(lam > 0)) max(lam)/min(lam) else NA
  }
  list(en = en, fold = fold)
}

# 1. rejection against interval count, the scaling the argument rests on
for (nint in c(4, 8, 16, 32)) for (kind in c("gmrf", "ar1")) {
  d <- draws(kind, nint)
  put(sprintf("rej_%s_%d", kind, nint), sprintf("%.2f", 100*mean(d$en > 400)))
}
# 2. the tail itself, at the design point
for (kind in c("gmrf", "ar1")) {
  d <- draws(kind, 8); e <- d$en[is.finite(d$en)]
  for (q in c(50, 90, 99, 99.9))
    put(sprintf("q%s_%s", sub("[.]", "p", q), kind), signif(quantile(e, q/100), 3))
  put(sprintf("fold_%s", kind), sprintf("%.1f", median(d$fold, na.rm = TRUE)))
  for (cap in c(4e2, 1e4, 1e8))
    put(sprintf("cap_%s_%g", kind, cap), sprintf("%.2f", 100*mean(d$en > cap)))
}
# 3. what the turnover reparameterisation and the horseshoe do, at 8 intervals
for (kind in c("gmrf_lm", "hsmrf")) {
  d <- draws(kind, 8)
  put(sprintf("rej_%s_8", kind), sprintf("%.2f", 100*mean(d$en > 400)))
  put(sprintf("fold_%s", kind), sprintf("%.1f", median(d$fold, na.rm = TRUE)))
}
put("ndraw", N)

keys <- sort(ls(out))
lines <- vapply(keys, function(k) paste0(k, "\t", get(k, envir = out)), "")
path <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ""
if (nzchar(path)) writeLines(lines, path) else writeLines(lines)
