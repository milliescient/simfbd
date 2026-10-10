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
g0 <- function(k) get(k, envir = out)

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
# how far the refusal rate moves over that range, which is what the page claims
for (kind in c("gmrf", "ar1")) {
  lo <- as.numeric(get(sprintf("rej_%s_4", kind), envir = out))
  hi <- as.numeric(get(sprintf("rej_%s_32", kind), envir = out))
  put(sprintf("rej_ratio_%s", kind), sprintf("%.1f", if (kind == "gmrf") hi/lo else lo/hi))
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

# 4. the closed form behind the scaling. log E[N] is the span times the AVERAGE net rate
# over the intervals, so what the screen bounds is a path average, not a single rate. For a
# walk that average does not concentrate as the grid is refined and for a stationary process
# it does, which is the whole asymmetry:
#
#   walk    Var(mean of x) = a^2 + s^2 (K-1)(2K-1) / 6K        -> grows like s^2 K / 3
#   AR(1)   Var(mean of x) = a^2 + (t^2/K^2)(K + 2 sum (K-k) phi^k)  -> a^2 + O(1/K)
#
# Both are exact, and checked against draws below.
vbar_rw <- function(K, a, s) a^2 + s^2 * (K-1) * (2*K-1) / (6*K)
vbar_ar <- function(K, a, t, phi) {
  k <- seq_len(K-1)
  a^2 + (t^2/K^2) * (K + 2*sum((K-k) * phi^k))
}
mc <- function(kind, K, n = 20000) {
  m <- replicate(n, mean(if (kind == "gmrf") rw(LM, ANCH, SIGMA, K)
                         else ar(LM, ANCH, TAU, PHI, K)))
  var(m)
}
vb <- list()
for (K in c(4, 8, 16, 32)) {
  ex <- c(gmrf = vbar_rw(K, ANCH, SIGMA), ar1 = vbar_ar(K, ANCH, TAU, PHI))
  for (kind in c("gmrf", "ar1")) {
    put(sprintf("vbar_%s_%d", kind, K), sprintf("%.2f", ex[[kind]]))
    vb[[sprintf("%s_%d", kind, K)]] <- c(exact = ex[[kind]], sim = mc(kind, K))
  }
}
put("vbar_gmrf_ratio", sprintf("%.1f", vbar_rw(32, ANCH, SIGMA) / vbar_rw(4, ANCH, SIGMA)))
put("vbar_ar1_ratio", sprintf("%.2f", vbar_ar(32, ANCH, TAU, PHI) / vbar_ar(4, ANCH, TAU, PHI)))
put("vbar_check", sprintf("%.3f", max(abs(vapply(vb, function(v) v[["sim"]]/v[["exact"]] - 1, 0)))))

# 5. per cell against per Myr. Scaling the increments by interval width is the obvious
# repair, and it is a real one: it makes the prior a property of the timeline, so refining
# the grid stops changing it. It does not make the walk usable, because what the screen
# bounds is the level and not the slope. Both arms below fix their per-Myr parameter at the
# design width, WREF, so refining the grid genuinely holds the process fixed.
WREF  <- 1.5
S0SQ  <- SIGMA^2 / WREF           # log-rate variance per Myr
T0SQ  <- TSIG^2 / WREF
THETA <- -log(PHI) / WREF         # mean reversion per Myr

mids <- function(d) cumsum(d) - d/2
bm <- function(m, anch, s0sq, d) {                 # Brownian motion on the timeline
  dt <- diff(mids(d))
  cumsum(c(rnorm(1, m, anch), rnorm(length(dt), 0, sqrt(s0sq*dt))))
}
ou <- function(m0, anch, tau, theta, d) {          # stationary OU on the timeline
  tm <- mids(d); m <- rnorm(1, m0, anch)
  x <- numeric(length(tm)); x[1] <- rnorm(1, m, tau)
  for (i in seq_along(tm)[-1]) {
    ph <- exp(-theta*(tm[i]-tm[i-1]))
    x[i] <- m + ph*(x[i-1]-m) + rnorm(1, 0, tau*sqrt(1-ph^2))
  }
  x
}
# The oldest interval reaches the origin, so it is wider than the rest and its lag is the
# one that is not the grid width.
draws_tl <- function(kind, nint, span = 12) {
  w <- span/nint; br <- w*seq_len(nint-1); en <- numeric(N)
  for (k in seq_len(N)) {
    age <- runif(1, AMIN, AMAX)
    d <- c(rep(w, nint-1), age - w*(nint-1))
    if (kind == "bm") { a <- bm(LM,ANCH,S0SQ,d);      tt <- bm(0,TANCH,T0SQ,d) }
    else              { a <- ou(LM,ANCH,TAU,THETA,d); tt <- ou(0,TANCH,TTAU,THETA,d) }
    en[k] <- expected_lineages(rev(exp(a)), rev(exp(a+tt)), age, c(0, age - rev(br)))
  }
  en
}
for (kind in c("bm", "ou")) for (nint in c(4, 8, 16, 32, 64))
  put(sprintf("rej_%s_%d", kind, nint),
      sprintf("%.2f", 100*mean(draws_tl(kind, nint) > 400)))
# and the same for the two per-cell arms at 64, so every row spans the same range
for (kind in c("gmrf", "ar1"))
  put(sprintf("rej_%s_64", kind), sprintf("%.2f", 100*mean(draws(kind, 64)$en > 400)))
for (kind in c("gmrf", "ar1", "bm", "ou")) {
  lo <- as.numeric(g0(sprintf("rej_%s_4", kind)))
  hi <- as.numeric(g0(sprintf("rej_%s_64", kind)))
  put(sprintf("span_%s", kind), sprintf("%.1f", if (hi > lo) hi/lo else lo/hi))
  put(sprintf("dir_%s", kind), if (hi > lo) "up" else "down")
}

# A readable summary first, then the key/value pairs the page substitutes.
g <- function(k) get(k, envir = out)
say <- c(
  sprintf("%d draws per configuration, ceiling 400 lineages, timeline span 12 Myr.", N),
  "",
  sprintf("  %-10s %12s %12s", "intervals", "random walk", "AR(1)"),
  sprintf("  %-10s %11s%% %11s%%", c(4, 8, 16, 32),
          c(g("rej_gmrf_4"), g("rej_gmrf_8"), g("rej_gmrf_16"), g("rej_gmrf_32")),
          c(g("rej_ar1_4"), g("rej_ar1_8"), g("rej_ar1_16"), g("rej_ar1_32"))),
  "",
  sprintf("  expected clade size at 8 intervals, p99.9: random walk %s, AR(1) %s",
          g("q99p9_gmrf"), g("q99p9_ar1")),
  "",
  "  variance of the log-rate path average, the quantity the screen bounds:",
  sprintf("  %-10s %12s %12s", "intervals", "random walk", "AR(1)"),
  sprintf("  %-10s %12s %12s", c(4, 8, 16, 32),
          c(g("vbar_gmrf_4"), g("vbar_gmrf_8"), g("vbar_gmrf_16"), g("vbar_gmrf_32")),
          c(g("vbar_ar1_4"), g("vbar_ar1_8"), g("vbar_ar1_16"), g("vbar_ar1_32"))),
  sprintf("  closed forms, worst relative error against %d draws: %s", N, g("vbar_check")),
  "",
  "  variance per cell against variance per Myr, refused draws:",
  sprintf("  %-10s %12s %12s %12s %12s", "intervals", "walk/cell", "walk/Myr", "AR(1)/cell", "OU/Myr"),
  sprintf("  %-10s %11s%% %11s%% %11s%% %11s%%", c(4, 8, 16, 32, 64),
          vapply(c(4,8,16,32,64), function(k) g(sprintf("rej_gmrf_%d", k)), ""),
          vapply(c(4,8,16,32,64), function(k) g(sprintf("rej_bm_%d", k)), ""),
          vapply(c(4,8,16,32,64), function(k) g(sprintf("rej_ar1_%d", k)), ""),
          vapply(c(4,8,16,32,64), function(k) g(sprintf("rej_ou_%d", k)), "")),
  sprintf("  over that range: walk/cell %s %sx, walk/Myr %s %sx, AR(1)/cell %s %sx, OU/Myr %s %sx",
          g("dir_gmrf"), g("span_gmrf"), g("dir_bm"), g("span_bm"),
          g("dir_ar1"), g("span_ar1"), g("dir_ou"), g("span_ou")))

keys <- sort(ls(out))
pairs <- vapply(keys, function(k) paste0(k, "\t", g(k)), "")
writeLines(c(say, "", "## values", pairs))
