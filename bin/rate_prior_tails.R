#!/usr/bin/env Rscript
# Why a random-walk prior on skyline rates cannot be simulated from.
#
# Both priors here live on the timeline, not on the grid: the walk's increment variance is
# per Myr and the stationary process's mean reversion is per Myr, so refining the grid at a
# fixed span converges instead of changing the prior. That makes the comparison a comparison
# of two processes. Expected clade size is exponential in the time-weighted average log rate,
# so what the simulator's screen bounds is that average. A diffusion holds no prior
# information about it and a stationary process does, and the gap widens with the length of
# the record rather than with how finely it is cut.
#
#   rate_prior_tails.R [out.tsv]        key<TAB>value, for the notebook page

set.seed(20261010)
N     <- as.integer(Sys.getenv("NDRAW", "20000"))
WIDTH <- 1.5            # cell width in Myr, held fixed as the span varies
SPAN  <- 12             # the design point
LM    <- -1.5           # log lambda anchor mean
ANCH  <- 0.4            # anchor sd, shared by both priors
S0SQ  <- 0.4^2 / WIDTH  # walk: log-rate variance per Myr
T0SQ  <- 0.2^2 / WIDTH  # the same, for log turnover
TAU   <- 0.5            # stationary sd of log lambda
TTAU  <- 0.3            # stationary sd of log turnover
TANCH <- 0.2            # turnover anchor sd
THETA <- 0.5            # stationary: mean reversion per Myr

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

# Cell midpoints, youngest first. d is the cell widths; the oldest cell reaches the origin
# and is wider than the rest, because the breakpoints cannot scale with an estimated origin.
mids <- function(d) cumsum(d) - d/2

# Brownian motion sampled at the midpoints: increment variance proportional to elapsed time.
bm <- function(m, anch, s0sq, d) {
  dt <- diff(mids(d))
  cumsum(c(rnorm(1, m, anch), rnorm(length(dt), 0, sqrt(s0sq*dt))))
}
# The same with horseshoe increments, after Magee et al. 2020, scaled the same way.
hs <- function(m, anch, s0sq, zeta, d) {
  dt <- diff(mids(d)); n <- length(dt)
  cumsum(c(rnorm(1, m, anch),
           rnorm(n, 0, sqrt(s0sq*dt) * abs(rcauchy(1,0,zeta)) * abs(rcauchy(n,0,1)))))
}
# Stationary OU sampled at the midpoints: correlation exp(-theta * lag in Myr).
ou <- function(m0, anch, tau, theta, d) {
  t <- mids(d); m <- rnorm(1, m0, anch)
  x <- numeric(length(t)); x[1] <- rnorm(1, m, tau)
  for (i in seq_along(t)[-1]) {
    ph <- exp(-theta*(t[i]-t[i-1]))
    x[i] <- m + ph*(x[i-1]-m) + rnorm(1, 0, tau*sqrt(1-ph^2))
  }
  x
}

# One batch under a named prior. The span sets the grid and the age prior together, so a
# longer record is more record and not a differently shaped one.
draws <- function(kind, span = SPAN, width = WIDTH) {
  nint <- round(span/width); br <- width*seq_len(nint-1)
  en <- numeric(N); fold <- numeric(N)
  for (k in seq_len(N)) {
    age <- runif(1, 1.25*span, 1.75*span)
    d <- c(rep(width, nint-1), age - width*(nint-1))
    if (kind == "ou") {
      a <- ou(LM, ANCH, TAU, THETA, d); tt <- ou(0, TANCH, TTAU, THETA, d)
      lam <- rev(exp(a)); mu <- rev(exp(a+tt))
    } else if (kind == "bm") {
      a <- bm(LM, ANCH, S0SQ, d); tt <- bm(0, TANCH, T0SQ, d)
      lam <- rev(exp(a)); mu <- rev(exp(a+tt))
    } else if (kind == "bm_lm") {           # lambda and mu walking independently
      a <- bm(LM, ANCH, S0SQ, d); b <- bm(LM, ANCH, S0SQ, d)
      lam <- rev(exp(a)); mu <- rev(exp(b))
    } else {                                # horseshoe increments
      a <- hs(LM, ANCH, S0SQ, 0.2, d); tt <- hs(0, TANCH, T0SQ, 0.1, d)
      lam <- rev(exp(a)); mu <- rev(exp(a+tt))
    }
    en[k] <- expected_lineages(lam, mu, age, c(0, age - rev(br)))
    fold[k] <- if (all(is.finite(lam)) && all(lam > 0)) max(lam)/min(lam) else NA
  }
  list(en = en, fold = fold)
}

SPANS <- c(6, 12, 24, 48, 96)

# 1. The claim: refusal against the length of the record, at fixed resolution.
for (sp in SPANS) for (kind in c("bm", "ou")) {
  d <- draws(kind, sp)
  put(sprintf("rej_%s_T%d", kind, sp), sprintf("%.2f", 100*mean(d$en > 400)))
}
# The ratio between the priors AT a span is the comparison. The trend within one prior is
# not, because the screen bounds T * xbar and so carries a T^2 that is nothing to do with
# which prior was chosen: refusal rises with the span under both. What belongs to the prior
# is Var(xbar) above the anchor, which is what the excess columns isolate.
for (sp in SPANS) {
  r <- as.numeric(c(g0(sprintf("rej_bm_T%d", sp)), g0(sprintf("rej_ou_T%d", sp))))
  put(sprintf("ratio_T%d", sp), sprintf("%.0f", r[1]/r[2]))
}

# 2. The closed forms, on a uniform grid of K cells spanning T.
#
#   walk  Var(xbar) = a^2 + s0^2 T (K-1)(2K-1) / 6K^2      -> a^2 + s0^2 T / 3
#   OU    Var(xbar) = a^2 + (tau^2/K^2) sum_ij exp(-theta dt |i-j|)
#                                                          -> a^2 + O(1/T)
#
# The walk's grows without bound in the span. The OU's falls to the anchor for any span.
# Both are checked against draws below.
vbar_bm <- function(K, T) ANCH^2 + S0SQ*T*(K-1)*(2*K-1)/(6*K^2)
vbar_ou <- function(K, T) {
  t <- (seq_len(K) - 0.5) * (T/K)
  ANCH^2 + TAU^2 * sum(exp(-THETA*abs(outer(t, t, "-")))) / K^2
}
vbar_bm_lim <- function(T) ANCH^2 + S0SQ*T/3
vbar_ou_lim <- function(T) ANCH^2 + (2*TAU^2/(THETA*T))*(1 - (1-exp(-THETA*T))/(THETA*T))

mc_vbar <- function(kind, T, n = 20000) {
  K <- round(T/WIDTH); d <- rep(T/K, K)
  var(replicate(n, mean(if (kind == "bm") bm(LM, ANCH, S0SQ, d)
                        else ou(LM, ANCH, TAU, THETA, d))))
}
chk <- c()
for (sp in SPANS) {
  K <- round(sp/WIDTH)
  for (kind in c("bm", "ou")) {
    ex <- if (kind == "bm") vbar_bm(K, sp) else vbar_ou(K, sp)
    put(sprintf("vbar_%s_T%d", kind, sp), sprintf("%.2f", ex))
    chk <- c(chk, mc_vbar(kind, sp)/ex - 1)
  }
}
put("vbar_check", sprintf("%.3f", max(abs(chk))))
for (sp in SPANS) {
  K <- round(sp/WIDTH)
  xb <- vbar_bm(K, sp) - ANCH^2; xo <- vbar_ou(K, sp) - ANCH^2
  put(sprintf("xs_bm_T%d", sp), sprintf("%.2f", xb))
  put(sprintf("xs_ou_T%d", sp), sprintf("%.2f", xo))
  put(sprintf("xs_ratio_T%d", sp), sprintf("%.0f", xb/xo))
}
put("vbar_bm_ratio", sprintf("%.1f", vbar_bm(round(max(SPANS)/WIDTH), max(SPANS)) /
                                     vbar_bm(round(min(SPANS)/WIDTH), min(SPANS))))
put("vbar_ou_ratio", sprintf("%.2f", vbar_ou(round(max(SPANS)/WIDTH), max(SPANS)) /
                                     vbar_ou(round(min(SPANS)/WIDTH), min(SPANS))))
put("vbar_bm_lim", sprintf("%.2f", vbar_bm_lim(SPAN)))
put("vbar_ou_lim", sprintf("%.2f", vbar_ou_lim(SPAN)))
put("anchor_var", sprintf("%.2f", ANCH^2))

# 3. Both priors are processes, not discretisations: refine at a fixed span and they settle.
for (nint in c(4, 8, 32)) for (kind in c("bm", "ou")) {
  d <- draws(kind, SPAN, SPAN/nint)
  put(sprintf("grid_%s_%d", kind, nint), sprintf("%.2f", 100*mean(d$en > 400)))
}

# 4. The tail at the design point, and how far a ceiling has to move to escape it.
for (kind in c("bm", "ou")) {
  d <- draws(kind, SPAN); e <- d$en[is.finite(d$en)]
  for (q in c(50, 90, 99, 99.9))
    put(sprintf("q%s_%s", sub("[.]", "p", q), kind), signif(quantile(e, q/100), 3))
  put(sprintf("fold_%s", kind), sprintf("%.1f", median(d$fold, na.rm = TRUE)))
  for (cap in c(4e2, 1e4, 1e8))
    put(sprintf("cap_%s_%g", kind, cap), sprintf("%.2f", 100*mean(d$en > cap)))
}

# 5. Two smaller findings, both at the design point.
for (kind in c("bm_lm", "hsbm")) {
  d <- draws(kind, SPAN)
  put(sprintf("rej_%s", kind), sprintf("%.2f", 100*mean(d$en > 400)))
  put(sprintf("fold_%s", kind), sprintf("%.1f", median(d$fold, na.rm = TRUE)))
}
put("ndraw", N)
put("theta", sprintf("%.1f", THETA))
put("width", sprintf("%.1f", WIDTH))

# A readable summary first, then the key/value pairs the page substitutes.
g <- function(k) as.character(get(k, envir = out))
say <- c(
  sprintf("%d draws per configuration, ceiling 400 lineages, cells %s Myr wide.", N, g("width")),
  sprintf("walk variance %.4f per Myr, stationary mean reversion %s per Myr, anchor sd %.1f.",
          S0SQ, g("theta"), ANCH),
  "",
  "  refused draws, and the variance of the log-rate time average, against record length:",
  sprintf("  %-10s %12s %12s %12s %12s", "span/Myr", "walk", "stationary", "Var walk", "Var stat"),
  sprintf("  %-10s %11s%% %11s%% %12s %12s", SPANS,
          vapply(SPANS, function(s) g(sprintf("rej_bm_T%d", s)), ""),
          vapply(SPANS, function(s) g(sprintf("rej_ou_T%d", s)), ""),
          vapply(SPANS, function(s) g(sprintf("vbar_bm_T%d", s)), ""),
          vapply(SPANS, function(s) g(sprintf("vbar_ou_T%d", s)), "")),
  sprintf("  walk-to-stationary refusal ratio by span: %s",
          paste(sprintf("%dMyr %sx", SPANS,
                vapply(SPANS, function(s) g(sprintf("ratio_T%d", s)), "")), collapse="  ")),
  sprintf("  the prior's own share, Var(xbar) above the anchor sd^2 of %s:", g("anchor_var")),
  sprintf("  %-10s %12s %12s %12s", "span/Myr", "walk", "stationary", "ratio"),
  sprintf("  %-10s %12s %12s %11sx", SPANS,
          vapply(SPANS, function(s) g(sprintf("xs_bm_T%d", s)), ""),
          vapply(SPANS, function(s) g(sprintf("xs_ou_T%d", s)), ""),
          vapply(SPANS, function(s) g(sprintf("xs_ratio_T%d", s)), "")),
  sprintf("  closed forms, worst relative error against %d draws: %s", N, g("vbar_check")),
  "",
  sprintf("  at a %d Myr span, refining the grid from 4 to 32 cells: walk %s%% to %s%%, stationary %s%% to %s%%",
          SPAN, g("grid_bm_4"), g("grid_bm_32"), g("grid_ou_4"), g("grid_ou_32")),
  "",
  sprintf("  expected clade size at the design point, p99.9: walk %s, stationary %s",
          g("q99p9_bm"), g("q99p9_ou")))

keys <- sort(ls(out))
pairs <- vapply(keys, function(k) paste0(k, "\t", g(k)), "")
writeLines(c(say, "", "## values", pairs))
