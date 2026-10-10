#!/usr/bin/env Rscript
# Why a random-walk prior on skyline rates cannot be simulated from.
#
# Both priors live on the timeline, not on the grid: the walk's increment variance is per Myr
# and the stationary process's mean reversion is per Myr, so refining the grid at a fixed
# span converges instead of changing the prior.
#
# The span is fixed by the design, so the axis here is how much the rates vary per unit time.
# The two processes are matched on that: the walk's variance per Myr equals the stationary
# process's infinitesimal variance 2*theta*tau^2, which makes them indistinguishable at short
# lags and different only in the restoring force. Expected clade size is exponential in the
# time-weighted average log rate, so what the screen bounds is that average, and the question
# is what each prior spends there to buy a given amount of rate variation.
#
#   rate_prior_tails.R [out.tsv]        key<TAB>value, for the notebook page

set.seed(20261010)
N     <- as.integer(Sys.getenv("NDRAW", "20000"))
WIDTH <- 1.5            # cell width in Myr, held fixed as the span varies
SPAN  <- 12             # the design point
LM    <- -1.5           # log lambda anchor mean
ANCH  <- 0.4            # anchor sd, shared by both priors
SD_DESIGN <- 0.5        # the design point: per-Myr sd of the log rate, both processes
TANCH <- 0.2            # turnover anchor sd
THETA <- 0.5            # stationary: mean reversion per Myr

out <- new.env(); put <- function(k, v) assign(k, v, envir = out)
g0 <- function(k) get(k, envir = out)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "prior_sim.R"))

# One batch under a named prior, at a given per-Myr sd of the log rate. Every arm reads the
# same SD, so the processes are matched on variation per unit time: the walk takes
# Var(x(t)) = a^2 + SD^2 t and the stationary process takes tau^2 = SD^2/(2 theta), which is
# the same infinitesimal variance. Turnover varies as much as lambda, which is how
# rates-prior.Rev ties them, since log_lambda and log_turn share tau_rates.
draws <- function(kind, sd, span = SPAN, width = WIDTH) {
  nint <- round(span/width); br <- width*seq_len(nint-1)
  tau <- sd/sqrt(2*THETA)
  en <- numeric(N); fold <- numeric(N)
  for (k in seq_len(N)) {
    age <- runif(1, 1.25*span, 1.75*span)
    d <- c(rep(width, nint-1), age - width*(nint-1))
    if (kind == "ou") {
      a <- ou(LM, ANCH, tau, THETA, d); tt <- ou(0, TANCH, tau, THETA, d)
      lam <- rev(exp(a)); mu <- rev(exp(a+tt))
    } else if (kind == "bm") {
      a <- bm(LM, ANCH, sd^2, d); tt <- bm(0, TANCH, sd^2, d)
      lam <- rev(exp(a)); mu <- rev(exp(a+tt))
    } else if (kind == "bm_lm") {           # lambda and mu walking independently
      a <- bm(LM, ANCH, sd^2, d); b <- bm(LM, ANCH, sd^2, d)
      lam <- rev(exp(a)); mu <- rev(exp(b))
    } else {                                # horseshoe increments, scaled the same way
      a <- hs(LM, ANCH, sd^2, 0.2, d); tt <- hs(0, TANCH, sd^2, 0.1, d)
      lam <- rev(exp(a)); mu <- rev(exp(a+tt))
    }
    en[k] <- expected_lineages(lam, mu, age, c(0, age - rev(br)))
    fold[k] <- if (all(is.finite(lam)) && all(lam > 0)) max(lam)/min(lam) else NA
  }
  list(en = en, fold = fold)
}

# The common axis: SD is the per-Myr sd of the log rate. The walk takes Var(x(t)) = a^2 +
# SD^2 t. The stationary process takes tau^2 = SD^2 / (2 theta), which gives it the same
# infinitesimal variance, so the two agree over short lags by construction.
SDS <- c(0.08, 0.12, 0.16, 0.20, 0.25, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00)
lab <- function(s) sprintf("%03d", round(100*s))

sweep1 <- function(kind, sd) {
  d <- draws(kind, sd)
  c(fold = median(d$fold, na.rm = TRUE), rej = 100*mean(d$en > 400))
}

# 1. Refusal against the rate variation each prior actually delivers.
tab <- list(bm = NULL, ou = NULL)
for (sd in SDS) for (kind in c("bm", "ou")) {
  v <- sweep1(kind, sd)
  put(sprintf("fold_%s_%s", kind, lab(sd)), sprintf("%.2f", v["fold"]))
  put(sprintf("rej_%s_%s", kind, lab(sd)), sprintf("%.2f", v["rej"]))
  tab[[kind]] <- rbind(tab[[kind]], c(sd = sd, v))
}
# The design has to hold refusal near or under 0.5%, so the usable amount of rate variation
# is where each prior crosses that. Interpolated between the bracketing draws.
BAR <- 0.5
for (kind in c("bm", "ou")) {
  m <- tab[[kind]]; i <- which(m[, "rej"] > BAR)[1]
  put(sprintf("bar_%s", kind),
      if (is.na(i) || i == 1) "NA" else sprintf("%.1f", approx(m[(i-1):i, "rej"], m[(i-1):i, "fold"], xout = BAR)$y))
}
put("bar", sprintf("%.1f", BAR))

# 2. The closed forms, on a uniform grid of K cells spanning T, at matched local variation.
#
#   walk  Var(xbar) - a^2 = SD^2 T (K-1)(2K-1) / 6K^2          -> SD^2 T / 3
#   OU    Var(xbar) - a^2 = (SD^2/2 theta) sum_ij exp(-theta dt|i-j|) / K^2
#
# Their ratio is free of SD: asking for more variation scales both by the same factor. So one
# number says how much more of it the walk spends on the quantity the screen bounds.
xs_bm <- function(s2, K, T) s2*T*(K-1)*(2*K-1)/(6*K^2)
xs_ou <- function(s2, K, T) {
  t <- (seq_len(K) - 0.5) * (T/K)
  (s2/(2*THETA)) * sum(exp(-THETA*abs(outer(t, t, "-")))) / K^2
}
KD <- round(SPAN/WIDTH)
put("xs_ratio", sprintf("%.1f", xs_bm(1, KD, SPAN) / xs_ou(1, KD, SPAN)))
for (K in c(4, 8, 16, 64))
  put(sprintf("xs_ratio_k%d", K), sprintf("%.1f", xs_bm(1, K, SPAN) / xs_ou(1, K, SPAN)))
chk <- c()
for (sd in c(0.2, 0.5)) {
  for (kind in c("bm", "ou")) {
    ex <- if (kind == "bm") xs_bm(sd^2, KD, SPAN) else xs_ou(sd^2, KD, SPAN)
    d <- rep(SPAN/KD, KD)
    mcv <- var(replicate(20000, mean(if (kind == "bm") bm(0, 0, sd^2, d)
                                     else ou(0, 0, sd/sqrt(2*THETA), THETA, d))))
    chk <- c(chk, mcv/ex - 1)
  }
}
put("vbar_check", sprintf("%.3f", max(abs(chk))))
put("anchor_var", sprintf("%.2f", ANCH^2))
put("sd_design", sprintf("%.2f", SD_DESIGN))

# 3. Both priors are processes, not discretisations: refine at a fixed span and they settle.
for (nint in c(4, 8, 32)) for (kind in c("bm", "ou")) {
  d <- draws(kind, SD_DESIGN, SPAN, SPAN/nint)
  put(sprintf("grid_%s_%d", kind, nint), sprintf("%.2f", 100*mean(d$en > 400)))
}

# 4. The tail at the design point, and how far a ceiling has to move to escape it.
for (kind in c("bm", "ou")) {
  d <- draws(kind, SD_DESIGN); e <- d$en[is.finite(d$en)]
  for (q in c(50, 90, 99, 99.9))
    put(sprintf("q%s_%s", sub("[.]", "p", q), kind), signif(quantile(e, q/100), 3))
  put(sprintf("fold_%s", kind), sprintf("%.1f", median(d$fold, na.rm = TRUE)))
  for (cap in c(4e2, 1e4, 1e8))
    put(sprintf("cap_%s_%g", kind, cap), sprintf("%.2f", 100*mean(d$en > cap)))
}

# 5. Two smaller findings. Both have to be read per unit of DELIVERED variation, not at a
# matched nominal scale, since an arm that simply varies the rates less will refuse less.
# So each is also scored against what the Gaussian walk refuses at the same fold.
fold_curve <- vapply(SDS, function(s) as.numeric(g0(sprintf("fold_bm_%s", lab(s)))), 0)
rej_curve  <- vapply(SDS, function(s) as.numeric(g0(sprintf("rej_bm_%s", lab(s)))), 0)
walk_at <- function(f) approx(fold_curve, rej_curve, xout = f, rule = 2)$y
for (kind in c("bm_lm", "hsbm")) {
  d <- draws(kind, SD_DESIGN)
  f <- median(d$fold, na.rm = TRUE); r <- 100*mean(d$en > 400)
  put(sprintf("rej_%s", kind), sprintf("%.1f", r))
  put(sprintf("fold_%s", kind), sprintf("%.1f", f))
  put(sprintf("walkat_%s", kind), sprintf("%.1f", walk_at(f)))
  put(sprintf("cost_%s", kind), sprintf("%.1f", r/walk_at(f)))
}
put("ndraw", N)
put("theta", sprintf("%.1f", THETA))
put("width", sprintf("%.1f", WIDTH))

# A readable summary first, then the key/value pairs the page substitutes.
g <- function(k) as.character(get(k, envir = out))
sd3 <- function(s) sprintf("%03d", round(100*s))
say <- c(
  sprintf("%d draws per configuration, ceiling 400 lineages, span %d Myr in %.1f Myr cells.",
          N, SPAN, WIDTH),
  sprintf("mean reversion %s per Myr, anchor sd %.1f. The walk and the stationary process are",
          g("theta"), ANCH),
  "matched on variation per unit time, so they agree over short lags by construction.",
  "",
  "  log-rate sd   |        random walk        |         stationary",
  "  per Myr       |      fold      refused    |      fold      refused",
  sprintf("  %-13s |  %8s %10s%%    |  %8s %10s%%", sprintf("%.2f", SDS),
          vapply(SDS, function(s) g(sprintf("fold_bm_%s", sd3(s))), ""),
          vapply(SDS, function(s) g(sprintf("rej_bm_%s", sd3(s))), ""),
          vapply(SDS, function(s) g(sprintf("fold_ou_%s", sd3(s))), ""),
          vapply(SDS, function(s) g(sprintf("rej_ou_%s", sd3(s))), "")),
  "",
  sprintf("  rate variation still affordable at a %s%% refusal budget: walk %sx, stationary %sx",
          g("bar"), g("bar_bm"), g("bar_ou")),
  sprintf("  at matched variation the walk puts %sx more variance on the screened average",
          g("xs_ratio")),
  sprintf("  that ratio by cell count: %s",
          paste(sprintf("K=%s %sx", c(4,8,16,64),
                vapply(c(4,8,16,64), function(k) g(sprintf("xs_ratio_k%d", k)), "")), collapse="  ")),
  sprintf("  closed forms, worst relative error against %d draws: %s", N, g("vbar_check")),
  "",
  sprintf("  refining the grid from 4 to 32 cells at a fixed span: walk %s%% to %s%%, stationary %s%% to %s%%",
          g("grid_bm_4"), g("grid_bm_32"), g("grid_ou_4"), g("grid_ou_32")),
  "",
  sprintf("  expected clade size at the design point, p99.9: walk %s, stationary %s",
          g("q99p9_bm"), g("q99p9_ou")))

keys <- sort(ls(out))
pairs <- vapply(keys, function(k) paste0(k, "\t", g(k)), "")
writeLines(c(say, "", "## values", pairs))
