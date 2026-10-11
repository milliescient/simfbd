#!/usr/bin/env Rscript
# Generators shared by the prior analyses, so two pages cannot disagree about what a walk or
# an OU process is. Sourced, not run. The caller supplies N, LM, ANCH and the rest.
#
# All three fields are defined on the timeline: the walk's increment variance is per Myr and
# the OU's mean reversion is per Myr, so refining the grid converges instead of changing the
# prior. A field is sampled at cell MIDPOINTS, and the oldest cell reaches the origin and is
# wider than the others, so its lag is not the grid width.

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

