######################################################
##                    Chapter 3                     ##
##      Treeless skyline FBDR simulation study      ##
##                Simulation script                 ##
##            Bruno do Rosario Petrucci             ##
######################################################

###
# load packages

# paleobuddy 
library(paleobuddy)

# The fork's closed-form draw carries CRAN's version string, so a reinstall replaces it
# silently and shows up only as a run that is orders of magnitude slower.
if (!exists("rexp.var") || !any(grepl("stepfun", deparse(body(rexp.var))))) {
  warning("paleobuddy has no stepfun fast path, so this will be far slower than it should be",
          call. = FALSE, immediate. = TRUE)
}

# ape
library(ape)




###
# write minor auxiliary functions

# smarter dir.create function
smart_dir_create <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

# equivalent to colSums for max
colMaxes <- function(df) {
  # make return
  res <- c()
  
  # iterate through columns
  for (i in 1:ncol(df)) {
    # add max to res
    res <- c(res, max(df[, i]))
  }
  
  # name res
  names(res) <- colnames(df)
  
  # return res
  return(res)
}

# same for min
colMins <- function(df) {
  # make return
  res <- c()
  
  # iterate through columns
  for (i in 1:ncol(df)) {
    # add max to res
    res <- c(res, min(df[, i]))
  }
  
  # name res
  names(res) <- colnames(df)
  
  # return res
  return(res)
}

# making time bins
make_bins <- function(age, unc) {
  # find range based on unc
  range <- switch(unc,
                  "mid" = {c(age / 6, age / 5)},
                  "low" = {c(age / 10, age / 8)},
                  "high" = {c(age / 4, age / 3)})
  
  # create bins vector
  bins <- c(0)
  
  # while bins doesn't have age
  while (!(age %in% bins)) {
    # if the highest bin is less than the max range, just add age
    if (max(bins) + range[2] > age) {
      bins <- c(bins, age)
    } else {
      # if not, draw a new bin
      bin <- runif(1, range[1], range[2])
      
      # if max(bins) + bin is higher than age, add age instead
      if (max(bins) + bin > age) {
        bins <- c(bins, age)
      } else {
        # add to bins
        bins <- c(bins, max(bins) + bin)
        
      }
    }
  }
  
  # return bins
  return(bins)
}

###
# SBC simulation constants
#
# Every generative parameter comes from the batch config.sh, exported into the
# environment by bin/run.sh. A missing key is fatal: a silent default is how the
# origin prior drifted from the value the analysis assumed.
cfg <- function(k) {
  v <- Sys.getenv(k)
  if (v == "") stop("config key ", k, " is not set; bin/run.sh must export it")
  as.numeric(v)
}

# a key a batch may leave unset, so older configs keep working
cfg_opt <- function(k, default) {
  v <- Sys.getenv(k)
  if (v == "") default else as.numeric(v)
}

NINTERVALS    <- as.integer(cfg("NINTERVALS"))
INTERVAL_WIDTH <- cfg("INTERVAL_WIDTH")
LMEAN <- cfg("LMEAN"); LSD <- cfg("LSD")
MMEAN <- cfg("MMEAN"); MSD <- cfg("MSD")
PMEAN <- cfg("PMEAN"); PSD <- cfg("PSD")
AGE_MIN <- cfg("AGE_MIN"); AGE_MAX <- cfg("AGE_MAX")
BIN_WIDTH <- cfg("BIN_WIDTH"); BIN_MAX <- cfg("BIN_MAX")

# Off the grid, the top occurrence bin runs past AGE_MAX, so no origin the prior can draw
# covers it and tau_1 is squeezed into what is left.
if (abs(AGE_MAX %% BIN_WIDTH) > 1e-9) {
  warning(sprintf("AGE_MAX %g is off the BIN_WIDTH %g grid, so the top bin extends past it",
                  AGE_MAX, BIN_WIDTH), call. = FALSE, immediate. = TRUE)
}

# Ceiling on lineages ever born. Unset means no ceiling, which is right when the origin
# is shallow. Deeper origins need one: the oldest interval absorbs the extra depth, so
# diversity is exponential in it and a draw in the upper tail of lambda never terminates.
MAX_LINEAGES <- cfg_opt("MAX_LINEAGES", Inf)
SIM_TIMEOUT <- cfg_opt("SIM_TIMEOUT", Inf)   # seconds one bd.sim call may take
RHO <- cfg_opt("RHO", 1)                     # chance an extant lineage is seen at the present
# The record stops here rather than at zero. Occurrences younger than PRESENT are discarded
# and a lineage still alive there is a boundary survivor, which rho then reports or misses.
PRESENT <- cfg_opt("PRESENT", 0)
# Chance a sampled interior occurrence is reported. 1 emits no thinned record, since it
# would duplicate the complete one. The extremes are always reported, which is what makes
# the rate estimable rather than confounded with psi.
REPORT_PROB <- cfg_opt("REPORT_PROB", 1)
if (PRESENT > 0 && abs(PRESENT %% BIN_WIDTH) > 1e-9) {
  warning(sprintf("PRESENT %g is off the BIN_WIDTH %g grid, so a bin straddles it",
                  PRESENT, BIN_WIDTH), call. = FALSE, immediate. = TRUE)
}
# Fossil-sampled taxa a replicate must hold to be kept, both bounds inclusive. Rejection is on
# a function of the record and redraws theta, so it leaves the posterior alone. At rho < 1 the
# extant singletons are added after this test, so it counts the fossil record only.
MIN_TAXA <- as.integer(cfg_opt("MIN_TAXA", 1))
# Order-1 GMRF on the log rates: anchor at the present, random walk back through time. Unset
# keeps the independent draw, so an older batch reproduces.
GMRF_SD <- cfg_opt("GMRF_SD", 0)
# Hierarchical variant: the walk's own standard deviation is drawn per replicate, so one batch
# spans clades that vary a lot through time and clades that barely vary, and the analysis
# estimates it. Setting GMRF_SD_MEAN switches this on, and GMRF_SD is then ignored.
#
# The hyperprior is lognormal on purpose. A heavy-tailed one puts real mass on walks steep
# enough to trip the MAX_LINEAGES screen, and that screen is a function of theta alone, so it
# truncates the prior the simulator draws from while the analysis keeps the untruncated one.
# The ranks would be non-uniform for that reason alone.
#
# lambda and mu carry the through-time signal; psi gets its own, usually narrower, so that a
# replicate's sampling regime stays a meaningful label to stratify on.
# Stationary alternative to the random walk: an OU process on log lambda, log turnover and
# log psi, sampled at the interval midpoints, with mean-reversion half-life OU_HALFLIFE in
# Myr and stationary sds LSD, TSD, PSD. Setting OU_HALFLIFE switches it on and GMRF_SD is
# then ignored.
#
# Two reasons to prefer it. A random walk's marginal prior widens with distance from the
# anchor, so the oldest interval gets a far vaguer prior than the youngest; an AR(1) is
# stationary, so every interval gets the same one. And parameterising on turnover rather
# than mu keeps net diversification tied to lambda, which is what stops a drifting
# lambda - mu compounding into clades that trip the MAX_LINEAGES screen. Measured, that
# screen fires on 18% of random-walk draws and effectively none of these.
#
# OU_HALFLIFE is in Myr, so the prior is a property of the timeline and not of the grid.
OU_HALFLIFE <- cfg_opt("OU_HALFLIFE", NA)
TMEAN  <- cfg_opt("TMEAN", 0)
TSD    <- cfg_opt("TSD", 0.3)
OU <- !is.na(OU_HALFLIFE)
OU_THETA <- if (OU) log(2) / OU_HALFLIFE else NA
# The stationary prior's own hierarchy, and the two things it has to do at once. The MEAN of
# each path is drawn per replicate, which sweeps regimes ACROSS replicates, and the
# stationary sd is drawn per replicate, which varies rates WITHIN one. Both are parameters
# the analysis estimates, so both are truth and both have to rank uniformly.
#
# The mean reversion is fixed, and is given as a half-life in Myr rather than as a lag-1
# correlation. A correlation per cell would make the prior a property of the grid: refining
# the grid at a fixed span would change the process. A half-life does not, so the interval
# count is a resolution choice and nothing else. Estimating it would be estimating a
# property of the discretisation, so it stays fixed.
#
# The intercept sds at 0 keep the fixed mean an older batch drew with.
LANCH <- cfg_opt("LANCH", 0)
TANCH <- cfg_opt("TANCH", 0)
PANCH <- cfg_opt("PANCH", 0)
OU_SD_MEAN  <- cfg_opt("OU_SD_MEAN", NA)
OU_SD_SD    <- cfg_opt("OU_SD_SD", 0)
OU_PSI_MEAN <- cfg_opt("OU_PSI_MEAN", NA)
OU_PSI_SD   <- cfg_opt("OU_PSI_SD", 0)
OU_HIER <- !is.na(OU_SD_MEAN)
if (OU_HIER && is.na(OU_PSI_MEAN)) stop("OU_SD_MEAN set without OU_PSI_MEAN")
if (OU_HIER && !OU) stop("OU_SD_MEAN set without OU_HALFLIFE")
GMRF_SD_MEAN  <- cfg_opt("GMRF_SD_MEAN", NA)
GMRF_SD_SD    <- cfg_opt("GMRF_SD_SD", 0)
GMRF_PSI_MEAN <- cfg_opt("GMRF_PSI_MEAN", NA)
GMRF_PSI_SD   <- cfg_opt("GMRF_PSI_SD", 0)
GMRF_HIER <- !is.na(GMRF_SD_MEAN)
if (GMRF_HIER && is.na(GMRF_PSI_MEAN)) stop("GMRF_SD_MEAN set without GMRF_PSI_MEAN")
MAX_TAXA <- cfg_opt("MAX_TAXA", Inf)
# "record" makes MIN_TAXA and MAX_TAXA count every reported taxon, the survivors seen at the
# present included; unset keeps the fossil-only count, so an older batch reproduces.
TAXA_COUNT <- Sys.getenv("TAXA_COUNT", "fossil")
if (!TAXA_COUNT %in% c("fossil", "record")) stop("TAXA_COUNT must be fossil or record")
# a floor on fossil-sampled taxa alone, for analyses that leave the present-only taxa out
MIN_FOSSIL_TAXA <- as.integer(cfg_opt("MIN_FOSSIL_TAXA", 0))

# Rate breakpoints (before the present) and fossil bins. Both are fixed rather than
# derived from age: the analysis estimates the origin, so anything it reads that
# scales with age hands it the answer. paleobuddy counts shift times forward from
# the start, so these are age - cov_breaks there.
# the intervals partition the process's own span, which starts at PRESENT rather than zero
cov_breaks <- PRESENT + INTERVAL_WIDTH * seq_len(NINTERVALS - 1)
cov_bins <- seq(0, BIN_MAX, BIN_WIDTH)

# The generative draw, in one place. Both the initial draw and the rejection redraw
# call this, so they cannot diverge.
# Returned oldest-first, the order paleobuddy's forward shifts need. The walk itself runs
# from the present backwards, so index 1 of the walk is the youngest interval.
# Lags between interval midpoints, in Myr, youngest first. The oldest interval runs from the
# last breakpoint to the origin, so its width is taken at the age prior's mean: making it the
# replicate's own age would put the origin in the rate prior, and the breakpoints are kept
# clear of age for that reason. Measured, the two differ by 0.02 points of refusal rate.
OU_WIDTHS <- c(rep(INTERVAL_WIDTH, NINTERVALS - 1),
               (AGE_MIN + AGE_MAX) / 2 - PRESENT - INTERVAL_WIDTH * (NINTERVALS - 1))
OU_LAGS <- diff(cumsum(OU_WIDTHS) - OU_WIDTHS / 2)
if (OU) {
  if (OU_WIDTHS[NINTERVALS] <= 0)
    stop("the grid spans past the mean age, so the oldest interval has no width")
  OU_PHI_VEC <- exp(-OU_THETA * OU_LAGS)
}

# One OU path on the log scale, youngest interval first, sampled on an irregular grid.
ou_path <- function(m, tau, n) {
  x <- numeric(n); x[1] <- rnorm(1, m, tau)
  if (n > 1) for (i in 2:n) {
    phi <- OU_PHI_VEC[i - 1]
    x[i] <- m + phi * (x[i - 1] - m) + rnorm(1, 0, tau * sqrt(1 - phi^2))
  }
  x
}

draw_rates <- function(m, s, n, walk_sd = GMRF_SD) {
  if (walk_sd <= 0) return(rlnorm(n, m, s))
  rev(exp(cumsum(c(rnorm(1, m, s), rnorm(n - 1, 0, walk_sd)))))
}

# Truth, named, in one place. It used to be assembled by positional appends, which put
# report_prob between two hyperparameters when both were on, and read the hyperparameters out
# of the global environment, so a replicate accepted without a rejection recorded the
# previous replicate's values. Names travel with the values now, so neither can recur.
truth_of <- function(th, lambda, mu, psi, age) {
  K <- length(lambda)
  v <- c(setNames(rev(lambda), paste0("lambda", seq_len(K))),
         setNames(rev(mu),     paste0("mu", seq_len(K))),
         setNames(rev(psi),    paste0("psi", seq_len(K))))
  if (th$lambda_a > 0) v <- c(v, lambda_a = th$lambda_a)
  if (REPORT_PROB < 1) v <- c(v, report_prob = th$report_prob)
  v <- c(v, age = age)
  if (OU) {
    if (LANCH > 0) v <- c(v, m_lambda = th$m_lambda)
    if (TANCH > 0) v <- c(v, m_turnover = th$m_turnover)
    if (PANCH > 0) v <- c(v, m_psi = th$m_psi)
    if (OU_HIER)   v <- c(v, tau_rates = th$tau_rates, tau_psi = th$tau_psi)
  } else if (GMRF_HIER) {
    v <- c(v, gmrf_sd_rates = th$sd_rates, gmrf_sd_psi = th$sd_psi)
  }
  v
}

draw_theta <- function(lambda_a) {
  if (OU) {
    # lambda and turnover share one stationary sd, psi gets its own, so that a replicate's
    # sampling regime stays a label worth stratifying on
    tau_rates <- if (OU_HIER) rlnorm(1, OU_SD_MEAN, OU_SD_SD) else LSD
    tau_psi   <- if (OU_HIER) rlnorm(1, OU_PSI_MEAN, OU_PSI_SD) else PSD
    m_lambda   <- if (LANCH > 0) rnorm(1, LMEAN, LANCH) else LMEAN
    m_turnover <- if (TANCH > 0) rnorm(1, TMEAN, TANCH) else TMEAN
    m_psi      <- if (PANCH > 0) rnorm(1, PMEAN, PANCH) else PMEAN
    # paleobuddy wants oldest-first, so reverse once, after lambda and mu are assembled
    ll <- ou_path(m_lambda, tau_rates, NINTERVALS)
    lt <- ou_path(m_turnover, tau_rates, NINTERVALS)
    lp <- ou_path(m_psi, tau_psi, NINTERVALS)
    return(list(age = runif(1, AGE_MIN, AGE_MAX),
                lambda = rev(exp(ll)), mu = rev(exp(ll + lt)), psi = rev(exp(lp)),
                m_lambda = m_lambda, m_turnover = m_turnover, m_psi = m_psi,
                tau_rates = tau_rates, tau_psi = tau_psi,
                lambda_a = if (lambda_a > 0) rlnorm(1, m_lambda, tau_rates) else lambda_a,
                report_prob = if (REPORT_PROB < 1) runif(1) else 1))
  }
  sd_rates <- if (GMRF_HIER) rlnorm(1, GMRF_SD_MEAN, GMRF_SD_SD) else GMRF_SD
  sd_psi   <- if (GMRF_HIER) rlnorm(1, GMRF_PSI_MEAN, GMRF_PSI_SD) else GMRF_SD
  list(age    = runif(1, AGE_MIN, AGE_MAX),
       lambda = draw_rates(LMEAN, LSD, NINTERVALS, sd_rates),
       mu     = draw_rates(MMEAN, MSD, NINTERVALS, sd_rates),
       psi    = draw_rates(PMEAN, PSD, NINTERVALS, sd_psi),
       sd_rates = sd_rates, sd_psi = sd_psi,
       lambda_a = if (lambda_a > 0) rlnorm(1, LMEAN, LSD) else lambda_a,
       report_prob = if (REPORT_PROB < 1) runif(1) else 1)
}

# Expected lineages ever born, 1 + integral of lambda(t)E[N(t)]: E[N] grows as exp(lambda-mu)
# within an interval, so both integrals are closed form. Screening on this before bd.sim is
# what keeps a runaway draw from hanging, and it is a function of theta alone, so a redraw
# leaves the accepted set a truncation of the prior rather than a conditioning on the data.
expected_lineages <- function(lambda, mu, age, shifts) {
  edges <- c(shifts, age)
  total <- 1
  n <- 1
  for (j in seq_along(lambda)) {
    dt <- edges[j + 1] - edges[j]
    r <- lambda[j] - mu[j]
    total <- total + lambda[j] * n * (if (abs(r) < 1e-12) dt else expm1(r * dt) / r)
    n <- n * exp(r * dt)
    if (!is.finite(total)) return(Inf)
  }
  total
}

# Expected size screens out the astronomical draws, but with rate shifts paleobuddy takes
# its general path, whose cost is not a function of tree size alone, so a rare draw still
# runs for minutes. NULL means it outran the clock and the caller should redraw. Which
# draws hit this depends on machine load, so a batch reproduces from its seeds only up to
# the replicates that time out; sim_timeouts reports how many that was.
sim_timeouts <- 0
# draws refused by the expected-clade screen. This one is a function of theta alone, so it
# truncates the prior and the ranks are non-uniform by exactly that much: the number belongs
# in the record of every batch, not just in a pilot.
sim_screened <- 0
sim_drawn <- 0
bd_sim_bounded <- function(n0, lambda, mu, tMax, lShifts, mShifts, nFinal) {
  call_it <- function() bd.sim(n0, lambda, mu, tMax, lShifts = lShifts,
                               mShifts = mShifts, nFinal = c(nFinal, Inf))
  if (!is.finite(SIM_TIMEOUT)) return(call_it())
  setTimeLimit(elapsed = SIM_TIMEOUT, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf))
  tryCatch(call_it(), error = function(e) {
    if (!grepl("elapsed time limit", conditionMessage(e))) stop(e)
    sim_timeouts <<- sim_timeouts + 1
    NULL
  })
}

# Reporting models, by the sim's model number. An SBC rep emits both from one
# simulated record, so they share a tree, a fossil set and a timeline.
sbc_models <- c(complete = 1, firstlast = 2)
if (REPORT_PROB < 1) sbc_models <- c(sbc_models, thinned = 3)

# Which of a species' occurrences get reported, under each model.
retain_occs <- function(occs, model, rp) {
  # complete: the whole record
  if (model == 1) return(occs)

  # first/last (model 2): the oldest and youngest. Keep both even when they fall
  # in the same bin, so a bracketed pair is reported as two (possibly identical)
  # occurrences and the count >= 2 likelihood marginalizes the interior within the
  # bin. Only a genuine single fossil yields count == 1.
  if (model == 2) {
    if (nrow(occs) > 1) return(occs[c(1, nrow(occs)), ])
    return(occs[1, ])
  }
  # thinned (model 3): both extremes always, each interior kept with probability rp
  if (nrow(occs) <= 2) return(occs)
  interior <- occs[-c(1, nrow(occs)), , drop = FALSE]
  keep <- runif(nrow(interior)) < rp
  rbind(occs[1, ], interior[keep, , drop = FALSE], occs[nrow(occs), ])
}

# Anagenetic speciation relabels the same tree rather than reshaping it: split each
# species at Poisson(lambda_a) times in its lifespan and reassign occurrences by time.
add_anagenesis <- function(sim, fossils, lambda_a) {
  n0 <- length(sim$TS)

  # segment boundaries per species, oldest first
  starts <- ends <- vector("list", n0)
  for (k in 1:n0) {
    te <- ifelse(is.na(sim$TE[k]), 0, sim$TE[k])
    ev <- sort(runif(rpois(1, lambda_a * (sim$TS[k] - te)), te, sim$TS[k]),
               decreasing = TRUE)
    starts[[k]] <- c(sim$TS[k], ev)
    ends[[k]] <- c(ev, te)
  }

  # new indices run in (species, segment) order
  n_seg <- sapply(starts, length)
  offset <- cumsum(c(0, n_seg))[1:n0]

  # segment of species k alive at time t
  seg_at <- function(k, t) {
    s <- which(t <= starts[[k]] & t >= ends[[k]])
    if (length(s) == 0) n_seg[k] else s[1]
  }

  N <- sum(n_seg)
  TS <- TE <- rep(NA_real_, N)
  PAR <- rep(NA_integer_, N)
  EXTANT <- rep(FALSE, N)

  for (k in 1:n0) {
    for (s in 1:n_seg[k]) {
      i <- offset[k] + s
      last <- (s == n_seg[k]) && sim$EXTANT[k]
      TS[i] <- starts[[k]][s]
      TE[i] <- ifelse(last, NA_real_, ends[[k]][s])
      EXTANT[i] <- last
      # a later segment descends from the previous one; the first keeps the original
      # parent, joined at whichever of that parent's segments was then alive
      PAR[i] <- if (s > 1) offset[k] + s - 1 else
        if (is.na(sim$PAR[k])) NA_integer_ else
          offset[sim$PAR[k]] + seg_at(sim$PAR[k], sim$TS[k])
    }
  }

  if (nrow(fossils) > 0) {
    k <- as.integer(sub("^t", "", fossils$Species))
    i <- offset[k] + mapply(seg_at, k, fossils$SampT)
    fossils$Species <- paste0("t", i)
    fossils$Extant <- EXTANT[i]
  }

  sim <- list(TE = TE, TS = TS, PAR = PAR, EXTANT = EXTANT)
  class(sim) <- "sim"

  list(SIM = sim, FOSSILS = fossils)
}

###
# write auxiliary simulation functions

# simulate one rep. lambda_a > 0 adds anagenetic speciation; it is redrawn with the
# other rates so the rejection loop keeps the joint draw intact.
simulate_rep <- function(rates, age, shifts, report_prob = 1,
                         model, unc, extant_singletons,
                         sbc = FALSE, lambda_a = 0, origin_sampled = FALSE,
                         theta = NULL) {
  th <- theta
  ## set parameters
  # number of initial species
  n0 <- 1
  
  # total simulation time
  tMax <- age
  
  # lambda
  lambda <- rates[[1]]
  
  # mu
  mu <- rates[[2]]
  
  # psi
  psi <- rates[[3]]
  
  # shifts
  lShifts <- shifts[[1]]
  mShifts <- shifts[[2]]
  pShifts <- shifts[[3]]
  
  # nFinal based on sbc
  # SBC fix: do NOT impose a diversification floor on the tree
  # (the old value 5 conditions on >=5 total species, which biases the rates).
  nFinal <- ifelse(sbc, 1, 10)
  
  # get bins
  if (sbc) {
    # fixed grid; binning up to age would put a bin edge on the origin
    bins <- cov_bins
  } else {
    # if not, make bins
    bins <- make_bins(age, unc)
  }
  
  ##
  # run simulations
  
  # create conditions boolean
  cond <- FALSE

  # Every rejection path redraws theta jointly, so none of them can diverge. Deriving every
  # local from the one draw in one place is what keeps a redraw from leaving part of the old
  # one behind: report_prob used to survive a rejection, which made the retained pair a
  # partial redraw rather than a draw from the joint.
  use_theta <- function(t) {
    th <<- t
    age <<- t$age; tMax <<- t$age
    lambda <<- t$lambda; mu <<- t$mu; psi <<- t$psi
    lambda_a <<- t$lambda_a; report_prob <<- t$report_prob
    lShifts <<- mShifts <<- pShifts <<- c(0, t$age - rev(cov_breaks))
    bins <<- cov_bins
  }
  redraw <- function() {
    sim_drawn <<- sim_drawn + 1
    use_theta(draw_theta(lambda_a))
  }

  # run until we fulfill conditions
  while (!cond) {
    # a draw whose expected tree exceeds the ceiling is refused before it is simulated
    if (sbc && is.finite(MAX_LINEAGES)) {
      while (expected_lineages(lambda, mu, tMax, lShifts) > MAX_LINEAGES) {
        sim_screened <<- sim_screened + 1
        redraw()
      }
    }

    # run BD simulation - make sure we get 10+ species
    sim <- bd_sim_bounded(n0, lambda, mu, tMax, lShifts, mShifts, nFinal)

    # a draw that overshot its expectation, or outran the clock, is refused before it costs
    # a fossil sampling pass
    if (sbc && (is.null(sim) || length(sim$TS) > MAX_LINEAGES)) {
      redraw()
      next
    }

    # run fossil sampling
    fossils <- suppressMessages(sample.clade(sim, psi, tMax,
                                             rShifts = pShifts,
                                             bins = bins,
                                             returnAll = TRUE))

    # split species at anagenetic events before any of them is reported
    if (lambda_a > 0) {
      ana <- add_anagenesis(sim, fossils, lambda_a)
      sim <- ana$SIM
      fossils <- ana$FOSSILS
    }

    # drop what the boundary hides, and recompute who is alive at it
    if (PRESENT > 0) {
      fossils <- fossils[fossils$SampT > PRESENT, , drop = FALSE]
      alive <- is.na(sim$TE) | sim$TE < PRESENT
      if (nrow(fossils) > 0) {
        fossils$Extant <- alive[as.integer(sub("^t", "", fossils$Species))]
      }
      sim$EXTANT <- alive
    }
    # SBC fix: if no species were sampled, fall through and let the condition
    # below redraw the whole replicate, rather than resampling fossils on the
    # same tree (which can loop forever for tiny trees once the diversification
    # floor is removed). The species loop below uses seq_along to handle 0 safely.

    # an SBC rep reports the same record three ways; an accuracy rep just once
    models <- if (sbc) sbc_models else model
    specs <- rep(list(data.frame(matrix(nrow = 0, ncol = 4))), length(models))
    names(specs) <- names(models)
    ranges <- data.frame(matrix(nrow = 0, ncol = 3))

    # loop through species (seq_along handles the 0-sampled case safely)
    sp_list <- unique(fossils$Species)
    for (i in seq_along(sp_list)) {
      # which species is this
      sp <- sp_list[i]
      
      # get vector of occurrences for that species
      occs <- fossils[fossils$Species == sp, -which(colnames(fossils) == "SampT")]
      
      # and get true occurrence times
      true_occs <- fossils[fossils$Species == sp, c("Species", "Extant", "SampT")]
      
      # check if it is extant
      ext_sp <- fossils$Extant[fossils$Species == sp][1]
      
      # report this species' record under each model
      for (m in seq_along(models)) {
        specs[[m]] <- rbind(specs[[m]], retain_occs(occs, models[m], report_prob))
      }

      # get true range
      range <- c(max(true_occs$SampT),
                 ifelse(sum(true_occs$Extant) > 0, PRESENT, min(true_occs$SampT)))
      
      # add to ranges
      ranges <- rbind(ranges, c(sp, range))
    }
    
    # reorder and rename columns
    for (m in seq_along(specs)) {
      specs[[m]] <- specs[[m]][, c(1, 4, 3, 2)]
      colnames(specs[[m]]) <- c("taxon", "min_age", "max_age", "status")

      # change status column to extant and extinct
      specs[[m]]$status <- c("extinct", "extant")[specs[[m]]$status + 1]
    }

    # An extant lineage is seen at the present with probability rho, and that observation is
    # the status flag rather than an occurrence. A survivor that is not seen has no present-day
    # observation at all: its record stops at its youngest fossil and reads as extinct. One
    # that left no fossil either is absent from the record entirely, below.
    rho_seen <- rep(TRUE, length(sim$TS))
    if (RHO < 1) rho_seen <- runif(length(sim$TS)) < RHO
    if (RHO < 1) {
      for (m in seq_along(specs)) {
        k <- as.integer(sub("^t", "", specs[[m]]$taxon))
        specs[[m]]$status[specs[[m]]$status == "extant" & !rho_seen[k]] <- "extinct"
      }
    }

    # the reporting models keep the same species, so any of them fixes the taxon set
    specimens <- specs[[1]]
    
    # name columns for ranges
    colnames(ranges) <- c("taxon", "first_age", "last_age")
    
    # conditions (depending if it's an SBC sim or not)
    if (sbc) {
      # SBC fix: accept on >=1 sampled species, matching the analysis'
      # condition="sampling". The previous 5-50 window conditions on the number
      # of sampled species, which no dnFBDRMatrix condition can correct for and
      # which biases the estimated rates (selecting for high-diversification trees).
      # t1 is the origin lineage (paleobuddy species 1 starts at tMax, and anagenesis leaves its
      # oldest segment at index 1). Requiring it makes the oldest SAMPLED birth equal age, which is
      # what the analysis puts its origin prior on. Rejection redraws theta below, so conditioning
      # on a function of the reported record leaves the posterior alone.
      # specimens holds only fossil-sampled taxa here; the extant singletons are added after this
      # loop, and rho = 1 means every extant lineage is sampled, so t1 counts either way
      ntax <- length(unique(specimens$taxon))
      if (TAXA_COUNT == "record")
        ntax <- length(union(unique(specimens$taxon), sprintf("t%d", which(sim$EXTANT & rho_seen))))
      cond <- ntax >= MIN_TAXA && ntax <= MAX_TAXA &&
              length(unique(specimens$taxon)) >= MIN_FOSSIL_TAXA &&
              ( origin_sampled == FALSE ||
                "t1" %in% specimens$taxon || isTRUE(sim$EXTANT[1]) )

      # if cond is false, redraw everything
      if (!cond) redraw()
    }
    else {
      cond <- length(unique(specimens$taxon)) > 10 &&
        length(unique(specimens$taxon)) < 500 &&
        length(unique(specimens$taxon)) / length(sim$TS) < 0.9 &&
        length(unique(specimens$taxon)) / length(sim$TS) > 0.2
    }
  }

  # A survivor seen at the present that left no fossil still belongs in the record, as a
  # present-day (0,0) tip named t{k} to match its sim lineage index. A lone tip reports
  # identically under every model. Survivors rho missed are simply not here.
  if (sbc) {
    reported <- specs[[1]]$taxon
    for (k in which(sim$EXTANT & rho_seen)) {
      nm <- paste0("t", k)
      if (nm %in% reported) next
      for (m in seq_along(specs)) {
        specs[[m]] <- rbind(specs[[m]], data.frame(taxon = nm, min_age = PRESENT,
                            max_age = PRESENT, status = "extant", stringsAsFactors = FALSE))
      }
      ranges <- rbind(ranges, data.frame(taxon = nm, first_age = PRESENT, last_age = PRESENT,
                                         stringsAsFactors = FALSE))
    }
  }

  # paleobuddy indexes intervals oldest-first because its shifts run forward from the origin.
  # Everything downstream indexes them youngest-first, so truth_of reverses, once.
  true_vals <- if (sbc) truth_of(th, lambda, mu, psi, age) else NULL

  # return sim, ranges and k
  return(list(SIM = sim, SPECIMENS = specs, RANGES = ranges, 
              TV = true_vals, BINS = bins))
}

# simulate one set
simulate_set <- function(n_key, reps, rates, age, base_dir,
                         model = 2, unc = FALSE, extant_singletons = FALSE,
                         diff_bins = 0, sbc = FALSE, lambda_a = 0, origin_sampled = FALSE) {
  # outside the SBC sets, rates are prepared as normal
  if (!sbc) {
    # create shifts list
    shifts <- vector("list", 3)
    
    # and number of stages per rate
    n_stages <- c()
    
    # iterate through rates to see how many shifts per rate
    for (i in 1:length(rates)) {
      # check number of stages for this rate
      n_stages <- c(n_stages, length(rates[[i]]))
      
      # make shifts vector
      if (n_stages[i] > 1) { 
        shifts[[i]] <- seq(0, age, age / n_stages[i])[-(n_stages[i] + 1)]
      }
    }
    
    # get bins to bin fossil data on
    bins <- seq(0, age, age / 3)
    if (max(n_stages) == 2) bins <- seq(0, age, age / 2)
    
    # if we're on the case where we need different bins
    if (diff_bins) {
      # set those bins
      bins <- seq(0, age, age / diff_bins)
    }
  }
  
  # create vectors for number of species, sampled, and percentage sampled
  n_sp <- c()
  n_sampled <- c()
  perc_sampled <- c()
  
  # create list for sim, fossil ranges, and k
  sims <- vector("list", reps)
  
  # reuse existing seeds when present so a re-run reproduces the same trees and
  # rates (only the retention/output changes); otherwise draw and save fresh ones
  if (file.exists(paste0(base_dir, "seeds.RData"))) {
    load(paste0(base_dir, "seeds.RData"))
  } else {
    # create seeds - reps*100 apart to ensure a lot of possible seeds
    if (sbc) seeds <- runif(reps, 0, 69 * reps * 100) else
      seeds <- runif(reps, (n_key - 1) * reps * 100, n_key * reps * 100)
    save(seeds, file = paste0(base_dir, "seeds.RData"))
  }
  
  # one named row per replicate, bound at the end so the names become the columns
  true_vals <- vector("list", reps)
  
  # iterate through reps
  for (rep in 1:reps) {
    # print some info
    print(paste0("key: ", n_key, " rep: ", rep, 
                 " seed: ", seeds[rep]))
    
    # set seed 
    set.seed(seeds[rep])
    
    # an SBC set draws its own values. AGE_MIN must exceed the oldest breakpoint so
    # that every replicate has lineages in all intervals.
    if (sbc) {
      th <- draw_theta(lambda_a)
      age <- th$age
      lambda <- th$lambda; mu <- th$mu; psi <- th$psi; lambda_a <- th$lambda_a

      # make rates list
      rates <- list(lambda = lambda,
                    mu = mu,
                    psi = psi)
      
      # shifts list
      shifts <- rep(list(c(0, age - rev(cov_breaks))), 3)
    }
    
    # run sim
    sim_rep <- simulate_rep(rates, age, shifts,
                            report_prob = if (sbc) th$report_prob else 1,
                            model = model, unc = unc,
                            extant_singletons = extant_singletons,
                            sbc = sbc, lambda_a = lambda_a,
                            origin_sampled = origin_sampled,
                            theta = if (sbc) th else NULL)
    
    # an SBC set shares one timeline and adds to true_vals
    if (sbc) {
      # write timeline. The bins no longer follow the rate intervals, so the
      # timeline is the breakpoints rather than the interior bin edges.
      smart_dir_create(paste0(base_dir, "/times"))
      write.table(t(cov_breaks), 
                  paste0(base_dir, "/times/times_", rep, ".tsv"),
                  col.names = FALSE, row.names = FALSE, 
                  quote = FALSE, sep = "\t")
      
      # add true values to true_vals
      true_vals[[rep]] <- sim_rep$TV
    }
    
    # get sim, specimens and ranges
    sim <- sim_rep$SIM
    specs <- sim_rep$SPECIMENS
    specimens <- specs[[1]]
    ranges <- sim_rep$RANGES

    # calculate numbers of interest
    n_sp <- c(n_sp, length(sim$TS))
    n_sampled <- c(n_sampled, length(unique(specimens$taxon)))
    perc_sampled <- c(perc_sampled, n_sampled[rep]/n_sp[rep])
    
    # save specimens. An SBC rep reports the same record under every model, so
    # each gets its own directory alongside the one shared timeline and true values.
    for (m in seq_along(specs)) {
      spec_dir <- if (sbc) paste0(base_dir, "/specimens_", names(specs)[m]) else
                                paste0(base_dir, "/specimens")
      smart_dir_create(spec_dir)
      write.table(specs[[m]], paste0(spec_dir, "/taxa_", rep, ".tsv"),
                  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
    
    # save true ranges
    smart_dir_create(paste0(base_dir, "/true_ranges"))
    write.table(ranges, paste0(base_dir, "true_ranges/ranges_", rep, ".tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # append to list of simulations
    sims[[rep]] <- sim
  }
  
  # save sim list
  save(sims, file = paste0(base_dir, "sim_list.RData"))
  
  # outside an SBC set
  if (!sbc) {
    # write timeline
    write.table(t(bins[-c(1, length(bins))]), 
                paste0(base_dir, "times.tsv"),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # get highest number of rates for this set
    max_len <- max(unlist(lapply(1:length(rates), 
                                 function(x) length(rates[[x]]))))
    
    # get make true rates by expanding rates that are not of len max_len
    true_rates <- lapply(1:length(rates), function(x) {
      if (length(rates[[x]]) < max_len) {
        rep(rates[[x]][1], max_len)
      } else {
        rates[[x]]
      }
    })
    
    # name them
    names(true_rates) <- c("lambda", "mu", "psi")
    
    # make it a data frame
    true_rates <- as.data.frame(true_rates)
    
    # write it to file
    write.table(true_rates, paste0(base_dir, "true_vals.tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    # truth_of named every column when it built the row, so binding them is the whole job
    true_vals <- do.call(rbind, true_vals)
    
    # save true_vals
    write.table(true_vals, 
                paste0(base_dir, "/true_vals.tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  # The inference reads these rather than rederiving them, so the two cannot disagree about
  # the prior's correlation structure.
  if (OU) {
    write.table(t(OU_LAGS), paste0(base_dir, "/ou_lags.tsv"),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  # make numbers into a data frame
  nums <- data.frame(n_sp = n_sp, n_sampled = n_sampled,
                     perc_sampled = perc_sampled)
  
  # return numbers of interest
  return(nums)
}


# create function for simulating an SBC set
simulate_batch <- function(reps, reps_dir, lambda_a = 0, origin_sampled = FALSE) {
  # run simulation set
  nums <- simulate_set("cov", reps, NULL, NULL, reps_dir,
                       2, FALSE, FALSE, 0, sbc = TRUE, lambda_a = lambda_a,
                       origin_sampled = origin_sampled)
  
  # make nums a data frame
  nums <- as.data.frame(nums)
  
  # name nums data frames
  colnames(nums) <- c("n_sp", "n_sampled", "perc_sampled")
  
  # return the data frames
  return(nums)
}

###
###
# run simulations
#
# One pass reports each simulated record under every reporting model, so this fills
# specimens_{complete,firstlast}/ alongside the one shared timeline,
# true_vals.tsv, true_ranges/ and sim_list.RData.
#
# seeds.RData is reused when present, so reruns of a batch stay comparable.
#
# Driven entirely by the batch config.sh, which bin/run.sh exports; BATCHDIR and NREPS
# come from there too. Sourcing with sim_functions_only set loads the functions only.
if (!exists("sim_functions_only")) {
  reps <- as.integer(cfg("NREPS"))
  reps_dir <- Sys.getenv("BATCHDIR")
  if (reps_dir == "") stop("BATCHDIR is not set; bin/run.sh must export it")
  if (!grepl("/$", reps_dir)) reps_dir <- paste0(reps_dir, "/")
  smart_dir_create(reps_dir)

  lambda_a <- cfg("LAMBDA_A")
  cat("simulating", reps, "reps ->", reps_dir, "\n")
  nums <- simulate_batch(reps, reps_dir, lambda_a = lambda_a,
                         origin_sampled = Sys.getenv("ORIGIN_SAMPLED") == "true")

  write.table(nums, paste0(reps_dir, "nums.tsv"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  # what the draws actually used, for the manifest to check against config.sh
  keys <- c("NINTERVALS","INTERVAL_WIDTH","LMEAN","LSD","MMEAN","MSD","PMEAN","PSD",
            "AGE_MIN","AGE_MAX","BIN_WIDTH","BIN_MAX","NREPS","LAMBDA_A")
  vals <- c(NINTERVALS, INTERVAL_WIDTH, LMEAN, LSD, MMEAN, MSD, PMEAN, PSD,
            AGE_MIN, AGE_MAX, BIN_WIDTH, BIN_MAX, reps, lambda_a)
  if (is.finite(MAX_LINEAGES)) {
    keys <- c(keys, "MAX_LINEAGES"); vals <- c(vals, MAX_LINEAGES)
  }
  if (is.finite(SIM_TIMEOUT)) {
    keys <- c(keys, "SIM_TIMEOUT"); vals <- c(vals, SIM_TIMEOUT)
    cat("draws refused for outrunning the clock:", sim_timeouts, "\n")
  }
  if (is.finite(MAX_LINEAGES)) {
    tot <- sim_drawn + reps
    cat(sprintf("draws refused by the expected-clade screen: %d of %d (%.2f%%)\n",
                sim_screened, tot, 100 * sim_screened / max(tot, 1)))
    writeLines(sprintf("screened\t%d\nscreened_of\t%d\nscreened_frac\t%s",
                       sim_screened, tot, signif(sim_screened / max(tot, 1), 4)),
               paste0(reps_dir, "screened.tsv"))
  }
  if (OU) {
    keys <- c(keys, "OU_HALFLIFE", "LANCH", "TANCH", "PANCH")
    vals <- c(vals, OU_HALFLIFE, LANCH, TANCH, PANCH)
    if (OU_HIER) {
      keys <- c(keys, "OU_SD_MEAN", "OU_SD_SD", "OU_PSI_MEAN", "OU_PSI_SD")
      vals <- c(vals, OU_SD_MEAN, OU_SD_SD, OU_PSI_MEAN, OU_PSI_SD)
    }
  }
  writeLines(sprintf("%s\t%s", keys, vals), paste0(reps_dir, "sim_params.tsv"))
  cat("done\n")
}
