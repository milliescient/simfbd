# The final SBC design

One batch of simulated histories, fitted four ways, answering three questions at once:
whether `dnFBDRP` is calibrated on a complete record, whether it is calibrated on a
first/last record, and how its estimates compare with PyRate's likelihood and with the
classical paleontological rate estimators on the same data. Every arm reads the same
`true_vals.tsv`, so every comparison is paired on identical truth.

## The prior

A stationary OU process on the log rates, sampled at the interval midpoints, with mean
reversion `OU_THETA` = 0.5 per Myr. It is fixed, because estimating it on eight intervals
would be estimating a property of the discretisation rather than of the process.

It is specified per unit time and not per interval, and that is the point. A lag-1
correlation attached to the cells makes the prior a property of the grid: refine the grid at
a fixed span and you get a different process, with more accumulated variance. Measured
against the expected-clade screen, a walk with variance per cell refuses 0.92% of draws at
four intervals and 19.80% at sixty-four, while the same walk with variance per Myr holds at
4.30% and 5.12%. Neither is usable, but only the first gets worse the finer you cut. The
same applies to a stationary field, and specifying theta per Myr is what makes the interval
count a resolution choice and nothing else. Refusal varies by less than 0.02 points over
theta between 0.4 and 0.6, so the value is not a delicate choice.

```sh
NINTERVALS=8
INTERVAL_WIDTH=1.5          # oldest shift at 10.5, spanned by 93% of records
LMEAN=-1.5; LSD=0.5         # LSD/MSD/PSD and MMEAN are inert in this arm, kept at the
MMEAN=-1.5; MSD=0.5         # medians their hyperpriors imply
PMEAN=-1.2; PSD=0.3
AGE_MIN=15; AGE_MAX=21
BIN_WIDTH=1.5; BIN_MAX=36
LAMBDA_A=0
ORIGIN_SAMPLED=false

OU_THETA=0.5                # per Myr, the one hyperparameter that stays fixed
LANCH=0.25                  # across-replicate level of log lambda
TMEAN=0; TANCH=0.2          # turnover, centred on mu = lambda
PANCH=0.85                  # across-replicate level of log psi: the sampling sweep
OU_SD_MEAN=-0.7; OU_SD_SD=0.5    # within-replicate sd, shared by lambda and turnover
OU_PSI_MEAN=-0.7; OU_PSI_SD=0.5  # within-replicate sd for psi

MIN_TAXA=10; MAX_TAXA=50; TAXA_COUNT=record; MIN_FOSSIL_TAXA=5
MAX_LINEAGES=20000; SIM_TIMEOUT=60
```

Why stationary rather than the usual GMRF: a random walk's marginal variance grows with
distance from its anchor, so the variance of the interval average that the expected-clade
screen bounds grows like `sigma^2 K / 3`, and the screen then refuses 3.7% of draws at eight
intervals and 14.8% at thirty-two. That screen is a function of theta alone, so it truncates
the prior the simulator draws from while the analysis keeps the untruncated one. A stationary
field's average concentrates instead, and the screen stops mattering. The derivation and the
measurements are notebook page 14.

Mean levels and within-replicate sds are drawn per replicate and estimated, so they are
truth like any other parameter: `m_lambda`, `m_turnover`, `m_psi`, `tau_rates`, `tau_psi`
appear in `true_vals.tsv` and are scored for rank uniformity alongside the rates.

## The grid, and why it has eight cells

There are two grids and they are different objects. The **rate grid** is the seven
breakpoints at 1.5 Myr intervals back from the present, giving eight rate intervals. The
**counting grid**, `BIN_WIDTH=1.5` out to `BIN_MAX=36`, is what the occurrence record is
binned into for the classical estimators. They share a width here and nothing else.

The rate grid is not equal-width. The breakpoints are fixed and do not scale with age,
because the analysis estimates the origin and anything it reads that scales with age hands
it the answer. The origin falls between 15 and 21, so the oldest interval runs from 10.5 to
wherever the origin is and is 4.5 to 10.5 Myr wide, against 1.5 for the other seven.

Theta per Myr handles that directly. The lags between interval midpoints are 1.5 Myr six
times and then 4.5, so the correlations are 0.472 six times and then 0.105: the oldest cell is
far from its neighbour and is coupled to it weakly, which is what an OU process on the
timeline says. The simulator writes those lags to `ou_lags.tsv` and the inference reads
them, so the two cannot disagree about the prior.

The oldest cell's width is taken at the age prior's mean rather than at each replicate's own
age. Using the replicate's age would be exact, but it would put the origin into the rate
prior, and this model has a history of origin couplings breaking calibration. Measured over
40,000 draws, the two differ by 0.015 points of refusal rate, so the coupling buys nothing.

**Eight cells, not more.** The prior is no longer a reason to choose otherwise, since theta
per Myr makes the count a pure resolution choice. The data is the reason. Measured over
400 replicates of the record:

| | |
|---|---|
| reported taxa per replicate | 17 median, 11 to 36 |
| occupied bins per replicate | 10 median |
| taxa per occupied bin | 3 median, 1 at p10 |
| bins holding exactly one taxon | 31% |
| three-timer taxa per replicate | 2 median, 0 at p10 |
| replicates with no three-timer at all | 33% |

The three-timer estimator needs a taxon sampled in three consecutive bins, and a third of
replicates already supply none. Halving the width would leave that arm measuring nothing,
and would double the rate parameters against a chain whose worst-mixing coordinate is
already the oldest interval. Because the prior is now grid-invariant, a finer run is
interpretable as the same process at higher resolution, so it stays available as a
robustness check rather than as the design.

A fully regular grid is also available and was tested: `dnFBDRP` accepts breakpoints older
than the origin, so a 1.5 Myr grid out to `AGE_MAX` with fourteen cells initialises and
runs. It is not used, because it puts the origin in the interior of the grid where it
crosses breakpoints, which is an untested regime for this likelihood and would need its own
validation to buy nothing the per-Myr parameterisation does not already give.

## Which filters are free and which are not

| filter | function of | effect |
|---|---|---|
| `MIN_TAXA`, `MAX_TAXA` on the reported record | the data | free |
| `MAX_LINEAGES` on the true lineage count | a latent variable | truncates |
| the expected-clade screen before simulating | theta | truncates |
| `SIM_TIMEOUT` | theta, in effect | truncates |

Conditioning on a function of the record leaves the posterior given that record alone, so a
floor and ceiling on reported taxa cost nothing as long as theta is redrawn jointly. The
other three truncate the prior, and the design is chosen to make them fire rarely: a
truncation of fraction `e` shifts a mean rank by about `e/2` against a standard error of
`0.289/sqrt(n)`, which is 0.0029 at n = 10000, so `e` has to stay near or below 0.5%.

Under the per-cell parameterisation, over the 111,811 draws its 10,000-replicate batch
took: **0.446% screened, 0 timeouts**. Being re-measured under theta per Myr. Raising
`MAX_LINEAGES` is what buys this, since it converts theta rejections into free data
rejections: 400 refuses 1.87%, 4000 refuses 0.81%, 20000 refuses 0.56% at the wider
`LANCH=0.4`, and 200000 starts producing timeouts instead. Every batch records
`screened_frac` in its manifest.

## What the design produces (1000 replicates)

| | |
|---|---|
| lambda fold within a replicate | 3.13 median |
| psi fold within a replicate | 3.19 median |
| lambda level across replicates | 2.5x, p10 to p90 |
| psi level across replicates | 9.3x, p10 to p90 |
| reported taxa | 17 median, 11 to 38 at p10/p90, bounded 10 to 50 |
| true lineages | 26 median, 124 at p99 |
| lineages sampled | 71% median, 43% at p10 |
| unsampled lineages, psi quartile 1 | 44% |
| unsampled lineages, psi quartile 4 | 15% |

The last two rows are why one batch serves both comparisons. The low-psi quarter is where
PyRate's complete-sampling assumption should bite, and the high-psi quarter is where the
three-timer estimator has the coverage it needs.

## The arms

| run | fits | reads | scored by |
|---|---|---|---|
| `fbdrp_complete` | `dnFBDRP`, `COND=time`, uniform origin | every occurrence | rank uniformity |
| `fbdrp_firstlast` | `dnFBDRP`, `COND=time`, uniform origin | first and last only | rank uniformity |
| `bds_pyrate` | `dnBDS`, `COND=sampling`, exponential origin, present-only dropped | every occurrence | accuracy |
| `bds_time` | `dnBDS`, `COND=time`, uniform origin | every occurrence | rank uniformity |
| `estimators` | `fnPerCapitaRates`, `fnBoundaryCrosserRates`, `fnThreeTimerRates` | binned counts | accuracy |

`bds_pyrate` is PyRate's likelihood exactly, and those three choices are also why it is not
a calibration arm: its origin prior is not the generating one, so its ranks are non-uniform
by construction and it has to be read as accuracy against known truth, stratified by psi.
`bds_time` is the control that separates a biased estimator from a miscalibrated likelihood.

## Chain length

100M generations, printing every 50,000, which retains 1500 draws after warmup. Measured on
30 replicates under the per-Myr prior: median per-replicate minimum ESS 708, 2 of 30 below
200, maximum Rhat 1.016. The per-cell prior gave 738 on 40 replicates, so decoupling the
oldest interval costs nothing in mixing.

What it does change is where the chain struggles. Under the per-cell prior the worst
coordinate was `lambda[1]` in every replicate. Under this one it is spread over `mu[6]` and
`psi[3]` to `psi[6]`, the middle and older intervals. Worth knowing before reading a low ESS
in one coordinate as a defect in that coordinate. The GMRF's figure at the same
length was 873, so the five estimated hyperparameters cost about 15% in the worst
coordinate rather than an order of magnitude, and `tau_rates` and `tau_psi` are not
themselves the bottleneck.

Those 1500 draws are autocorrelated, so summarize with `SBC_THIN_RANKS=2`. Thinning the
ranks is the fix for the tail-ESS warning, not a longer chain, and it costs nothing. If the
low-ESS replicates have to be dropped, ESS is a function of the record and the sampler seed
alone, so that filter is free in the sense of the table above.

## One binary

Every arm runs on one `rb`, built from `fbdr-estimators` = `fbdr` merged with
`paleo-rate-estimators`, because no earlier build has both `dnBDS(origin_prior=)` and the
three estimator functions. Build it with `BUILD_DIR=<dir> ./build.sh`, never bare make, and
point `RBIN` at it.
