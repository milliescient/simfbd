# The final SBC design

One batch of simulated histories, fitted four ways, answering three questions at once:
whether `dnFBDRP` is calibrated on a complete record, whether it is calibrated on a
first/last record, and how its estimates compare with PyRate's likelihood and with the
classical paleontological rate estimators on the same data. Every arm reads the same
`true_vals.tsv`, so every comparison is paired on identical truth.

## The prior

A stationary AR(1) on the log rates, which on an equal-width grid is an OU process observed
at the interval boundaries: with lag-1 correlation `phi` on intervals of width `dt` the
mean reversion is `-log(phi)/dt`. `phi` is fixed, because estimating it on eight intervals
would be estimating a property of the discretisation rather than of the process.

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

OU_PHI=0.5                  # the one hyperparameter that stays fixed
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

Measured at the settings above, over 11,013 draws: **0.45% screened, 0 timeouts**. Raising
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

## One binary

Every arm runs on one `rb`, built from `fbdr-estimators` = `fbdr` merged with
`paleo-rate-estimators`, because no earlier build has both `dnBDS(origin_prior=)` and the
three estimator functions. Build it with `BUILD_DIR=<dir> ./build.sh`, never bare make, and
point `RBIN` at it.
