# Pilot: same process as batch1 at 200 reps, for a quick look. skyline FBD range process, budding speciation, extended tree.
#
# This file is the only definition of the generating process. One pass emits every
# reporting model from one record, so the analyses under runs/ share a timeline and
# one true_vals.tsv. They inherit these hyperpriors and may not redeclare them.

# skyline structure: NINTERVALS rates, breakpoints at INTERVAL_WIDTH before the present
NINTERVALS=3
INTERVAL_WIDTH=3

# rate hyperpriors: lambda, mu ~ Lognormal(LMEAN, LSD); psi ~ Lognormal(PMEAN, PSD)
LMEAN=-2.0; LSD=0.5
MMEAN=-2.0; MSD=0.5
PMEAN=-0.125; PSD=0.5

# origin ~ Uniform(AGE_MIN, AGE_MAX). AGE_MIN must exceed the oldest breakpoint
# ((NINTERVALS-1)*INTERVAL_WIDTH) so every replicate has lineages in all intervals.
AGE_MIN=9; AGE_MAX=15

# occurrence age bins, kept independent of the origin: anything that scales with age
# hands the analysis the origin it is meant to estimate
BIN_WIDTH=3; BIN_MAX=30

# anagenetic speciation rate; 0 is budding only
LAMBDA_A=0

# require the origin lineage to be sampled
ORIGIN_SAMPLED=false

NREPS=200
