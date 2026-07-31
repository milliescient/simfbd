# Skyline FBD-range SBC: what's behind the miscalibration

*Draft for the FBDR group — June, 2026. Summarizes a diagnostic dig into the
skyline coverage SBC that looked mis-calibrated (ψ especially).*

## TL;DR
- **The skyline FBD-range likelihood is correct.** Verified two independent ways
  (a from-scratch Python reimplementation that matches RevBayes to ~1e-4, and an
  SBC pipeline built on the *validated* FossilSim/TreeSim simulators that calibrates
  when the analysis is specified consistently).
- **The miscalibration comes from the validation *pipeline*, not the model** — a bug
  in the SBC diagnostic script plus a few simulation/analysis specification
  mismatches.
- **paleobuddy's core looks fine** (time direction, budding, Poisson sampling).
- One genuine *residual* (the oldest interval's rates) is a **tooling gap**: there is
  no age-conditioned skyline tree simulator, so the simulation conditions on
  something the likelihood can't exactly correct for.

## How we checked the likelihood
1. **Independent Python reimplementation** of the `dnFBDRMatrix` complete-data
   log-density (the `PiecewiseConstantRates` math). It matches `rb`'s
   `lnProbability` to ~1e-4 across constant & skyline, equal & unequal rates, 1–3
   species, every interval-crossing pattern — including the γ coexistence term and
   `condition="time"`. So the C++ faithfully implements the published density; no
   coding error in the crossing terms.
2. **Independent SBC** using **FossilSim + TreeSim** (not paleobuddy) through the
   *same* `rb` analysis. With matched specification it calibrates (0/9 rank-SBC
   failures); the failures only appear when we deliberately misspecify.

## The actual causes

**1. BUG in `skyfbdr_SBC.Rmd` — skyline interval index not reversed.**
paleobuddy runs forward-time (t=0 = origin), so its true index 1 = **oldest**
interval; RevBayes `timeline` indexes 1 = **youngest**. Bruno's `post_processing.R`
correctly reverses (`c(3,2,1,…)`); David's Rmd compares `lambda[1]` to true
`lambda1` in natural order. Under the wrong mapping the SBC fails 9/9 (λ-oldest
p≈6e-60); under the correct mapping 6/9 — so a chunk of the alarming picture was a
**diagnostic artifact**, and which stage looked worst was scrambled. *(Fixed; the
reusable diagnostics take a `direct|reversed` flag.)*

**2. MISSPEC — `complete=false` used on complete data.** `complete=false` integrates
over κ *unobserved* fossils; it's the model for **stratigraphic-range** data where
the per-species count is dropped. If every fossil is kept (count known) but analyzed
`complete=false`, the model invents phantom fossils → **ψ over-estimated 40–90%**.
Evidence: identical FossilSim sims, ψ rank-SBC **0/3 fail** (exact ages +
`complete=true`) vs **3/3 fail, p≈1e-59** (all fossils + `complete=false`). For
genuine range data (oldest+youngest kept) `complete=false` *is* appropriate and
calibrates (modulo #4).

**3. MISSPEC — conditioning.** Analyses used `condition="time"`, but the data are
selected (we only keep sampled/surviving clades). The matched choice is
`condition="sampling"` (÷ P(≥1 observation)) or `condition="survival"` (÷ P(≥1
extant)), depending on how the dataset is defined. `time` + selection biases all
rates. **This is worth re-testing on the real coverage data** — just flip the flag.

**4. RESIDUAL — deepest-interval rates = a simulator limitation, not the model.**
Taxa-conditioned simulation (TreeSim `sim.rateshift.taxa`, and paleobuddy needing
extant tips / skipping μ>λ draws) selects for **deep-lineage survival**, which shifts
the *true* oldest-interval λ above the prior (accepted mean 0.187 vs prior 0.153;
P(λ>μ in oldest)=0.85 vs 0.5). The analysis uses the unconditioned prior, so the
oldest interval is mildly mis-calibrated. **No `rb` conditioning exactly matches
"n extant tips,"** so this can't be fully corrected without an **age-conditioned
skyline tree simulator** — which doesn't appear to exist (TreeSim only does
taxa-conditioned rate-shift trees; paleobuddy is age-conditioned but is the thing we
were testing). *This is the one open methodological question — @Tanja, is an
age-conditioned skyline FBD simulator tractable?*

**5. paleobuddy core checks out.** `make.rate` (time direction), `bd.sim`
(budding/asymmetric), and `sample.time` (Poisson sampling) all match the model on
inspection — not implicated as a simulator bug. (The earlier "paleobuddy looked
worse than FossilSim" gap is explained by the diagnostic bug + the
specification/conditioning differences above, not by a paleobuddy error.)

## Recommendations
1. Fix the index reversal in `skyfbdr_SBC.Rmd` (or use the new flagged diagnostics).
2. Match the data form to the flag: `complete=true` for full occurrence data;
   `complete=false` only for stratigraphic ranges.
3. Use `condition="sampling"` (or `"survival"`) to match how datasets are selected —
   not `"time"`. Re-run the real coverage SBC with this; expect it to calibrate
   (modulo #4).
4. Drop arbitrary species-count windows (they're a silent conditioning).
5. Treat the oldest-interval residual as a known limitation pending an
   age-conditioned simulator.

## What we built (in `simfbd`, refactored)
- A reusable SBC + misspecification pipeline on FossilSim/TreeSim covering the
  current `dnFBDRMatrix` model set: `bin/{sim.r, infer.Rev, pyrate.sh, run.sh,
  sbc_diagnostics.{py,R}}`, one analysis per `sims/<name>/`, `docs/`.
- The independent Python likelihood (cross-check).
- `docs/MISSPECIFICATIONS.md` — the detailed catalog behind this summary.

*Caveat on scope: the SBC numbers above are from our FossilSim re-simulations, not a
re-run of the original paleobuddy coverage set. The likelihood verification is
independent of any simulator. The single highest-value next step is to re-run the
existing coverage analysis with the corrected `condition` (#3) and the fixed
diagnostic (#1).*
