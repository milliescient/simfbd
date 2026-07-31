# Where the existing FBD-range validation pipelines misspecify the model

Catalogue of issues found while diagnosing the skyline FBD-range SBC
miscalibration (June 2026). Each entry is tagged:

- **BUG** — a genuine error in a script; fix it.
- **MISSPEC** — the model/likelihood is correct, but the pipeline applies it to
  data generated under different assumptions (so SBC fails through no fault of the
  model).
- **LIMITATION** — a tooling gap we currently cannot remove.

Evidence comes from an independent SBC pipeline built on the *validated* FossilSim/
TreeSim simulators (this pipeline), cross-checked against an independent Python
reimplementation of the dnFBDRMatrix likelihood (matches rb to ~1e-4).

---

### 1. Skyline interval index not reversed in David's SBC — **BUG**
`skyfbdr_simstudy/skyfbdr_SBC.Rmd` compares rb `lambda[1]` to true `lambda1` in
natural order. But paleobuddy runs forward-time (t=0 = origin), so its true index
1 = **oldest** interval, while RevBayes' `timeline` indexes 1 = **youngest**. They
are reversed. Bruno's `post_processing.R` correctly reverses (`c(3,2,1,...)`);
David's Rmd does not.
- *Effect:* end intervals (1,3) compared to the wrong stage → spuriously inflates
  the apparent miscalibration (direct mapping fails 9/9 SBC vs 6/9 reversed; λ3
  catastrophic p=6e-60 only under direct).
- *Fix:* reverse the true columns before matching (patched in the cloned Rmd; the
  reusable, flagged version is `sbc_diagnostics.R` / `.py` with `MAPPING`).
- Note simfbd writes youngest-first, matching rb → use `direct`.

### 2. `condition="time"` on sampled-only data — **MISSPEC**
Both pipelines keep only species that left fossils, then analyse with
`condition="time"`, which conditions on the origin time but **not** on the process
leaving any sample. The matched choice is `condition="sampling"` (rb subtracts
`log(1 - p_sampled)`). Pairing "observe only sampled lineages" with `"time"` is a
conditioning mismatch that biases all rates.
- *Fix:* use `condition="sampling"`. (Default in `config.sh`.) Worth re-running
  Bruno's real analysis with this flag.

### 3. Sampled-species count window (5–50) — **MISSPEC**
Rejecting sims outside a species-count window conditions on sample size, which the
likelihood does not model (any window is a silent conditioning).
- *Fix:* require only `>= NMIN` (enough to have data), no upper cap; pair with
  `condition="sampling"`.

### 4. `complete=false` requires INCOMPLETE data — **MISSPEC (usage rule)**
`complete=false` integrates over κ *additional, unrecorded* fossils. It is correct
only for **stratigraphic-range** data where the per-species count is dropped
(e.g. keep oldest+youngest occurrence). If every fossil is kept (count known),
`complete=false` is misspecified and the model invents phantom fossils → **ψ
over-estimated**.
- *Evidence:* identical FossilSim sims, ψ rank-SBC: **0/3 fail** with
  exact-ages+`complete=true`; **3/3 fail, ψ bias +40–90%** with all-fossils-binned
  +`complete=false`. The only differences were exact→binned and complete flag.
- *Rule:* `complete`→`complete=true`; `incomplete`→`complete=false`. Encoded in the
  pipeline (`MODEL` ↔ `complete`).

### 5. `rho=0` while data contain extant species — **MISSPEC (minor)**
Coverage sims include extant species (`status=extant`) but were analysed with
`rho=0` (all lineages treated as extinct). Tested: re-running paleobuddy data with
`rho=1`+status changed only the present-ward interval (λ3 bias 0.14→0.03 but ψ3
0.17→0.30) — mean |bias| unchanged (0.160 vs 0.163). So not the main driver, but
still a misspecification.
- *Fix:* set `rho` to match the simulation (1 if extant lineages are observed at
  the present); ensure `status` is honoured.

### 6. Deepest-interval λ from taxa-conditioned simulation — **LIMITATION**
There is no age-conditioned *skyline* tree simulator (TreeSim only offers
taxa-conditioned `sim.rateshift.taxa`; paleobuddy is age-conditioned but is the
thing under test). Conditioning on n extant tips (and skipping μ>λ draws that
cannot reach the target) selects for deep-lineage survival, which shifts the true
oldest-interval λ above the analysis prior (accepted λ_oldest mean 0.187 vs prior
0.153; P(λ>μ in oldest)=0.85 vs 0.5). The analysis uses the unconditioned prior →
oldest-interval λ mildly under-estimated (p≈2e-3, bias −0.09).
- *Status:* unavoidable with current tooling (simfbd's own sim.r has the same
  taxa-conditioning + 5 s timeout). Ask whether an age-conditioned skyline FBD
  simulator is feasible (Tanja). Everything else calibrates.

### 7. `binary` ("model 3" presence/absence) not in current rb build — **NOTE**
Older templates (`m3_template.Rev`) pass `binary=true` to `dnFBDRMatrix`, but the
current `fbdr`-branch build has no such argument. The current model set is
`{complete, incomplete} x {FBDR, BDS} x {constant, skyline}` with
`complete`, `rho`, `condition`, `resample`. Presence/absence is not testable here
without restoring that code path.

---

## Net conclusion
The skyline FBD-range **likelihood is correct** (independent FossilSim data are
well-calibrated under the matched configuration, and the Python reimplementation
matches rb term-for-term). The reported miscalibration came from pipeline
mis-application — chiefly the `complete=false`-on-complete-data misspecification
(#4) and conditioning mismatches (#2,#3) — plus the David-Rmd index bug (#1) that
distorted the *diagnosis*. paleobuddy's core (time direction, budding, Poisson
sampling) checked out; it is not implicated as a simulator bug.
