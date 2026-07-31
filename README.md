# simfbd — FBD-range simulation & SBC validation

Simulation-based calibration and misspecification analyses for the RevBayes
`dnFBDRMatrix` model set.

## Layout

A **batch** owns the generating process and every hyperprior. One simulation pass emits
all reporting models from one record, so the analyses under it share a timeline and one
`true_vals.tsv`, and comparisons between them are paired on identical truth.

An **analysis** under `runs/` chooses only what to read and how to fit it. It inherits the
batch's hyperpriors and may not redeclare them.

```
bin/                 shared scripts; one of each, no per-analysis copies
  sim.R              paleobuddy simulator; emits every reporting model in one pass
  infer.Rev          rb inference; every prior arrives as a variable, none hardcoded
  run.sh             driver: bash bin/run.sh sims/<batch>[/runs/<name>]
  summarize.R        SBC scoring and figures, driven by the run's manifest
  env.sh             toolchain paths; monitor.sh, kill.sh generic helpers
sims/batch1/
  config.sh          the generating process and every hyperprior
  manifest.tsv       what was actually used, written by bin/run.sh
  specimens_complete/  specimens_firstlast/  specimens_truncated/
  times/  true_ranges/  true_vals.tsv  sim_list.RData  seeds.RData
  runs/fbdrp_firstlast/
    config.sh        which reporting model to read, and how to fit it
    manifest.tsv     + the batch hash it ran against, and the rb identity
    output/  aux/  results/
legacy/              the previous TreeSim pipeline, preserved verbatim
```

## Running

```sh
bash bin/run.sh sims/batch1                      # simulate the batch
bash bin/run.sh sims/batch1/runs/fbdrp_firstlast # simulate if needed, then infer
Rscript bin/summarize.R sims/batch1/runs/fbdrp_firstlast
```

`sims/smoke/` is the same thing at 6 replicates and short chains, for checking the
pipeline end to end in about a minute.

## The two key axes

- **REPORTING** = `complete` | `firstlast` | `truncated` — which record the analysis
  reads. All three come from the same simulated history.
- **INFER** = `complete` | `incomplete` — how `dnFBDRMatrix` treats it
  (`complete=true` / `false`).

Matching them is validation; mismatching them is a misspecification analysis. Each entry
in `docs/MISSPECIFICATIONS.md` corresponds to a runnable analysis.

## Why the parameters are tracked this way

The hyperpriors have to reach two places: the simulator that draws the true values, and
the rb analysis that puts priors on them. When they are written down twice they drift,
and the drift is silent because both halves still run. That has already happened twice
here: once when the origin prior was `U(9,15)` in the simulator and `U(5,15)` on record,
and once when a copied analysis read one batch's data while being scored against another
batch's truth.

So `config.sh` is the single definition. `bin/run.sh` sources it, exports it to
`bin/sim.R`, and writes the same values into the per-replicate `aux/run_<rep>.Rev` that
`bin/infer.Rev` reads. Neither script contains a prior of its own. Four rules keep it
that way:

- **A missing key is fatal.** `bin/sim.R` has no defaults, because a default is what lets
  a renamed or forgotten key produce a plausible run.
- **An analysis may not set a batch-owned key.** `bin/run.sh` refuses if one appears in a
  `runs/*/config.sh`, so there is never a second definition to disagree with.
- **The simulator reports what it drew with.** `sim_params.tsv` is compared against the
  config, so a key that fails to reach it is caught rather than assumed.
- **An edited config is refused, not mixed.** Each `manifest.tsv` records the hash of the
  config that produced its data, and a run records the batch hash it ran against.
  Changing either after the fact stops the run instead of quietly mixing parameters
  across replicates. `FORCE=1` overrides.

`manifest.tsv` is what makes an output self-describing: it names the batch, the reporting
model, the resolved hyperpriors, the config hashes, the script hashes and the rb binary.
`bin/summarize.R` reads it rather than re-deriving anything from directory names, so a
result cannot be scored against the wrong truth.

## Scope

`MODEL=fbdr` is wired. `bds` and `pyrate` remain as config keys and are rejected with a
clear error until reconnected; the previous TreeSim-based pipeline is preserved verbatim
under `legacy/`.

## Toolchain (this machine)

- `rb`: pinned build at `/research/phyloworks/rb-sim/projects/cmake/build/rb`,
  overridable with `RBIN=...`; the manifest records which one ran.
- R + paleobuddy: micromamba env `/research/phyloworks/mm/root/envs/fsim`
