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
  simfbd             the front door: simfbd extend | fit | status | summarize | ranks | --version
  run.sh             driver behind `simfbd fit`: bash bin/run.sh <batch>[/runs/<name>]
  extend.sh          behind `simfbd extend`: bash bin/extend.sh <batch> <how many>
  manifest.py        writes, reads and checks manifest.yml, and says what every config key means
  summarize.R        SBC scoring and figures, driven by the run's manifest
  bias.R             dnBDS against dnFBDRP on identical data, by sampled fraction
  env.sh             toolchain paths; monitor.sh, kill.sh generic helpers
sims/batch1/
  config.sh          the generating process and every hyperprior
  manifest.yml       what was actually used and how many replicates exist, written by bin/run.sh
  specimens_complete/  specimens_firstlast/
  times/  true_ranges/  true_vals.tsv  sim_list.RData  seeds.RData
  runs/fbdrp_firstlast/
    config.sh        which reporting model to read, and how to fit it
    manifest.yml     + the batch hash it ran against, and the rb identity
    output/  aux/  results/
bin/pyrate/          PyRate and BDS inference scripts; sampling/ holds their prior.rev
```

## Running

```sh
simfbd extend sims/batch1 200                # simulate the batch's first 200 replicates
simfbd fit sims/batch1/runs/fbdrp_firstlast  # infer on every replicate not yet fitted
simfbd summarize sims/batch1/runs/fbdrp_firstlast
simfbd extend sims/batch1 800                # later: 800 more, and the batch has 1000
simfbd status sims/batch1 sims/batch1/runs/fbdrp_firstlast
```

`simfbd` is `bin/simfbd`: call it by its path, or link it onto `PATH`, and it finds the rest of
`bin/` beside itself. The scripts behind it can still be called directly.

A batch is a directory holding a `config.sh`, and it can be anywhere on disk. A relative path is
looked for where the command is typed and then under this repository, so `sims/batch1` means what
it always has. Nothing in a manifest records where a batch is: a run finds its batch as the
directory two above it. So a study can be kept, moved or handed over as a directory, with this
repository at a stated commit to read it.

`simfbd --version` prints that commit, in the words a manifest records it in, ending in `-dirty`
when the code has uncommitted changes.

### Checking a batch or a run

A manifest records what a notebook's lock records, in the same words. `uses` names the code that
ran as `name@version`, and `inputs` and `outputs` give each path with its hash:

```yaml
uses: ["simfbd@f0c6fe5", "rb@28e5b2"]
inputs:
  specimens_complete: a108a52a6513    # read from the batch
outputs:
  output: caf92bd48090                # written by this run
```

`simfbd status` recomputes the hashes and reports each path as `current`, `changed` or `missing`.
For a run it also checks the edge: what the run read must be what the batch says it wrote, and a
run fitted before its batch was extended is reported `stale`. It exits non-zero if anything is
changed, missing or stale. The hash is `digest()` in `bin/manifest.py`: SHA-1 over each file's
relative name and bytes, files in sorted order, first 12 hex digits. It is the hash sloop
computes, so a page that pins one of these paths records the same value.

A batch's manifest keeps the output hashes of every extension, and only the last can still be
checked, since an extension changes the files. A batch or run written before hashes were recorded
is reported as `unrecorded`, and one with a `manifest.tsv` has nothing to check.

A batch's `config.sh` does not say how many replicates it has; the manifest counts them.
`bin/extend.sh` reads the batch's own config, so new replicates cannot be drawn under
different settings, and it refuses if the config or `bin/sim.R` has changed since the earlier
ones. Those are never redrawn or rewritten: an extension simulates only what is missing,
appends rows to `true_vals.tsv` and `nums.tsv`, and adds its per-replicate files beside the
others. Anything that pinned a hash of `true_vals.tsv`, a specimens directory,
`sim_list.RData` or the manifest will see that hash move, because those did grow.

A batch simulated before this has `NREPS` in its config and a `manifest.tsv`. It runs as it
always did, its manifest is read by the same scripts, and it cannot be extended.

`sims/smoke/` is the same thing at 6 replicates and short chains, for checking the
pipeline end to end in about a minute.

Batches that differ only in the origin prior pool into a range of tree sizes. A deep
origin wants `MAX_LINEAGES`, a ceiling on lineages ever born, because it leaves the oldest
interval unbounded and a draw in the upper tail of lambda then never finishes simulating.

`SIM_MAX_DRAWS` is the second ceiling, on the work one simulation may do. It counts the
waiting times the birth-death simulator draws, across every tree it tries, and refuses the
draw when the count is passed. The count follows from the random stream, so the same seeds
refuse the same draws on any machine and a batch reproduces from `seeds.RData`. It replaces
`SIM_TIMEOUT`, which counted seconds and so refused whatever the machine was too busy for;
a new batch that sets `SIM_TIMEOUT` is refused. A tree of n lineages takes about 3n draws,
so keep the ceiling several times `MAX_LINEAGES`: 100000 cuts a runaway off after roughly
what 20 seconds used to, 200000 after roughly 60. The simulator prints how many draws were
refused and the most any finished simulation used.

That ceiling is the one filter that costs something. It is applied by refusing the draw
before it is simulated, so it is a function of theta alone: it truncates the prior the
simulator draws from while the analysis keeps the untruncated one, and the SBC ranks are
non-uniform by exactly that much. A floor and ceiling on the number of reported taxa is a
function of the record and is free, because conditioning on a function of the record leaves
the posterior given that record alone. Every batch now records `screened_frac`, the fraction
of draws the ceiling refused, in its manifest.

How much it costs is a question of resolution: a truncation of fraction e shifts a mean rank
by roughly e/2, against a standard error of 0.289/sqrt(n). At n = 10000 that is 0.0029, so e
has to stay near or below 0.5% for the truncation to be invisible. Whether a given prior can
meet that is a question for the study that chooses it; `screened_frac` is what lets it be
measured rather than assumed.

## The two key axes

- **REPORTING** = `complete` | `firstlast` — which record the analysis
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
- **An edited config is refused, not mixed.** Each manifest records the hash of the
  config that produced its data, and a run records the batch hash it ran against.
  Changing either after the fact stops the run instead of quietly mixing parameters
  across replicates. `FORCE=1` overrides for a run, and for a batch with a `manifest.tsv`.

`manifest.yml` is what makes an output self-describing: it names the batch, the reporting
model, the resolved hyperpriors, the config hashes, the script hashes and the rb binary.
`bin/summarize.R` reads it rather than re-deriving anything from directory names, so a
result cannot be scored against the wrong truth.

It is also where a config key is explained. `config.sh` keeps short shell names because
`bin/run.sh` sources it, so the manifest writes each setting under a name that says what it
is, with the config key beside it:

```yaml
rate_priors:
  fossil_sampling:
    log_mean: -0.125     # PMEAN
    log_sd: 0.5          # PSD
```

The tables behind that are `BATCH_SETTINGS` and `RUN_SETTINGS` in `bin/manifest.py`. A new config key gets a line there,
with a name a reader can understand without opening `bin/sim.R`. Scripts read a manifest
through `bin/manifest.py get` and `flat`, which answer by config key for either format.

## Scope

`MODEL=fbdr` is wired. `bds` and `pyrate` remain as config keys and are rejected with a
clear error until reconnected. Their scripts are in `bin/pyrate/` and `sampling/`, and
`bin/pyrate.sh` reads the same `taxon/min_age/max_age/status` specimen format `bin/sim.R`
already writes, so reconnecting them needs no separate simulator. The previous TreeSim
pipeline is in git history at 8563117.

## Toolchain (this machine)

- `rb`: pinned build, set in `bin/env.sh` and overridable with `RBIN=...`; the manifest
  records the path, the md5 and the binary's own git describe, so it names what actually ran.
- R + paleobuddy: micromamba env `/research/phyloworks/mm/root/envs/fsim`
