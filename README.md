# simfbd — FBD-range simulation & SBC validation

Simulation-based calibration and misspecification analyses for the RevBayes
`dnFBDRMatrix` model set.

## Install

    pip install -e /path/to/simfbd        # or: uv tool install /path/to/simfbd

That puts `simfbd` on the path. It runs R and RevBayes, which are not Python's to install.
Each is found from `$RBIN` or `$RSCRIPT`, else `~/.config/simfbd/tools.yml`, else `PATH`:

```yaml
rb:                               # one build, or a list of them
  - /opt/revbayes/build/rb
  - /opt/revbayes/build-older/rb
rscript: /opt/R/bin/Rscript
```

`simfbd doctor` says which this machine would use, which rule supplied each, and whether the R
packages are there. Every `extend` and `fit` prints the same two lines before it starts, since
a machine with no tools file quietly falls back to whatever is on `PATH`. `simfbd uses` prints
what would run as `name@version`, and `simfbd uses rb` one version.

A config can say which versions it is to be run with, and is then refused where they are not
those, however the tools were found:

```yaml
uses: [rb@28e5b2]                 # in a run's config: any of simfbd, rb
uses: [paleobuddy@1.1.0, R@4.5.3] # in a batch's: any of simfbd, R, paleobuddy
```

### Which rb analyzes a run

A machine can hold several builds of rb, and each is known by the commit its binary carries,
which `simfbd doctor` lists. `simfbd fit` picks one for the run, by the first of these that
applies:

1. the version the run's config states under `uses`. If no build here is that one, the run is
   refused and the builds that are here are listed.
2. the build that analyzed the run's earlier replicates, read from its manifest, so a run that
   is topped up after its batch is extended keeps one rb throughout. If that build is gone the
   run is refused; `--force` goes on with the default and the manifest then names the new one.
3. the first build the tools file lists, else `rb` on `PATH`.

`$RBIN` overrides all three for one command, and is still held to what the config states.

So which rb analyzes a run is simfbd's to settle, not a notebook's. A notebook's own check of rb
sees the one its scripts get by calling `rb`, which need not be the one simfbd picks. The run's
manifest says which it was, as `uses: ["simfbd@...", "rb@..."]`, and a notebook can show that.

## Layout

A **batch** owns the generating process and every hyperprior. One simulation pass emits
all reporting models from one record, so the runs under it share a timeline and one
`true_vals.tsv`, and comparisons between them are paired on identical truth.

A **run** under a batch's `runs/` chooses only what to read and how to analyze it. It inherits
the batch's hyperpriors and has no way to set them.

```
a-batch/
  config.yml         the generating process and every hyperprior, written by you
  manifest.yml       the record of what was simulated from it, written by simfbd
  specimens_complete/  specimens_firstlast/
  times/  true_ranges/  true_vals.tsv  sim_list.RData  seeds.RData
  runs/fbdrp_firstlast/
    config.yml       which record to read, and how to analyze it
    manifest.yml     + what it read and wrote, and which code and which rb ran
    output/  aux/  results/
```

A batch is a directory anywhere on disk, and nothing in a manifest records where it is: a
run finds its batch as the directory two above it. So a study can be kept, moved or handed
over as a directory, with simfbd at a stated commit to read it.

```
src/simfbd/          the package
  cli.py             the simfbd command
  settings.py        every setting: the name it is written under and what it means
  defaults/          batch.yml and run.yml: every setting, with the value it takes when left out
  config.py          reads and checks a config.yml
  batch.py  fit.py   simulate replicates, analyze a run
  manifest.py        writes, reads and checks a manifest.yml
  ranks.py           mean rank per parameter
  scripts/           what simfbd runs or carries: the simulator sim.R, the Rev templates,
                     summarize.R, and the shell driver for batches made before configs were YAML
```

What it carries is only what it runs. Analyses of particular batches are in the notebook's
`scripts/`, and the classical-estimator coverage sweeps are with the Warnock probe.

## Running

```sh
simfbd extend sims/batch1 200                # simulate the batch's first 200 replicates
simfbd fit sims/batch1/runs/fbdrp_firstlast  # analyze every replicate not yet analyzed
simfbd summarize sims/batch1/runs/fbdrp_firstlast
simfbd extend sims/batch1 800                # later: 800 more, and the batch has 1000
simfbd status sims/batch1 sims/batch1/runs/fbdrp_firstlast
```

## Config and manifest

A config is what you ask for and a manifest is what was done. You write the config and
simfbd never edits it; simfbd writes the manifest and you never edit it. A setting is written
once, in the config, and the manifest holds only the record: the hash of the config it was
made from, how many replicates there are, which code ran, the hashes of what it read and
wrote, and what the filter refused. So a change to `config.yml` is always a decision, and a
change to `manifest.yml` is always an event.

```yaml
# a batch's config.yml
timeline:
  rate_intervals: 3
  interval_width: 3
origin_age:
  min: 9
  max: 15
  origin_lineage_must_be_sampled: false
rate_priors:
  speciation: {log_mean: -2.0, log_sd: 0.5}
  extinction: {log_mean: -2.0, log_sd: 0.5}
  fossil_sampling: {log_mean: -0.125, log_sd: 0.5}
  anagenetic_speciation: 0
record:
  occurrence_bins: {width: 3, oldest_edge: 30}
simulation_limits:
  max_lineages: 250
```

```yaml
# a run's config.yml
inference:
  model: fbdr
  record_read: complete
  record_treated_as: complete
  conditioned_on: time
  extant_sampling_probability: 1.0
mcmc:
  generations: 200000
  sample_every: 1000
  chains_at_once: 6
```

A config holds only what it needs to say. Every setting there is, with what it means and the
value it takes when a config leaves it out, is in `defaults/batch.yml` and `defaults/run.yml`
in the package, and `simfbd config default batch` prints one: copy it to start a config. Those
files are where the defaults come from, not a description of them. A name that is not a
setting is refused, with the nearest one that is, and so is a run setting in a batch's config
or the other way round. Times are in whatever unit the rates are per.

A batch's config does not say how many replicates it has: `simfbd extend` adds them and the
manifest counts them. New replicates cannot be drawn under different settings, since they
come from the batch's own config, and `extend` refuses if that config or the simulator has
changed since the earlier ones. Those are never redrawn or rewritten: an extension simulates
only what is missing, appends rows to `true_vals.tsv` and `nums.tsv`, and adds its
per-replicate files beside the others. Anything that pinned a hash of `true_vals.tsv`, a
specimens directory or the manifest will see that hash move, because those did grow.

### Checking a batch or a run

A manifest records what a notebook's lock records, in the same words. `uses` names the code
that ran as `name@version`, and `inputs` and `outputs` give each path with its hash:

```yaml
uses: ["simfbd@f0c6fe5", "rb@28e5b2"]    # a batch's names R and paleobuddy too
inputs:
  config.yml: d26610a2e386               # this run's config
  batch:
    config.yml: 60319e066c3d             # the batch's
    specimens_complete: 5d415912e147
  simfbd:
    infer-fbdrp.Rev: 12eb464f40e7        # the template simfbd carried
outputs:
  output: 0b1b35e37934                   # written by this run
```

Everything a batch or a run read is under `inputs`, and everything is hashed one way. The hash
is `digest()` in `manifest.py`: SHA-1 over one line per file, its relative name and the SHA-1
of its bytes, files in sorted order, first 12 hex digits. It is the hash sloop computes, so a
page that pins one of these paths, a config included, records the same value.

`simfbd status` recomputes the hashes and reports each path as `current`, `changed` or
`missing`. For a run it also checks the edge: what the run read must be what the batch says
it wrote, and a run analyzed before its batch was extended is reported `stale`. It exits
non-zero if anything is changed, missing or stale. A script of simfbd's own that has changed
since is reported as `differs` and not counted, since the commit under `uses` already names
what ran.

A batch's manifest keeps the output hashes of every extension, and only the last can still
be checked, since an extension changes the files.

`simfbd --version` prints the commit of the code, ending in `-dirty` when it has uncommitted
changes. Run from a checkout it asks git; an installed copy carries the commit it was built
from.

### Batches made before configs were YAML

Such a batch holds a `config.sh` of short capitalised keys, and a `manifest.tsv` or a
`manifest.yml` that carries its settings as well as its record, since they have nowhere
readable to be but there. `simfbd fit`, `summarize`, `ranks` and `status` still work on it: the shell
driver it was made with is carried in `scripts/` and is called for it and for nothing else.
It cannot be regenerated, and one with a `manifest.tsv` cannot be extended. A new batch takes
a `config.yml`, and `simfbd config convert old/config.sh` writes one from an old config,
saying what it left out.

## The two ceilings on a simulation

Batches that differ only in the origin prior pool into a range of tree sizes. A deep origin
wants `simulation_limits.max_lineages`, a ceiling on lineages ever born, because it leaves
the oldest interval unbounded and a draw in the upper tail of the speciation rate then never
finishes simulating.

`simulation_limits.max_waiting_time_draws` is the second ceiling, on the work one simulation
may do. It counts the waiting times the birth-death simulator draws, across every tree it
tries, and refuses the draw when the count is passed. The count follows from the random
stream, so the same seeds refuse the same draws on any machine and a batch reproduces from
`seeds.RData`. A tree of n lineages takes about 3n draws, so keep it several times
`max_lineages`: 100000 cuts a runaway off after roughly 20 seconds on the machine it was
measured on, 200000 after roughly 60. The simulator prints how many draws were refused and
the most any finished simulation used.

The first ceiling is the one filter that costs something. It is applied by refusing the draw
before it is simulated, so it is a function of theta alone: it truncates the prior the
simulator draws from while the analysis keeps the untruncated one, and the SBC ranks are
non-uniform by exactly that much. A floor and ceiling on the number of reported taxa,
`accepted_replicates`, is a function of the record and is free, because conditioning on a
function of the record leaves the posterior given that record alone. Every batch records
`prior_filter`, the fraction of draws the ceiling refused, in its manifest.

How much it costs is a question of resolution: a truncation of fraction e shifts a mean rank
by roughly e/2, against a standard error of 0.289/sqrt(n). At n = 10000 that is 0.0029, so e
has to stay near or below 0.5% for the truncation to be invisible. Whether a given prior can
meet that is a question for the study that chooses it; `prior_filter` is what lets it be
measured rather than assumed.

## The two key axes

- **`inference.record_read`** = `complete` | `firstlast`: which record the analysis reads. All
  come from the same simulated history.
- **`inference.record_treated_as`** = `complete` | `incomplete`: how the model treats it.

Matching them is validation; mismatching them is a misspecification analysis.

## Dating resolution

A batch can report its record on other occurrence grids beside its own, so that analyses differing
only in how finely the fossils are dated are paired on the same fossils:

```yaml
record:
  occurrence_bins:
    width: 1.5
    oldest_edge: 36
    also_widths: [0, 0.75, 3, 6]     # 0 is the exact ages
    also_offsets: [0, 0.25]          # fractions of a bin
```

Each width at each offset is written as `specimens_<model>_w<width>`, with `_o<offset>` when
the offset is not 0, and a run reads one with `record_read: complete_w3_o0.25`. These grids
stop at the present and at `origin_age.max`, so an offset grid has no bin past the origin
prior. The batch's own grid is written as it always was.

How an analysis treats the ages it reads is the run's choice, `inference.occurrence_ages`:
`reported` hands rb the bins, `midpoint` puts each occurrence at its bin's middle, and `drawn`
at a uniform draw within it. `inference.age_draws: 10` makes ten such draws of each replicate
and runs a chain on each, which is PyRate's practice, and pools them: `output/rep_<n>.log` is
then the length of one chain, and after burn-in takes each sample from the next draw in turn,
so each draw counts the same. The chains themselves are under `aux/draw_<k>/`. A midpoint
run is refused where a midpoint is a rate-shift time, since rb cannot start there.

## Why the parameters are tracked this way

The hyperpriors have to reach two places: the simulator that draws the true values, and
the rb analysis that puts priors on them. When they are written down twice they drift,
and the drift is silent because both halves still run. That has already happened twice
here: once when the origin prior was `U(9,15)` in the simulator and `U(5,15)` on record,
and once when a copied analysis read one batch's data while being scored against another
batch's truth.

So the batch's `config.yml` is the single definition. simfbd hands it to the simulator and
writes the same values into the per-replicate `aux/run_<rep>.Rev` the Rev template reads.
Neither script contains a prior of its own. Four rules keep it that way:

- **A missing setting is fatal.** The simulator has no defaults for what a batch must say,
  because a default is what lets a renamed or forgotten setting produce a plausible run.
- **A run cannot set what the batch owns.** A batch setting in a run's config is refused,
  so there is never a second definition to disagree with.
- **The simulator reports what it drew with.** `sim_params.tsv` is compared against the
  config, so a setting that fails to reach it is caught rather than assumed.
- **An edited config is refused, not mixed.** Each manifest records the hash of the config
  that produced its data, and a run records the batch hash it ran against. Changing either
  after the fact stops the run instead of quietly mixing parameters across replicates.
  `simfbd fit --force` overrides for a run.

A config and its manifest together make an output self-describing: the settings, the config
hashes, the script hashes and the rb binary. `summarize.R` reads them rather than re-deriving
anything from directory names, so a result cannot be scored against the wrong truth.

The names are in `settings.py`, one line a setting: the name it is written under, what it
means, and the name the simulator and the Rev templates know it by. That last is how a value
reaches them and nothing a person writes or reads. A new setting gets a line there.

## What is not reproducible

A simulation reproduces: the same seeds give the same batch, byte for byte. An analysis does too,
on one machine with one binary: each replicate's chain is seeded with `mcmc.seed` plus the
replicate's number, so running it again gives the same log. Across machines or builds of rb a
chain may differ, and its summaries agree to within Monte Carlo error.

A run under a batch made before configs were YAML is analyzed by the shell driver, which sets no
seed, so those analyses differ from one time to the next.

## Scope

`inference.model: fbdr`, `bds` and `fbdsp` are wired. PyRate is not: its scripts are in
`src/simfbd/scripts/pyrate/`, and `pyrate.sh` beside them reads the same
`taxon/min_age/max_age/status` specimen format the simulator already writes, so reconnecting
it needs no separate simulator. The previous TreeSim pipeline is in git history at 8563117.
