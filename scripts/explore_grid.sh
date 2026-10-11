#!/bin/bash
# Exploratory sweep for the final SBC design. Simulation only, no inference.
# Varies interval count, sampling rate and speciation rate; reports how often the
# classical estimators are defined and how big the records get.
set -u
SIMS="${SIMS:-sims}"      # where the batches go
RB=/research/phyloworks/rb-estimators/projects/cmake/build-estimators/rb
NREPS=${NREPS:-40}
for nint in 5 10; do
  if [ "$nint" = 5 ]; then w=3; else w=1.5; fi
  for pmean in -1.5 -0.5; do
    for lmean in -1.5 -0.9; do
      name="exp-i${nint}-p${pmean#-}-l${lmean#-}"
      rm -rf "$SIMS/$name"; mkdir -p "$SIMS/$name"
      cat > "$SIMS/$name/config.yml" <<CFG
timeline:
  rate_intervals: $nint
  interval_width: $w
origin_age:
  min: 15
  max: 21
  origin_lineage_must_be_sampled: false
rate_priors:
  speciation: {log_mean: $lmean, log_sd: 0.5}
  extinction: {log_mean: -1.5, log_sd: 0.5}
  fossil_sampling: {log_mean: $pmean, log_sd: 0.5}
  anagenetic_speciation: 0
rates_through_time:
  random_walk:
    step_sd: 0.5
record:
  occurrence_bins: {width: 1.5, oldest_edge: 36}
accepted_replicates:
  min_taxa: 10
  max_taxa: 400
  taxa_counted: record
  min_fossil_taxa: 5
simulation_limits:
  max_lineages: 400
  max_waiting_time_draws: 100000
CFG
      echo "### $name"
      simfbd extend "$SIMS/$name" "$NREPS" > "$SIMS/$name/sim.log" 2>&1 \
        && awk 'END{}' /dev/null \
        || { echo "  simulation failed, see sims/$name/sim.log"; continue; }
      python3 - "$SIMS/$name" <<'SZ'
import glob, sys, os
fs = glob.glob(os.path.join(sys.argv[1], "specimens_complete", "*.tsv"))
tx = []; oc = []
for f in fs:
    rows = open(f).read().strip().split("\n")[1:]
    tx.append(len({r.split("\t")[0] for r in rows})); oc.append(len(rows))
if tx: print(f"  median {sorted(tx)[len(tx)//2]} taxa, {sorted(oc)[len(oc)//2]} occurrences")
SZ
      simfbd script trad_coverage.py "$SIMS/$name" "$RB" complete | tail -5
    done
  done
done
