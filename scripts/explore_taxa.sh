#!/bin/bash
# Does the taxon window change what the classical estimators can do, and at how many
# intervals? Simulation only. The psi sweep stays on throughout.
set -u
SIMS="${SIMS:-sims}"      # where the batches go
RB=/research/phyloworks/rb-estimators/projects/cmake/build-estimators/rb
NREPS=${NREPS:-100}
for nint in 6 10; do
  if [ "$nint" = 6 ]; then w=2.5; agemin=15; agemax=21
  else w=1.5; agemin=15; agemax=21; fi
  for win in "5 30" "5 50" "10 50" "10 100"; do
    set -- $win; lo=$1; hi=$2
    name="tx-i${nint}-${lo}to${hi}"
    rm -rf "$SIMS/$name"; mkdir -p "$SIMS/$name"
    cat > "$SIMS/$name/config.yml" <<CFG
timeline:
  rate_intervals: $nint
  interval_width: $w
origin_age:
  min: $agemin
  max: $agemax
  origin_lineage_must_be_sampled: false
rate_priors:
  speciation: {log_mean: -1.5, log_sd: 0.4}
  extinction: {log_mean: -1.5, log_sd: 0.4}
  fossil_sampling: {log_mean: -1.2, log_sd: 0.85}
  anagenetic_speciation: 0
rates_through_time:
  random_walk:
    step_sd: 0.4
record:
  occurrence_bins: {width: 1.5, oldest_edge: 36}
accepted_replicates:
  min_taxa: $lo
  max_taxa: $hi
  taxa_counted: record
  min_fossil_taxa: 5
simulation_limits:
  max_lineages: 400
  max_waiting_time_draws: 100000
CFG
    if simfbd extend "$SIMS/$name" "$NREPS" > "$SIMS/$name/sim.log" 2>&1; then
      simfbd script explore_report.py "$SIMS/$name" "$RB"
    else
      echo "  $name: simulation failed"; tail -2 "$SIMS/$name/sim.log"
    fi
  done
done
