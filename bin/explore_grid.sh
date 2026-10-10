#!/bin/bash
# Exploratory sweep for the final SBC design. Simulation only, no inference.
# Varies interval count, sampling rate and speciation rate; reports how often the
# classical estimators are defined and how big the records get.
set -u
RB=/research/phyloworks/rb-estimators/projects/cmake/build-estimators/rb
NREPS=${NREPS:-40}
for nint in 5 10; do
  if [ "$nint" = 5 ]; then w=3; else w=1.5; fi
  for pmean in -1.5 -0.5; do
    for lmean in -1.5 -0.9; do
      name="exp-i${nint}-p${pmean#-}-l${lmean#-}"
      rm -rf "sims/$name"; mkdir -p "sims/$name"
      cat > "sims/$name/config.sh" <<CFG
NINTERVALS=$nint
INTERVAL_WIDTH=$w
LMEAN=$lmean; LSD=0.5
MMEAN=-1.5; MSD=0.5
PMEAN=$pmean; PSD=0.5
AGE_MIN=15; AGE_MAX=21
BIN_WIDTH=1.5; BIN_MAX=36
LAMBDA_A=0
ORIGIN_SAMPLED=false
GMRF_SD=0.5
NREPS=$NREPS
MIN_TAXA=10
MAX_TAXA=400
TAXA_COUNT=record
MIN_FOSSIL_TAXA=5
MAX_LINEAGES=400
SIM_TIMEOUT=20
CFG
      echo "### $name"
      bash bin/run.sh "sims/$name" > "sims/$name/sim.log" 2>&1 \
        && awk 'END{}' /dev/null \
        || { echo "  simulation failed, see sims/$name/sim.log"; continue; }
      python3 - "sims/$name" <<'SZ'
import glob, sys, os
fs = glob.glob(os.path.join(sys.argv[1], "specimens_complete", "*.tsv"))
tx = []; oc = []
for f in fs:
    rows = open(f).read().strip().split("\n")[1:]
    tx.append(len({r.split("\t")[0] for r in rows})); oc.append(len(rows))
if tx: print(f"  median {sorted(tx)[len(tx)//2]} taxa, {sorted(oc)[len(oc)//2]} occurrences")
SZ
      python3 bin/trad_coverage.py "sims/$name" "$RB" complete | tail -5
    done
  done
done
