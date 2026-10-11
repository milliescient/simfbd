#!/bin/bash
# Does the taxon window change what the classical estimators can do, and at how many
# intervals? Simulation only. The psi sweep stays on throughout.
set -u
RB=/research/phyloworks/rb-estimators/projects/cmake/build-estimators/rb
NREPS=${NREPS:-100}
for nint in 6 10; do
  if [ "$nint" = 6 ]; then w=2.5; agemin=15; agemax=21
  else w=1.5; agemin=15; agemax=21; fi
  for win in "5 30" "5 50" "10 50" "10 100"; do
    set -- $win; lo=$1; hi=$2
    name="tx-i${nint}-${lo}to${hi}"
    rm -rf "sims/$name"; mkdir -p "sims/$name"
    cat > "sims/$name/config.sh" <<CFG
NINTERVALS=$nint
INTERVAL_WIDTH=$w
LMEAN=-1.5; LSD=0.4
MMEAN=-1.5; MSD=0.4
PMEAN=-1.2; PSD=0.85
AGE_MIN=$agemin; AGE_MAX=$agemax
BIN_WIDTH=1.5; BIN_MAX=36
LAMBDA_A=0
ORIGIN_SAMPLED=false
GMRF_SD=0.4
MIN_TAXA=$lo
MAX_TAXA=$hi
TAXA_COUNT=record
MIN_FOSSIL_TAXA=5
MAX_LINEAGES=400
SIM_MAX_DRAWS=100000
CFG
    if bash bin/extend.sh "sims/$name" "$NREPS" > "sims/$name/sim.log" 2>&1; then
      python3 bin/explore_report.py "sims/$name" "$RB"
    else
      echo "  $name: simulation failed"; tail -2 "sims/$name/sim.log"
    fi
  done
done
