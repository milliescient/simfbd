#!/bin/bash
# PyRate inference for one replicate. Reads the SAME model-agnostic simulator
# output the rb methods use (DATADIR/specimens/taxa_<rep>.tsv), converts it to
# PyRate input format (the only PyRate-specific step — it lives here, not in
# sim.r), runs PyRate over the same interval timeline, writes OUTDIR/pyrate_<rep>.log.
#
# PYRATE_HOME defaults to the checkout under /research/phyloworks/tools, which has its
# own venv: PyRate 3.1.3 pins numpy < 2 and the shared one is on 2.x.
#
# args: <rep> <DATADIR> <OUTDIR> <NINTERVALS> <INTERVAL_WIDTH> [SPECIMENS]
set -u
rep="$1"; DATADIR="$2"; OUTDIR="$3"; NINTERVALS="$4"; IW="$5"; SPECIMENS="${6:-specimens_complete}"
PYRATE_HOME="${PYRATE_HOME:-/research/phyloworks/tools/PyRate}"
PYRATE_PY="${PYRATE_PY:-/research/phyloworks/tools/pyrate-venv/bin/python}"
RSCRIPT="${RSCRIPT:-/research/phyloworks/mm/root/envs/fsim/bin/Rscript}"
taxa="$DATADIR/$SPECIMENS/taxa_${rep}.tsv"
[ -f "$taxa" ] || { echo "no $taxa"; exit 1; }
work="$OUTDIR/pyrate_work_$rep"; mkdir -p "$work"

[ -f "$PYRATE_HOME/PyRate.py" ] || { echo "PyRate not found at \$PYRATE_HOME=$PYRATE_HOME"; exit 1; }
[ -x "$PYRATE_PY" ]                 || { echo "no python at \$PYRATE_PY=$PYRATE_PY"; exit 1; }

# interval breakpoints (interior) — same timeline the rb analysis uses
: > "$work/epochs.txt"
for ((j=1; j<NINTERVALS; j++)); do echo "$((j*IW))" >> "$work/epochs.txt"; done

# convert the shared taxa file -> PyRate input. extract.ages wants Species/Status/MinT/MaxT,
# not the simulator's taxon/status/min_age/max_age, so rename before handing it over.
"$RSCRIPT" -e "
  source('$PYRATE_HOME/pyrate_utilities.r')
  d <- read.table('$taxa', header=TRUE, sep='\t')
  write.table(data.frame(Species=d\$taxon, Status=d\$status, MinT=d\$min_age, MaxT=d\$max_age),
              '$work/occ.txt', sep='\t', row.names=FALSE, quote=FALSE)
  extract.ages(file='$work/occ.txt', random=FALSE, replicates=1)
" > "$work/extract.out" 2>&1
[ -f "$work/occ_PyRate.py" ] || { echo "extract.ages failed; see $work/extract.out"; exit 1; }

# BDS with the same breakpoints for the birth-death rates and for preservation, and a
# homogeneous Poisson process within each interval -- the model dnBDS implements.
"$PYRATE_PY" "$PYRATE_HOME/PyRate.py" "$work/occ_PyRate.py" \
  -fixShift "$work/epochs.txt" -qShift "$work/epochs.txt" -mHPP -A 0 \
  -n "${PYRATE_GENS:-1000000}" -s 1000 -p 1000000 -b 0 \
  -out "_$rep" -wd "$work" > "$work/pyrate.out" 2>&1

mv "$work"/pyrate_mcmc_logs/*_BDS_mcmc.log "$OUTDIR/pyrate_${rep}.log" 2>/dev/null \
  || { echo "PyRate wrote no mcmc log; see $work/pyrate.out"; exit 1; }
