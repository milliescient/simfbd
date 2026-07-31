#!/bin/bash
# PyRate inference for one replicate. Reads the SAME model-agnostic simulator
# output the rb methods use (DATADIR/specimens/taxa_<rep>.tsv), converts it to
# PyRate input format (the only PyRate-specific step — it lives here, not in
# sim.r), runs PyRate over the same interval timeline, writes OUTDIR/pyrate_<rep>.log.
#
# REQUIRES PyRate installed; set PYRATE_HOME (default /usr/local/bin/PyRate).
# Not testable on this machine (PyRate not installed) — wiring/structure only.
#
# args: <rep> <DATADIR> <OUTDIR> <NINTERVALS> <INTERVAL_WIDTH>
set -u
rep="$1"; DATADIR="$2"; OUTDIR="$3"; NINTERVALS="$4"; IW="$5"
PYRATE_HOME="${PYRATE_HOME:-/usr/local/bin/PyRate}"
RSCRIPT="${RSCRIPT:-Rscript}"
taxa="$DATADIR/specimens/taxa_${rep}.tsv"
work="$OUTDIR/pyrate_work_$rep"; mkdir -p "$work"

if [ ! -x "$PYRATE_HOME/PyRate.py" ]; then
  echo "PyRate not found at \$PYRATE_HOME=$PYRATE_HOME — skipping rep $rep"; exit 1
fi

# interval breakpoints (interior) — same timeline the rb analysis uses
: > "$work/epochs.txt"
for ((j=1; j<NINTERVALS; j++)); do echo "$((j*IW))" >> "$work/epochs.txt"; done

# convert the shared taxa file -> PyRate input (writes <taxa>_PyRate.py beside it)
"$RSCRIPT" -e "source('$PYRATE_HOME/pyrate_utilities.r'); extract.ages(file='$taxa', random=FALSE, replicates=1)" \
  > "$work/extract.out" 2>&1

# run PyRate (fixed-shift), same style as the original comparison
"$PYRATE_HOME/PyRate.py" "${taxa%.tsv}_PyRate.py" -fixShift "$work/epochs.txt" \
  -cauchy 0 0 -n 1000000 -s 100 -p 1000000 -out "_$rep" -wd "$work" > "$work/pyrate.out" 2>&1

mv "$work"/pyrate_mcmc_logs/*_BDS_mcmc.log "$OUTDIR/pyrate_${rep}.log" 2>/dev/null
