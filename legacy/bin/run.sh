#!/bin/bash
# Run one analysis: bash bin/run.sh sims/<name>
# Reads sims/<name>/config.sh. ONE model-agnostic simulator (bin/sim.r) generates
# the data; the inference depends on MODEL:
#   MODEL=fbdr   -> bin/infer.Rev (dnFBDRMatrix, BDS=false)
#   MODEL=bds    -> bin/infer.Rev (dnFBDRMatrix, BDS=true)
#   MODEL=pyrate -> bin/pyrate.sh (PyRate; converts the data to PyRate format)
#   SIM   = complete | incomplete   (how data is GENERATED/recorded)
#   INFER = complete | incomplete   (how the rb likelihood treats it)
# SIM==INFER => validation; SIM!=INFER => misspecification analysis.
set -u
ANA="${1:?usage: run.sh sims/<name>}"; ANA="${ANA%/}"
BIN="$(cd "$(dirname "$0")" && pwd)"
source "$BIN/env.sh"
source "$ANA/config.sh"
DATADIR="$ANA/data"; OUTDIR="$ANA/output"
mkdir -p "$DATADIR" "$OUTDIR" "$ANA/aux" "$ANA/results"
rm -f "$DATADIR"/truevals_*.tsv

export SIM NINTERVALS INTERVAL_WIDTH RHO LMEAN LSD PMEAN PSD AGE_MIN AGE_MAX NMIN DATADIR
case "$INFER" in complete) COMPLETE=true;; incomplete) COMPLETE=false;; *) echo "INFER must be complete|incomplete"; exit 1;; esac
case "$MODEL" in
  fbdr)   METHOD=rb; BDS=false ;;
  bds)    METHOD=rb; BDS=true  ;;
  pyrate) METHOD=pyrate ;;
  *) echo "MODEL must be fbdr|bds|pyrate"; exit 1 ;;
esac
[ "$NINTERVALS" -gt 1 ] && SKY=true || SKY=false     # NINTERVALS>1 => skyline
echo "sim=$(basename "$ANA"): MODEL=$MODEL SIM=$SIM INFER=$INFER(complete=$COMPLETE) skyline=$SKY cond=$COND rho=$RHO reps=$NREPS"
[ "$SIM" != "$INFER" ] && echo "  >> MISSPECIFICATION analysis (SIM != INFER)"

run_one() {
  local rep="$1"
  "$RSCRIPT" "$BIN/sim.r" "$rep" "$rep" >/dev/null 2>"$ANA/aux/sim_$rep.err" || return   # shared, model-agnostic
  if [ "$METHOD" = "pyrate" ]; then
    "$BIN/pyrate.sh" "$rep" "$DATADIR" "$OUTDIR" "$NINTERVALS" "$INTERVAL_WIDTH" > "$ANA/aux/pyrate_$rep.out" 2>&1
  else
    # Single self-contained file: set the config, then source the shared template.
    # (RevBayes' CLI sources only the FIRST file arg and treats the rest as args[],
    # so the old `rb aux.Rev infer.Rev` two-file form no longer runs infer.Rev.)
    cat > "$ANA/aux/run_$rep.Rev" <<EOF
rep <- "$rep"
SKYLINE <- $SKY
COMPLETE <- $COMPLETE
BDS <- $BDS
COND <- "$COND"
RHO <- $RHO
GENS <- $GENS
PRINTGEN <- $PRINTGEN
DATADIR <- "$DATADIR"
OUTDIR <- "$OUTDIR"
source("$BIN/infer.Rev")
EOF
    "$RBIN" "$ANA/aux/run_$rep.Rev" < /dev/null > "$ANA/aux/rb_$rep.out" 2>&1
  fi
}
export -f run_one
export BIN RSCRIPT RBIN ANA SKY COMPLETE BDS COND RHO GENS PRINTGEN DATADIR OUTDIR METHOD NINTERVALS INTERVAL_WIDTH
seq 1 "$NREPS" | xargs -P "$NCORES" -I {} bash -c 'run_one "$@"' _ {}

# assemble true values (youngest-first index => DIRECT mapping)
ni=$NINTERVALS
hdr="rep"; for p in lambda mu psi; do for i in $(seq 1 "$ni"); do hdr="$hdr\t$p$i"; done; done; hdr="$hdr\torigin\tnsamp"
{ echo -e "$hdr"; cat "$DATADIR"/truevals_*.tsv 2>/dev/null | sort -n; } > "$DATADIR/true_vals.tsv"
echo "logs: $(ls "$OUTDIR"/*_*.log 2>/dev/null | wc -l)/$NREPS"

# diagnostics: SBC for the rb methods (rate logs); PyRate logs are summarised by bin/summarize.R across analyses
if [ "$METHOD" = "rb" ]; then
  "$PY" "$BIN/sbc_diagnostics.py" "$DATADIR" "$OUTDIR" "${MAPPING:-direct}" | tee "$ANA/results/sbc_summary.txt"
  mv "$OUTDIR"/sbc_diagnostics.png "$ANA/results/" 2>/dev/null || true
else
  echo "PyRate logs in $OUTDIR; compare across analyses with bin/summarize.R" | tee "$ANA/results/sbc_summary.txt"
fi
echo "results in $ANA/results/"
