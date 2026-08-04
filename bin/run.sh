#!/bin/bash
# Run one batch or one analysis.
#
#   bash bin/run.sh sims/<batch>                  # simulate the batch only
#   bash bin/run.sh sims/<batch>/runs/<name>      # simulate if needed, then infer
#
# A batch owns the generating process and every hyperprior. One simulation pass emits
# all reporting models from one record, so the analyses under it share a timeline and
# one true_vals.tsv and their comparisons are paired.
#
# An analysis under runs/ chooses only what to read and how to fit it. It inherits the
# batch's hyperpriors and may not redeclare them, and bin/run.sh feeds the same values
# to bin/sim.R and to bin/infer.Rev, so the generating and fitting priors cannot drift.
#
# Both levels write a manifest.tsv recording the resolved parameters, the config hash,
# the script hashes and the rb identity, and a run records the batch hash it ran against.
# A config edited after its data or output exists is refused rather than mixed; FORCE=1
# overrides.
set -u

TARGET="${1:?usage: run.sh sims/<batch>[/runs/<name>]}"; TARGET="${TARGET%/}"
BIN="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$BIN/.." && pwd)"
cd "$ROOT" || exit 1
source "$BIN/env.sh"

# the hyperpriors a batch owns; an analysis that sets one of these is an error
PRIOR_KEYS="NINTERVALS INTERVAL_WIDTH LMEAN LSD MMEAN MSD PMEAN PSD AGE_MIN AGE_MAX BIN_WIDTH BIN_MAX NREPS LAMBDA_A ORIGIN_SAMPLED"
# owned by the batch too, but optional, so configs written before it keep working
OPT_KEYS="MAX_LINEAGES SIM_TIMEOUT RHO MIN_TAXA MAX_TAXA GMRF_SD PRESENT"

hash_of() { sha1sum "$1" | cut -c1-12; }

# A content hash says whether a script still matches, but cannot restore it. The commit can,
# so record both: --dirty marks the case where the commit does not describe what ran.
commit_of() { git -C "$1" describe --always --dirty --tags 2>/dev/null || echo none; }

# bin/ is tracked, so a batch generated from a modified tree records a commit that will not
# reproduce it. Warn once rather than refuse, since a run mid-experiment is still worth having.
warn_if_dirty() {
  git -C "$ROOT" diff --quiet HEAD -- "$BIN" 2>/dev/null && return 0
  echo "warning: $BIN has uncommitted changes, so simfbd_commit will not reproduce this" >&2
}

die() { echo "error: $*" >&2; exit 1; }

# resolve batch vs analysis from the path shape
case "$TARGET" in
  */runs/*) RUNDIR="$TARGET"; BATCHDIR="${TARGET%%/runs/*}" ;;
  *)        RUNDIR=""; BATCHDIR="$TARGET" ;;
esac
[ -f "$BATCHDIR/config.sh" ] || die "no $BATCHDIR/config.sh"

# ---- batch config: the single definition of the generating process ----
set -a; source "$BATCHDIR/config.sh"; set +a
for k in $PRIOR_KEYS; do
  [ -n "${!k:-}" ] || die "$BATCHDIR/config.sh does not set $k"
done
BATCH_RHO="${RHO:-}"          # the generating value, before a run config shadows it
BATCH_CFG_HASH="$(hash_of "$BATCHDIR/config.sh")"
export BATCHDIR

# ---- simulate the batch, unless its data already matches this config ----
BATCH_MANIFEST="$BATCHDIR/manifest.tsv"
have_batch=false
if [ -f "$BATCH_MANIFEST" ]; then
  prev="$(awk -F'\t' '$1=="config_hash"{print $2}' "$BATCH_MANIFEST")"
  if [ "$prev" = "$BATCH_CFG_HASH" ]; then
    have_batch=true
  elif [ "${FORCE:-0}" != "1" ]; then
    die "$BATCHDIR/config.sh changed since its data was generated ($prev -> $BATCH_CFG_HASH).
       Re-simulating would mix parameters across replicates. Delete the data or set FORCE=1."
  fi
fi

if [ "$have_batch" = false ]; then
  echo "simulating $BATCHDIR: ${NREPS} reps, ${NINTERVALS} intervals of width ${INTERVAL_WIDTH}"
  warn_if_dirty
  "$RSCRIPT" "$BIN/sim.R" || die "simulation failed"
  {
    printf 'kind\tbatch\n'
    printf 'batch\t%s\n' "$(basename "$BATCHDIR")"
    printf 'config_hash\t%s\n' "$BATCH_CFG_HASH"
    printf 'sim_script_hash\t%s\n' "$(hash_of "$BIN/sim.R")"
    printf 'simfbd_commit\t%s\n' "$(commit_of "$ROOT")"
    printf 'created\t%s\n' "$(date -Is)"
    # batches simulated before this key was added hold true_vals oldest-first
    printf 'true_vals_order\tyoungest_first\n'
    for k in $PRIOR_KEYS; do printf '%s\t%s\n' "$k" "${!k}"; done
    for k in $OPT_KEYS; do printf '%s\t%s\n' "$k" "${!k:-none}"; done
  } > "$BATCH_MANIFEST"

  # the simulator reports what it actually drew with; disagreement means a key did not reach it
  if [ -f "$BATCHDIR/sim_params.tsv" ]; then
    while IFS=$'\t' read -r k v; do
      want="${!k:-}"
      [ "$(awk -v a="$v" 'BEGIN{print a+0}')" = "$(awk -v a="$want" 'BEGIN{print a+0}')" ] \
        || die "simulator used $k=$v but config.sh says $want"
    done < "$BATCHDIR/sim_params.tsv"
  fi
  echo "  -> $BATCH_MANIFEST"
fi

[ -n "$RUNDIR" ] || { echo "batch ready: $BATCHDIR"; exit 0; }

# ---- analysis config: what to read and how to fit it ----
[ -f "$RUNDIR/config.sh" ] || die "no $RUNDIR/config.sh"
for k in $PRIOR_KEYS; do
  grep -Eq "^[[:space:]]*$k=" "$RUNDIR/config.sh" \
    && die "$RUNDIR/config.sh sets $k, which the batch owns. Remove it, or make a new batch."
done
for k in MAX_LINEAGES SIM_TIMEOUT; do
  grep -Eq "^[[:space:]]*$k=" "$RUNDIR/config.sh" \
    && die "$RUNDIR/config.sh sets $k, which the batch owns. Remove it, or make a new batch."
done
set -a; source "$RUNDIR/config.sh"; set +a
ORIGIN_PRIOR="${ORIGIN_PRIOR:-uniform}"
# rho is the generating parameter and the fitted one. A batch that sets it wins, and a run
# that disagrees is refused rather than silently fitting a rho the data was not drawn under.
# compared as numbers, so 0 and 0.0 are not a disagreement
if [ -n "$BATCH_RHO" ] && ! awk -v a="$RHO" -v b="$BATCH_RHO" 'BEGIN{exit !(a+0==b+0)}'; then
  die "$RUNDIR/config.sh fits RHO=$RHO but $BATCHDIR generated under RHO=$BATCH_RHO."
fi
SURVIVORS="${SURVIVORS:-true}"
for k in MODEL REPORTING INFER COND RHO GENS PRINTGEN NCORES; do
  [ -n "${!k:-}" ] || die "$RUNDIR/config.sh does not set $k"
done
case "$MODEL" in
  fbdr) TEMPLATE="$BIN/infer.Rev" ;;       # dnFBDRP: accounts for unsampled lineages
  bds)  TEMPLATE="$BIN/infer-bds.Rev" ;;   # dnBDS: complete lineage sampling, PyRate's assumption
  fbdsp) TEMPLATE="$BIN/infer-fbdsp.Rev" ;; # dnFBDSP: the tree process, the only one with lambda_a
  *)    die "MODEL=$MODEL is not wired; use fbdr, bds or fbdsp" ;;
esac
case "$INFER" in
  complete)   COMPLETE=true ;;
  incomplete) COMPLETE=false ;;
  *) die "INFER must be complete|incomplete" ;;
esac
SPECIMENS="specimens_$REPORTING"
[ -d "$BATCHDIR/$SPECIMENS" ] || die "batch has no $SPECIMENS (reporting model $REPORTING)"
[ "$NINTERVALS" -gt 1 ] && SKY=true || SKY=false

RUN_CFG_HASH="$(hash_of "$RUNDIR/config.sh")"
RUN_MANIFEST="$RUNDIR/manifest.tsv"
if [ -f "$RUN_MANIFEST" ] && [ "${FORCE:-0}" != "1" ]; then
  pb="$(awk -F'\t' '$1=="batch_config_hash"{print $2}' "$RUN_MANIFEST")"
  pr="$(awk -F'\t' '$1=="config_hash"{print $2}' "$RUN_MANIFEST")"
  [ "$pb" = "$BATCH_CFG_HASH" ] || die "$RUNDIR ran against batch config $pb, now $BATCH_CFG_HASH. Its output would mix. Delete output/ or set FORCE=1."
  [ "$pr" = "$RUN_CFG_HASH" ]   || die "$RUNDIR/config.sh changed since its output was written ($pr -> $RUN_CFG_HASH). Delete output/ or set FORCE=1."
fi

OUTDIR="$RUNDIR/output"; AUXDIR="$RUNDIR/aux"
mkdir -p "$OUTDIR" "$AUXDIR" "$RUNDIR/results"
echo "$(basename "$RUNDIR"): batch=$(basename "$BATCHDIR") reporting=$REPORTING infer=$INFER(complete=$COMPLETE) cond=$COND rho=$RHO reps=$NREPS"
warn_if_dirty

# a finished rep's log has this many lines; reruns fill gaps rather than redo work
COMPLETE_LINES=$(( GENS / PRINTGEN + 2 ))
: > "$RUNDIR/failures.log"

run_one() {
  local rep="$1"
  local log="$OUTDIR/rep_$rep.log"
  if [ -f "$log" ] && [ "$(wc -l < "$log")" = "$COMPLETE_LINES" ]; then return 0; fi
  # one self-contained file per rep: the config as Rev variables, then the shared template
  # the oldest bin floor across taxa: PyRate's max FA, and the tightest lower bound the
  # data puts on the origin. Only used when ORIGIN_PRIOR=exponential.
  local maxfa
  maxfa=$(awk -F'\t' 'NR>1 { if ($2+0 > m[$1]) m[$1]=$2+0 }
                       END { x=0; for (t in m) if (m[t]>x) x=m[t]; printf "%.6f", (x>0 ? x : 1) }' \
          "$BATCHDIR/$SPECIMENS/taxa_${rep}.tsv")
  cat > "$AUXDIR/run_$rep.Rev" <<EOF
rep <- "$rep"
ORIGIN_PRIOR <- "$ORIGIN_PRIOR"
MAXFA <- $maxfa
SKYLINE <- $SKY
COMPLETE <- $COMPLETE
COND <- "$COND"
RHO <- $RHO
GENS <- $GENS
PRINTGEN <- $PRINTGEN
BATCHDIR <- "$BATCHDIR"
SPECIMENS <- "$SPECIMENS"
OUTDIR <- "$OUTDIR"
LMEAN <- $LMEAN
LSD <- $LSD
MMEAN <- $MMEAN
MSD <- $MSD
PMEAN <- $PMEAN
PSD <- $PSD
AGE_MIN <- $AGE_MIN
AGE_MAX <- $AGE_MAX
SURVIVORS <- $SURVIVORS
GMRF_SD <- ${GMRF_SD:-0}
LAMBDA_A <- ${LAMBDA_A:-0}
PRESENT <- ${PRESENT:-0}
source("$TEMPLATE")
EOF
  "$RBIN" "$AUXDIR/run_$rep.Rev" < /dev/null > "$AUXDIR/rb_$rep.out" 2>&1
  if [ ! -f "$log" ] || [ "$(wc -l < "$log")" != "$COMPLETE_LINES" ]; then
    echo "$rep :: $(grep -m1 -i error "$AUXDIR/rb_$rep.out" || echo 'no log written')" >> "$RUNDIR/failures.log"
  fi
  return 0
}
export -f run_one
export ORIGIN_PRIOR
export RBIN BIN TEMPLATE BATCHDIR SPECIMENS OUTDIR AUXDIR RUNDIR SKY COMPLETE COND RHO GENS PRINTGEN \
       COMPLETE_LINES LMEAN LSD MMEAN MSD PMEAN PSD AGE_MIN AGE_MAX SURVIVORS GMRF_SD

seq 1 "$NREPS" | xargs -P "$NCORES" -I {} bash -c 'run_one "$@"' _ {}

{
  printf 'kind\trun\n'
  printf 'run\t%s\n' "$(basename "$RUNDIR")"
  printf 'batch\t%s\n' "$(basename "$BATCHDIR")"
  printf 'batch_dir\t%s\n' "$BATCHDIR"
  printf 'batch_config_hash\t%s\n' "$BATCH_CFG_HASH"
  printf 'config_hash\t%s\n' "$RUN_CFG_HASH"
  printf 'infer_script_hash\t%s\n' "$(hash_of "$TEMPLATE")"
  printf 'survivors\t%s\n' "$SURVIVORS"
  printf 'simfbd_commit\t%s\n' "$(commit_of "$ROOT")"
  rb_bin="$(command -v "$RBIN" || echo "$RBIN")"
  printf 'rb_path\t%s\n' "$rb_bin"
  printf 'rb_md5\t%s\n' "$(md5sum "$rb_bin" 2>/dev/null | cut -c1-12 || echo unknown)"
  # the binary carries its own git describe, so this names what actually ran even after the
  # source tree moves on. An md5 identifies a build, not a commit: rb embeds its build date.
  printf 'rb_commit\t%s\n' "$(strings "$rb_bin" 2>/dev/null | grep -m1 -E '^[A-Za-z0-9._-]+-[0-9]+-g[0-9a-f]{6,}$' || echo unknown)"
  printf 'created\t%s\n' "$(date -Is)"
  printf 'specimens\t%s\n' "$SPECIMENS"
  for k in MODEL REPORTING INFER COND RHO GENS PRINTGEN ORIGIN_PRIOR; do printf '%s\t%s\n' "$k" "${!k}"; done
  # the inherited hyperpriors, so a run manifest describes its own fit without the batch
  for k in $PRIOR_KEYS; do printf '%s\t%s\n' "$k" "${!k}"; done
  for k in $OPT_KEYS; do printf '%s\t%s\n' "$k" "${!k:-none}"; done
} > "$RUN_MANIFEST"

nfail=$(wc -l < "$RUNDIR/failures.log")
echo "logs: $(ls "$OUTDIR"/rep_*.log 2>/dev/null | wc -l)/$NREPS, $nfail failed"
echo "manifest: $RUN_MANIFEST"
echo "now run: $RSCRIPT bin/summarize.R $RUNDIR"
