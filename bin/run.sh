#!/bin/bash
# Run one batch or one analysis.
#
#   bash bin/run.sh <batch>                  # simulate the batch only
#   bash bin/run.sh <batch>/runs/<name>      # simulate if needed, then infer
#
# bin/simfbd is the front door: `simfbd fit <batch>/runs/<name>` is this script.
#
# A batch owns the generating process and every hyperprior. One simulation pass emits
# all reporting models from one record, so the analyses under it share a timeline and
# one true_vals.tsv and their comparisons are paired.
#
# An analysis under runs/ chooses only what to read and how to fit it. It inherits the
# batch's hyperpriors and may not redeclare them, and bin/run.sh feeds the same values
# to bin/sim.R and to bin/infer.Rev, so the generating and fitting priors cannot drift.
#
# Both levels write a manifest.yml recording the resolved parameters, the config hash,
# the script hashes and the rb identity, and a run records the batch hash it ran against.
# A config edited after its data or output exists is refused rather than mixed.
#
# A batch does not say how many replicates it has: bin/extend.sh adds them, to a new batch
# or to one that already has some, and the manifest counts them. A batch simulated before
# that carries NREPS in its config.sh and a manifest.tsv, and is run as it always was, but
# cannot be extended.
set -u

TARGET="${1:?usage: run.sh <batch>[/runs/<name>]}"; TARGET="${TARGET%/}"
BIN="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$BIN/.." && pwd)"
# A batch is a directory anywhere. A relative path is looked for where the command was typed
# and then under the repository, which is where sims/<batch> has always been found; a batch
# inside the repository is still named relative to it, so what is printed and recorded for one
# does not depend on where the command was typed.
case "$TARGET" in
  /*) ;;
  *) if [ -f "$PWD/${TARGET%%/runs/*}/config.sh" ]; then TARGET="$PWD/$TARGET"; else TARGET="$ROOT/$TARGET"; fi ;;
esac
TARGET="$(realpath -m "$TARGET")"
case "$TARGET" in "$ROOT"/*) TARGET="${TARGET#"$ROOT"/}" ;; esac
cd "$ROOT" || exit 1
source "$BIN/env.sh"

# the hyperpriors a batch owns; an analysis that sets one of these is an error
PRIOR_KEYS="NINTERVALS INTERVAL_WIDTH LMEAN LSD MMEAN MSD PMEAN PSD AGE_MIN AGE_MAX BIN_WIDTH BIN_MAX NREPS LAMBDA_A ORIGIN_SAMPLED"
# owned by the batch too, but optional, so configs written before it keep working
OPT_KEYS="MAX_LINEAGES SIM_TIMEOUT SIM_MAX_DRAWS RHO MIN_TAXA MAX_TAXA GMRF_SD GMRF_SD_MEAN GMRF_SD_SD GMRF_PSI_MEAN GMRF_PSI_SD OU_THETA OU_DREF OU_PHI_A OU_PHI_B TMEAN TSD LANCH TANCH PANCH OU_SD_MEAN OU_SD_SD OU_PSI_MEAN OU_PSI_SD PRESENT REPORT_PROB TAXA_COUNT MIN_FOSSIL_TAXA"

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

# one recorded value, from a manifest.yml or a manifest.tsv alike
mf() { "$PY" "$BIN/manifest.py" get "$1" "$2"; }

# resolve batch vs analysis from the path shape
case "$TARGET" in
  */runs/*) RUNDIR="$TARGET"; BATCHDIR="${TARGET%%/runs/*}" ;;
  *)        RUNDIR=""; BATCHDIR="$TARGET" ;;
esac
[ -f "$BATCHDIR/config.sh" ] || die "no $BATCHDIR/config.sh"

# ---- batch config: the single definition of the generating process ----
BATCH_MANIFEST="$BATCHDIR/manifest.yml"; LEGACY=false
if [ -f "$BATCHDIR/manifest.tsv" ] && [ ! -f "$BATCH_MANIFEST" ]; then
  BATCH_MANIFEST="$BATCHDIR/manifest.tsv"; LEGACY=true
fi
EXTEND="${EXTEND:-0}"         # replicates to add, set by bin/extend.sh
case "$EXTEND" in ''|*[!0-9]*) die "EXTEND must be a whole number, not '$EXTEND'" ;; esac
if [ "$LEGACY" = false ]; then
  unset NREPS
  grep -Eq "(^|;)[[:space:]]*NREPS=" "$BATCHDIR/config.sh" \
    && die "$BATCHDIR/config.sh sets NREPS, which a batch no longer declares.
       Remove it, then: bash bin/extend.sh $BATCHDIR <how many>"
  grep -Eq "(^|;)[[:space:]]*SIM_TIMEOUT=" "$BATCHDIR/config.sh" \
    && die "$BATCHDIR/config.sh sets SIM_TIMEOUT, a wall-clock limit that keeps a batch from
       reproducing from its seeds. Set SIM_MAX_DRAWS instead; see the README for a value."
fi
set -a; source "$BATCHDIR/config.sh"; set +a
BATCH_CFG_HASH="$(hash_of "$BATCHDIR/config.sh")"
SIM_HASH="$(hash_of "$BIN/sim.R")"
export BATCHDIR

# ---- simulate the batch, or add to it, unless its data already matches this config ----
have_batch=false
if [ -f "$BATCH_MANIFEST" ]; then
  prev="$(mf "$BATCHDIR" config_hash)"
  if [ "$prev" = "$BATCH_CFG_HASH" ]; then
    have_batch=true
  elif [ "$LEGACY" = false ]; then
    die "$BATCHDIR/config.sh changed since its replicates were simulated ($prev -> $BATCH_CFG_HASH).
       Adding to them or redrawing them would mix parameters. Make a new batch."
  elif [ "${FORCE:-0}" != "1" ]; then
    die "$BATCHDIR/config.sh changed since its data was generated ($prev -> $BATCH_CFG_HASH).
       Re-simulating would mix parameters across replicates. Delete the data or set FORCE=1."
  fi
fi

REPS_HAVE=0
if [ "$LEGACY" = true ]; then
  [ "$EXTEND" = 0 ] || die "$BATCHDIR predates bin/extend.sh and cannot be extended: its
       replicate count is part of the config its runs were checked against. Make a new batch."
  simulate=true; [ "$have_batch" = true ] && simulate=false
else
  [ "$have_batch" = true ] && REPS_HAVE="$(mf "$BATCHDIR" NREPS)"
  NREPS=$(( REPS_HAVE + EXTEND ))
  [ "$NREPS" -gt 0 ] || die "$BATCHDIR has no replicates yet: bash bin/extend.sh $BATCHDIR <how many>"
  if [ "$EXTEND" -gt 0 ] && [ "$REPS_HAVE" -gt 0 ]; then
    was="$(mf "$BATCHDIR" sim_script_hash)"
    # same config, different simulator is still a mixed batch, and nothing else would say so
    [ "$was" = "$SIM_HASH" ] || [ "${ALLOW_SIM_CHANGE:-0}" = "1" ] \
      || die "bin/sim.R has changed since replicates up to $REPS_HAVE were simulated ($was -> $SIM_HASH).
       Make a new batch, or set ALLOW_SIM_CHANGE=1 to add to this one; the manifest records both."
  fi
  simulate=false; [ "$EXTEND" -gt 0 ] && simulate=true
fi
export NREPS REPS_HAVE
for k in $PRIOR_KEYS; do
  [ -n "${!k:-}" ] || die "$BATCHDIR/config.sh does not set $k"
done
BATCH_RHO="${RHO:-}"          # the generating value, before a run config shadows it

if [ "$simulate" = true ]; then
  echo "simulating $BATCHDIR: reps $(( REPS_HAVE + 1 )) to ${NREPS}, ${NINTERVALS} intervals of width ${INTERVAL_WIDTH}"
  warn_if_dirty
  "$RSCRIPT" "$BIN/sim.R" || die "simulation failed"
  if [ "$LEGACY" = false ]; then
    "$PY" "$BIN/manifest.py" write-batch "$BATCHDIR" first=$(( REPS_HAVE + 1 )) last="$NREPS" \
      config_hash="$BATCH_CFG_HASH" sim_script_hash="$SIM_HASH" \
      simfbd_commit="$(commit_of "$ROOT")" || die "could not write $BATCH_MANIFEST"
  else
  {
    printf 'kind\tbatch\n'
    printf 'batch\t%s\n' "$(basename "$BATCHDIR")"
    printf 'config_hash\t%s\n' "$BATCH_CFG_HASH"
    printf 'sim_script_hash\t%s\n' "$SIM_HASH"
    printf 'simfbd_commit\t%s\n' "$(commit_of "$ROOT")"
    printf 'created\t%s\n' "$(date -Is)"
    # batches simulated before this key was added hold true_vals oldest-first
    printf 'true_vals_order\tyoungest_first\n'
    for k in $PRIOR_KEYS; do printf '%s\t%s\n' "$k" "${!k}"; done
    for k in $OPT_KEYS; do printf '%s\t%s\n' "$k" "${!k:-none}"; done
    # what the theta screen actually refused, measured rather than configured
    [ -f "$BATCHDIR/screened.tsv" ] && cat "$BATCHDIR/screened.tsv"
  } > "$BATCH_MANIFEST"
  fi

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
for k in MAX_LINEAGES SIM_TIMEOUT SIM_MAX_DRAWS NREPS; do
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
# dnBDS only: leave out taxa seen only alive at the present, which PyRate's input never holds
DROP_PRESENT_ONLY="${DROP_PRESENT_ONLY:-false}"
for k in MODEL REPORTING INFER COND RHO GENS PRINTGEN NCORES; do
  [ -n "${!k:-}" ] || die "$RUNDIR/config.sh does not set $k"
done
case "$MODEL" in
  fbdr) TEMPLATE="$BIN/infer-fbdrp.Rev" ;;       # dnFBDRP: accounts for unsampled lineages
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
RUN_MANIFEST="$RUNDIR/manifest.yml"
[ "$LEGACY" = true ] && RUN_MANIFEST="$RUNDIR/manifest.tsv"
if [ -f "$RUN_MANIFEST" ] && [ "${FORCE:-0}" != "1" ]; then
  pb="$(mf "$RUNDIR" batch_config_hash)"
  pr="$(mf "$RUNDIR" config_hash)"
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
DROP_PRESENT_ONLY <- $DROP_PRESENT_ONLY
GMRF_SD <- ${GMRF_SD:-0}
GMRF_SD_MEAN <- ${GMRF_SD_MEAN:--1}
GMRF_SD_SD <- ${GMRF_SD_SD:-0}
GMRF_PSI_MEAN <- ${GMRF_PSI_MEAN:--1}
GMRF_PSI_SD <- ${GMRF_PSI_SD:-0}
BINDIR <- "$BIN"
OU_THETA <- ${OU_THETA:--1}
OU_DREF <- ${OU_DREF:-1.5}
OU_PHI_A <- ${OU_PHI_A:--1}
OU_PHI_B <- ${OU_PHI_B:--1}
TMEAN <- ${TMEAN:-0}
TSD <- ${TSD:-0.3}
LANCH <- ${LANCH:-0}
TANCH <- ${TANCH:-0}
PANCH <- ${PANCH:-0}
OU_SD_MEAN <- ${OU_SD_MEAN:--1}
OU_SD_SD <- ${OU_SD_SD:-0}
OU_PSI_MEAN <- ${OU_PSI_MEAN:--1}
OU_PSI_SD <- ${OU_PSI_SD:-0}
LAMBDA_A <- ${LAMBDA_A:-0}
PRESENT <- ${PRESENT:-0}
REPORT_PROB <- ${REPORT_PROB:-1}
source("$TEMPLATE")
EOF
  # SLOTS caps chains across every run.sh sharing SLOT_DIR, so concurrent runs fill the cores
  # between them: each chain holds one slot lock for as long as rb runs
  if [ -n "${SLOTS:-}" ]; then
    local fd s
    while true; do
      for s in $(seq 1 "$SLOTS"); do
        exec {fd}>"$SLOT_DIR/slot_$s"
        flock -n "$fd" && break 2
        exec {fd}>&-
      done
      sleep 10
    done
  fi
  "$RBIN" "$AUXDIR/run_$rep.Rev" < /dev/null > "$AUXDIR/rb_$rep.out" 2>&1
  [ -n "${SLOTS:-}" ] && exec {fd}>&-
  if [ ! -f "$log" ] || [ "$(wc -l < "$log")" != "$COMPLETE_LINES" ]; then
    echo "$rep :: $(grep -m1 -i error "$AUXDIR/rb_$rep.out" || echo 'no log written')" >> "$RUNDIR/failures.log"
  fi
  return 0
}
export -f run_one
export ORIGIN_PRIOR
export RBIN BIN TEMPLATE BATCHDIR SPECIMENS OUTDIR AUXDIR RUNDIR SKY COMPLETE COND RHO GENS PRINTGEN \
       COMPLETE_LINES LMEAN LSD MMEAN MSD PMEAN PSD AGE_MIN AGE_MAX SURVIVORS GMRF_SD \
       GMRF_SD_MEAN GMRF_SD_SD GMRF_PSI_MEAN GMRF_PSI_SD OU_THETA OU_DREF OU_PHI_A OU_PHI_B TMEAN TSD \
       LANCH TANCH PANCH OU_SD_MEAN OU_SD_SD OU_PSI_MEAN OU_PSI_SD DROP_PRESENT_ONLY \
       SLOTS SLOT_DIR

seq 1 "$NREPS" | xargs -P "$NCORES" -I {} bash -c 'run_one "$@"' _ {}

rb_bin="$(command -v "$RBIN" || echo "$RBIN")"
rb_md5="$(md5sum "$rb_bin" 2>/dev/null | cut -c1-12 || echo unknown)"
# the binary carries its own git describe, so this names what actually ran even after the
# source tree moves on. An md5 identifies a build, not a commit: rb embeds its build date.
rb_commit="$(strings "$rb_bin" 2>/dev/null | grep -m1 -E '^[A-Za-z0-9._-]+-[0-9]+-g[0-9a-f]{6,}$' || echo unknown)"
if [ "$LEGACY" = false ]; then
  "$PY" "$BIN/manifest.py" write-run "$RUNDIR" batch_dir="$BATCHDIR" \
    batch_config_hash="$BATCH_CFG_HASH" config_hash="$RUN_CFG_HASH" replicates="$NREPS" \
    specimens="$SPECIMENS" simfbd_commit="$(commit_of "$ROOT")" \
    infer_script_hash="$(hash_of "$TEMPLATE")" rb_path="$rb_bin" rb_md5="$rb_md5" \
    rb_commit="$rb_commit" || die "could not write $RUN_MANIFEST"
else
{
  printf 'kind\trun\n'
  printf 'run\t%s\n' "$(basename "$RUNDIR")"
  printf 'batch\t%s\n' "$(basename "$BATCHDIR")"
  printf 'batch_dir\t%s\n' "$BATCHDIR"
  printf 'batch_config_hash\t%s\n' "$BATCH_CFG_HASH"
  printf 'config_hash\t%s\n' "$RUN_CFG_HASH"
  # the scoring needs it and must not have to guess; a batch predating the key is oldest-first
  tvo="$(awk -F'\t' '$1=="true_vals_order"{print $2}' "$BATCH_MANIFEST" | head -1)"
  printf 'true_vals_order\t%s\n' "${tvo:-oldest_first}"
  printf 'infer_script_hash\t%s\n' "$(hash_of "$TEMPLATE")"
  printf 'survivors\t%s\n' "$SURVIVORS"
  printf 'drop_present_only\t%s\n' "$DROP_PRESENT_ONLY"
  printf 'simfbd_commit\t%s\n' "$(commit_of "$ROOT")"
  printf 'rb_path\t%s\n' "$rb_bin"
  printf 'rb_md5\t%s\n' "$rb_md5"
  printf 'rb_commit\t%s\n' "$rb_commit"
  printf 'created\t%s\n' "$(date -Is)"
  printf 'specimens\t%s\n' "$SPECIMENS"
  for k in MODEL REPORTING INFER COND RHO GENS PRINTGEN ORIGIN_PRIOR; do printf '%s\t%s\n' "$k" "${!k}"; done
  # the inherited hyperpriors, so a run manifest describes its own fit without the batch
  for k in $PRIOR_KEYS; do printf '%s\t%s\n' "$k" "${!k}"; done
  for k in $OPT_KEYS; do printf '%s\t%s\n' "$k" "${!k:-none}"; done
} > "$RUN_MANIFEST"
fi

nfail=$(wc -l < "$RUNDIR/failures.log")
echo "logs: $(ls "$OUTDIR"/rep_*.log 2>/dev/null | wc -l)/$NREPS, $nfail failed"
echo "manifest: $RUN_MANIFEST"
echo "now run: simfbd summarize $RUNDIR"
