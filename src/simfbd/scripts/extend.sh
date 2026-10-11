#!/bin/bash
# Add replicates to a batch.
#
#   simfbd extend <batch> 200      # a new batch gets its first 200
#   simfbd extend <batch> 800      # and now it has 1000
#
# The batch's own config.sh is the only source of settings, so the new replicates cannot be
# drawn under different ones, and a config or a simulator that has changed since the earlier
# replicates is refused. Those are never touched: an extension simulates only what is missing,
# appends to true_vals.tsv and nums.tsv, and writes its per-replicate files beside the others.
#
# Runs under the batch pick the new replicates up the next time they are run, since a run
# skips every replicate whose log is already complete.
set -u
BATCH="${1:?usage: extend.sh sims/<batch> <how many>}"; BATCH="${BATCH%/}"
N="${2:?usage: extend.sh sims/<batch> <how many>}"
case "$BATCH" in */runs/*) echo "error: $BATCH is a run; extend its batch" >&2; exit 1 ;; esac
case "$N" in ''|*[!0-9]*|0) echo "error: how many must be a positive whole number, not '$N'" >&2; exit 1 ;; esac
EXTEND="$N" exec bash "$(dirname "$0")/run.sh" "$BATCH"
