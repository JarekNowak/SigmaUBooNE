#!/bin/bash
# Run on lapx (in screen/tmux or with nohup): keep fetching one or more samples while their
# prestage proceeds at FNAL. Each cycle renews the bearer token (silent renewal from the vault
# token, valid ~1 week after one browser login) and runs fetch_sample.sh on every key; it stops
# when every listed file is local, or after MAXCYCLES.
#   nohup ./fetch_loop.sh dirt_rhc_run3b dv_hyp_rhc_CV > fetch_loop.log 2>&1 &
# Environment: SLEEP (s between cycles, default 1800), MAXCYCLES (default 200), PAR (default 4)
TDIR="$(cd "$(dirname "$0")" && pwd)"
source "$TDIR/sample_lookup.sh"
GETTOKEN="$TDIR/../../production/nuwro_prod/get_token.sh"
SLEEP="${SLEEP:-1800}"; MAXCYCLES="${MAXCYCLES:-200}"; PAR="${PAR:-4}"
[ $# -ge 1 ] || { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
for c in $(seq 1 "$MAXCYCLES"); do
  echo "=== cycle $c $(date '+%F %T')"
  "$GETTOKEN" < /dev/null > /dev/null 2>&1 || echo "WARNING: token renewal failed (log in again: $GETTOKEN)"
  done_all=1
  for k in "$@"; do
    "$TDIR/fetch_sample.sh" "$k" -P "$PAR" | grep -v '^$'
    d="$DEST_ROOT/$k"
    [ -s "$d/pnfs.list" ] && [ "$(wc -l < "$d/files.list")" -eq "$(wc -l < "$d/pnfs.list")" ] || done_all=0
  done
  [ $done_all = 1 ] && { echo "all complete"; exit 0; }
  sleep "$SLEEP"
done
