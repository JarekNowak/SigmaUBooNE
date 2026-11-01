#!/bin/bash
# Run on lapx: copy a prestaged sample of samples.tsv from FNAL dCache to local disk over xrootd.
# Resumable (files already here are skipped; a .part file is renamed only after a complete copy).
#
#   ./fetch_sample.sh <key> [-N n] [-P par] [-M max] [--unblind] [--relist]
#     -N n       the prestaged sub-definition with limit n (gpvm_prestage.sh <key> n); default: the full one.
#                If the sub-definition does not exist yet, the sample's own dimensions are used (only
#                files already on disk at FNAL are copied).
#     -P par     parallel copies (default 4)
#     -M max     stop after max new copies (0 = only build the file list and report)
#     --relist   rebuild the pnfs list from SAM (default: reuse <dest>/pnfs.list if present)
#     --unblind  allow beam-on keys (blind analysis: only after unblinding is approved)
# Destination: $DEST_ROOT/<key> (default /data/uboone/SigmaFiles/transfer/<key>); writes
#   pnfs.list  (SAM -> dCache paths), files.list (complete local files), .xrdcp_errors.log
#
# Credentials: a uboone bearer token, from ../../production/nuwro_prod/get_token.sh (renew every ~3 h;
# a long fetch should be run in a loop that renews it, see README). Files not yet on disk at FNAL are
# reported as tape-only and skipped: prestage them first (gpvm_prestage.sh), then re-run.
set -u
TDIR="$(cd "$(dirname "$0")" && pwd)"
source "$TDIR/sample_lookup.sh"
CONTAINER="/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest"
N=""; PAR=4; MAX=-1; UNBLIND=0; RELIST=0; KEY=""
while [ $# -gt 0 ]; do
  case "$1" in
    -N) N="$2"; shift ;; -P) PAR="$2"; shift ;; -M) MAX="$2"; shift ;;
    --unblind) UNBLIND=1 ;; --relist) RELIST=1 ;;
    -h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) KEY="$1" ;;
  esac; shift
done
[ -n "$KEY" ] || { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
lookup "$KEY" || { echo "unknown key $KEY"; exit 1; }
if [ "$S_BLIND" = 1 ] && [ $UNBLIND = 0 ]; then echo "REFUSED: $KEY is beam-on data (blind analysis); needs --unblind"; exit 1; fi

DEST="$DEST_ROOT/$KEY"; mkdir -p "$DEST"
# keep >= 5 TB free on /data
free=$(df -B1T --output=avail /data | tail -1 | tr -d ' ')
[ "$free" -ge 5 ] || { echo "ERROR: only ${free} TB free on /data"; exit 1; }
exec 8> "$DEST/.fetch.lock"; flock -n 8 || { echo "ERROR: another fetch is copying into $DEST"; exit 1; }

SUBDEF="$(subdef "$KEY" "$N")"; RAW="$S_DIMS${N:+ with limit $N}"
[ -s "$DEST/pnfs.list" ] && [ $RELIST = 0 ] && SUBDEF="" && RAW=""
if [ -z "${BEARER_TOKEN_FILE:-}" ]; then
  for tf in "/run/user/$(id -u)/bt_u$(id -u)" "/tmp/bt_u$(id -u)"; do [ -s "$tf" ] && { export BEARER_TOKEN_FILE="$tf"; break; }; done
fi
[ -n "${BEARER_TOKEN_FILE:-}" ] || echo "WARNING: no bearer token found (run get_token.sh)"
echo "$(date +%T) $KEY ($S_TIER) -> $DEST"
apptainer exec --env BEARER_TOKEN_FILE="${BEARER_TOKEN_FILE:-}" \
  -B /cvmfs -B /data -B /home -B /run/user -B /tmp "$CONTAINER" \
  bash "$TDIR/fetch_worker.sh" "$DEST" "$SUBDEF" "$RAW" "$MAX" "$PAR"
while read -r p; do f="$DEST/$(basename "$p")"; [ -s "$f" ] && echo "$f"; done < "$DEST/pnfs.list" | sort > "$DEST/files.list"
echo "$(date +%T) local: $(wc -l < "$DEST/files.list") of $(wc -l < "$DEST/pnfs.list") -> $DEST/files.list"
