#!/bin/bash
# Run on lapx: HyperonProduction (hypana ntuples) on every transferred reco2 file of a sample of
# samples.tsv that has no ntuple yet. One SLURM array task = one file; idempotent (re-run to pick
# up failures or newly transferred files).
#   ./submit_hypana_sample.sh <key> [-n N] [--dry-run]
# Input:  <reco2 dir>/files.list (fetch_sample.sh)  e.g. /data/uboone/dirt/reco2/dirt_rhc_run3b
# Output: $DATA_ROOT/<category>/hypana/<key>/<base>_hypana.root (+ .log), slurm logs in slurm_logs/
# Config: HYPANA_FCL (default run_HyperonProduction.fcl, the overlay campaign config; beam-off data
#         needs run_HyperonProduction_data.fcl, see production/ext_hypana/README.md).
# Environment: HYP_MAXRUN (concurrent tasks, default 100), CAMPAIGN_NICE (default 10000).
set -euo pipefail
TDIR="$(cd "$(dirname "$0")" && pwd)"
source "$TDIR/sample_lookup.sh"
RUNNER=/home/t2k/nowak/MicroBooNE/production/hyperon_ana_local/run_hypana_file.sh
KEY=""; NMAX=0; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in -n) NMAX="$2"; shift ;; --dry-run) DRY=1 ;; -h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;; *) KEY="$1" ;; esac; shift
done
[ -n "$KEY" ] || { grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
lookup "$KEY" || { echo "unknown key $KEY"; exit 1; }
[ "$S_TIER" = reco2 ] || { echo "$KEY is $S_TIER: run our reco2 first"; exit 1; }
if [ "$S_BLIND" = 1 ]; then echo "REFUSED: $KEY is beam-on data (blind analysis)"; exit 1; fi
IN="$(destdir "$KEY")"; OUT="$DATA_ROOT/$(category "$KEY")/hypana/$KEY"
[ -s "$IN/files.list" ] || { echo "no $IN/files.list (fetch first)"; exit 1; }
mkdir -p "$OUT/slurm_logs" "$OUT/manifests"
FCL="${HYPANA_FCL:-run_HyperonProduction.fcl}"

# not yet done and not in a live array of this sample
live=$(squeue -h -u "$USER" -n "hyp_$KEY" -o %F | sort -u)
inflight=$(for a in $live; do awk -v a="$a" '$1==a{print $2}' "$OUT/manifests/arrays.tsv" 2>/dev/null | xargs -r cat; done | sort -u)
todo=$(mktemp)
while read -r f; do
  b=$(basename "$f" .root)
  [ -s "$OUT/${b}_hypana.root" ] && continue
  grep -qxF "$f" <<< "$inflight" && continue
  echo "$f"
done < "$IN/files.list" > "$todo"
[ "$NMAX" -gt 0 ] && { head -n "$NMAX" "$todo" > "$todo.t"; mv "$todo.t" "$todo"; }
N=$(wc -l < "$todo")
echo "$KEY: $N of $(wc -l < "$IN/files.list") reco2 files need hypana ($FCL) -> $OUT"
[ "$N" -gt 0 ] && [ $DRY = 0 ] || { rm -f "$todo"; exit 0; }

for (( s=1; s<=N; s+=1000 )); do         # MaxArraySize 1001
  m="$OUT/manifests/hyp_$(date +%Y%m%dT%H%M%S)_$s.list"; sed -n "${s},$((s+999))p" "$todo" > "$m"; n=$(wc -l < "$m")
  id=$(sbatch --parsable --nice="${CAMPAIGN_NICE:-10000}" --job-name="hyp_$KEY" --partition=main \
        --cpus-per-task=1 --mem=8G --time=4:00:00 --array=0-$((n-1))%${HYP_MAXRUN:-100} \
        --output="$OUT/slurm_logs/%x_%A_%a.out" \
        --wrap="IN=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID+1))p\" $m); echo host=\$(hostname) input=\$IN; \
apptainer exec --env HYPANA_FCL=$FCL -B /cvmfs -B /data -B /home \
/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest bash $RUNNER \"\$IN\" $OUT")
  echo "$id $m" >> "$OUT/manifests/arrays.tsv"
  echo "array $id: $n files"
  sleep 1
done
rm -f "$todo"
