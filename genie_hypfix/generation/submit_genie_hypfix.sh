#!/bin/bash
# GENIE Sigma0 QE samples with the FIXED hyperon cross section (genie::HypQELTracePXSec,
# hypfix/README.md), as SLURM arrays of chunks run by gen_genie_chunk_hypfix.sh.
# Same generator settings, chunk size and chunk counts as the original samples
# (submit_genie.sh): RHC 3050 x 500, FHC 1200 x 500; new sample names with _hypfix.
#
# At most MAXRUN jobs run at once over the WHOLE submission: the arrays (<= 1000 tasks each,
# MaxArraySize=1001) are chained with afterany dependencies and each carries %MAXRUN.
# -d JOBID makes the first array wait for JOBID (e.g. FHC after the last RHC array).
#
# Usage: ./submit_genie_hypfix.sh -b rhc|fhc [-k FIRST (0)] [-n CHUNKS (3050 rhc, 1200 fhc)]
#                                 [-e EVENTS_PER_CHUNK (500)] [-m MAXRUN (100)] [-d JOBID] [--dry-run]
# Prints the job id of the last array (for -d of the next submission).
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIRST=0; NCHUNK=""; NEV=500; DRY=0; MAXRUN=100; BEAM=""; DEP=""
while [ $# -gt 0 ]; do
  case "$1" in -b) BEAM="$2"; shift 2 ;; -k) FIRST="$2"; shift 2 ;; -m) MAXRUN="$2"; shift 2 ;; -n) NCHUNK="$2"; shift 2 ;;
    -e) NEV="$2"; shift 2 ;; -d) DEP="$2"; shift 2 ;; --dry-run) DRY=1; shift ;;
    *) echo "unknown option $1" >&2; exit 1 ;; esac
done
case "$BEAM" in
  rhc) NAME=genie_uboone_NuMI_RHC_sigmazero_run3b_2026_hypfix; NCHUNK=${NCHUNK:-3050} ;;
  fhc) NAME=genie_uboone_NuMI_FHC_sigmazero_run1_2026_hypfix;  NCHUNK=${NCHUNK:-1200} ;;
  *) echo "-b must be rhc or fhc" >&2; exit 1 ;;
esac
FCL=genie_lambdaccqe_numi_${BEAM}_hypfix.fcl
[ -f "$SCRIPTDIR/$FCL" ] || { echo "missing $FCL" >&2; exit 1; }
[ -s /data/uboone/genie_samples/hypfix_splines/gxspl-FNALsmall-hypfix.xml ] || { echo "missing hypfix splines" >&2; exit 1; }
[ -s "$SCRIPTDIR/hypfix/lib/libGHypFix.so" ] || { echo "missing libGHypFix.so (hypfix/build.sh)" >&2; exit 1; }
BASE=/data/uboone/genie_samples/$NAME
CONTAINER=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
mkdir -p "$BASE/chunks" "$BASE/logs"
LAST=$(( FIRST + NCHUNK - 1 ))
old=$(grep -oE '^chunks=[0-9]+' "$BASE/sample.conf" 2>/dev/null | cut -d= -f2); [ "${old:-0}" -gt $((LAST+1)) ] && LASTN=$old || LASTN=$((LAST+1))
printf 'name=%s\nchunks=%s\nevents_per_chunk=%s\nfcl=%s\nselect="primary hyperon 3212"\nxsec=genie::HypQELTracePXSec (hypfix)\n' \
  "$NAME" "$LASTN" "$NEV" "$FCL" > "$BASE/sample.conf"
prev="$DEP"
for (( off=FIRST; off<=LAST; off+=1000 )); do
  end=$(( off + 999 )); [ $end -gt $LAST ] && end=$LAST
  todo=$(for ((i=off; i<=end; i++)); do [ -e "$BASE/chunks/$(printf c%04d $i)/DONE" ] || echo $((i - off)); done | paste -sd,)
  [ -n "$todo" ] || { echo "chunks $off..$end all DONE" >&2; continue; }
  cmd=(sbatch --parsable --nice="${CAMPAIGN_NICE:-10000}" --exclude="${CAMPAIGN_EXCLUDE-py-batch[13-16]}"
       --job-name=geniefix${BEAM} --partition=main --cpus-per-task=1 --mem=2G --time=48:00:00
       --array="${todo}%${MAXRUN}" --output="$BASE/logs/%x_%A_%a.out")
  [ -n "$prev" ] && cmd+=(--dependency=afterany:$prev)
  cmd+=(--wrap "c=\$(printf c%04d \$((SLURM_ARRAY_TASK_ID + $off))); apptainer exec -B /cvmfs -B /data -B /home $CONTAINER bash $SCRIPTDIR/gen_genie_chunk_hypfix.sh $BASE/chunks/\$c $NEV $FCL")
  echo "chunks $off..$end: ${cmd[*]}" | cut -c1-300 >&2
  if [ $DRY -eq 1 ]; then prev="DRY$off"; else prev=$("${cmd[@]}") || exit 1; echo "  -> job $prev" >&2; fi
done
echo "$prev"
