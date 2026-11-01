#!/bin/bash
# Run on a uboone gpvm (uboonegpvmXX.fnal.gov), NOT on lapx: creates a SAM sub-definition for a
# sample of samples.tsv and prestages it from tape to dCache, so fetch_sample.sh on lapx can copy it
# (a lapx user token cannot recall from tape: "Staging not allowed").
#
#   ./gpvm_prestage.sh <key> [N]       sub-definition <SAMUSER>_sigma_<key>[_N] (first N files), prestage
#   ./gpvm_prestage.sh -b R            every key of block R (R, F, E, O; D needs --unblind)
#   ./gpvm_prestage.sh -s <key> [N]    only create the sub-definition and print its summary
#   ./gpvm_prestage.sh -l              list the keys
# Options: --unblind (allow beam-on keys), -p N (prestage parallelism, default 4)
#
# Setup on the gpvm first:
#   source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh; setup sam_web_client
#   htgettoken -a htvaultprod.fnal.gov -i uboone     # or kx509, as the gpvm requires
# Prestaging runs in the background (nohup); logs in ./prestage_logs/<subdef>.log. A recall of
# several TB takes hours to days; dCache keeps the copies on disk for a limited time, so copy them
# soon after (fetch_sample.sh), in priority order, rather than prestaging everything at once.
set -u
source "$(dirname "$0")/sample_lookup.sh"
UNBLIND=0; PAR=4; ONLYDEF=0; BLOCK=""
args=()
while [ $# -gt 0 ]; do
  case "$1" in
    --unblind) UNBLIND=1 ;; -p) PAR="$2"; shift ;; -s) ONLYDEF=1 ;; -b) BLOCK="$2"; shift ;;
    -l) awk -F'\t' '!/^#/ && NF>=5 {printf "%-28s %s %s %s\n",$1,$2,($3==1?"BLIND":"     "),$4}' "$TABLE"; exit 0 ;;
    -h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) args+=("$1") ;;
  esac; shift
done
command -v samweb >/dev/null || { echo "ERROR: samweb not set up (see the header)"; exit 1; }
mkdir -p prestage_logs

one() {   # key [N]
  local key="$1" n="${2:-}"
  lookup "$key" || { echo "unknown key $key (./gpvm_prestage.sh -l)"; return 1; }
  if [ "$S_BLIND" = 1 ] && [ $UNBLIND = 0 ]; then echo "SKIP $key: beam-on data, blind (needs --unblind)"; return 0; fi
  local def; def=$(subdef "$key" "$n")
  local dims="$S_DIMS"; [ -n "$n" ] && dims="$dims with limit $n"
  if samweb describe-definition "$def" >/dev/null 2>&1; then echo "exists: $def"
  else samweb create-definition "$def" "$dims" && echo "created: $def = $dims" || return 1; fi
  echo "  $(samweb list-files --summary "defname: $def" | tr '\n' ' ')"
  [ $ONLYDEF = 1 ] && return 0
  nohup samweb prestage-dataset --defname="$def" --parallel="$PAR" > "prestage_logs/$def.log" 2>&1 &
  echo "  prestaging in the background (pid $!), log prestage_logs/$def.log"
}

if [ -n "$BLOCK" ]; then for k in $(list_keys "$BLOCK"); do one "$k"; done
elif [ ${#args[@]} -ge 1 ]; then one "${args[@]}"
else grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 1; fi
