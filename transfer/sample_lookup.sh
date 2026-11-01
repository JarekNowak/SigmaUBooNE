#!/bin/bash
# Shared by the transfer scripts: look a sample up in samples.tsv.
#   source sample_lookup.sh; lookup <key>   -> sets S_BLOCK S_BLIND S_TIER S_DIMS (returns 1 if unknown)
#   list_keys [block]                        -> prints the keys (of one block)
TDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TABLE="$TDIR/samples.tsv"
SAMUSER="${SAMUSER:-janow}"                         # FNAL user name: prefix of the sub-definitions
DEST_ROOT="${DEST_ROOT:-/data/uboone/SigmaFiles/transfer}"

lookup() {
  local line; line=$(awk -F'\t' -v k="$1" '!/^#/ && $1==k' "$TABLE")
  [ -n "$line" ] || return 1
  IFS=$'\t' read -r _ S_BLOCK S_BLIND S_TIER S_DIMS <<< "$line"
}
list_keys() { awk -F'\t' -v b="${1:-}" '!/^#/ && NF>=5 && (b=="" || $2==b) {print $1}' "$TABLE"; }
subdef() { echo "${SAMUSER}_sigma_$1${2:+_$2}"; }   # sub-definition name (optional limit suffix)
