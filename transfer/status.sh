#!/bin/bash
# Run on lapx: per sample of samples.tsv, local files against the pnfs list (and SAM, with -s).
#   ./status.sh          local counts only (fast)
#   ./status.sh -s       also the SAM file count and size of each sample (slow; container)
source "$(dirname "$0")/sample_lookup.sh"
SAM=0; [ "${1:-}" = "-s" ] && SAM=1
printf "%-28s %-5s %-5s %8s %8s %9s\n" key block blind local listed localTB
for k in $(list_keys); do
  lookup "$k"; d="$(destdir "$k")"
  nl=$( [ -s "$d/files.list" ] && wc -l < "$d/files.list" || echo 0)
  np=$( [ -s "$d/pnfs.list" ] && wc -l < "$d/pnfs.list" || echo -)
  tb=$( [ -d "$d" ] && du -s --block-size=1G "$d" 2>/dev/null | awk '{printf "%.2f", $1/1000}' || echo 0)
  printf "%-28s %-5s %-5s %8s %8s %9s\n" "$k" "$S_BLOCK" "$([ "$S_BLIND" = 1 ] && echo BLIND)" "$nl" "$np" "$tb"
done
if [ $SAM = 1 ]; then
  apptainer exec -B /cvmfs -B /home /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest bash -c '
    source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh >/dev/null 2>&1; setup sam_web_client >/dev/null 2>&1
    awk -F"\t" "!/^#/ && NF>=5 {print \$1 \"\t\" \$5}" "'"$TABLE"'" | while IFS=$'"'\t'"' read -r k dims; do
      echo "$k | $(samweb -e uboone list-files --summary "$dims" | tr "\n" " ")"; done'
fi
