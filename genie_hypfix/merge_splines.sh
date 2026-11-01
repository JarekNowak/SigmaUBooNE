#!/bin/bash
# Production spline file + the HypQELTracePXSec splines -> gxspl-FNALsmall-hypfix.xml.
# The Pais splines stay in the file (different key, unused once ModelConfiguration points at
# HypQELTracePXSec). Usage: merge_splines.sh <outdir>
O="$1"
SRC=/cvmfs/uboone.opensciencegrid.org/products/genie_xsec/v3_00_04_ub2/NULL/G1810a0211a-k250-e1000/data/gxspl-FNALsmall.xml
OUT="$O/gxspl-FNALsmall-hypfix.xml"
for f in "$O"/parts/tg_*.log; do grep -q 'rc=0' "$f" || { echo "ERROR: $f did not finish cleanly" >&2; exit 1; }; done
{
  sed '/<\/genie_tune>/,$d' "$SRC"
  for f in "$O"/parts/tg_*.xml; do sed -n '/<spline /,/<\/spline>/p' "$f"; done
  sed -n '/<\/genie_tune>/,$p' "$SRC"
} > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
echo "Pais splines:             $(grep -c 'PaisQELLambdaPXSec' "$OUT")"
echo "HypQELTracePXSec splines: $(grep -c 'HypQELTracePXSec' "$OUT")"
echo "total splines:            $(grep -c '<spline ' "$OUT") (source $(grep -c '<spline ' "$SRC"))"
