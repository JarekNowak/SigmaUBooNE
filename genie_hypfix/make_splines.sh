#!/bin/bash
# Make HypQELTracePXSec splines (LambdaCCQE channels) for every target in targets.txt and the three
# antineutrino flavours, with the knot grid of the production file (250 knots, up to 1000 GeV).
# Usage (inside the SL7 container): make_splines.sh <outdir> <njobs>
O="$1"; NJ="${2:-8}"
H="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$H/../../setup_env.sh" >/dev/null 2>&1
source "$H/env_fix.sh"
cd "$O" || exit 1
split -n l/$NJ -d targets.txt parts/tg_
for f in parts/tg_*; do
  t=$(paste -sd, "$f")
  ( gmkspl -p -12,-14,-16 -t "$t" -n 250 -e 1000 --event-generator-list LambdaCCQE \
           --tune G18_10a_02_11a -o "$f.xml" > "$f.log" 2>&1; echo "rc=$?" >> "$f.log" ) &
done
wait
echo ALLDONE
