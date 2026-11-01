#!/bin/bash
# Generate ONE chunk of GENIE QE-hyperon events with the hyperon cross-section fix
# (genie::HypQELTracePXSec, hypfix/README.md) -- otherwise identical to gen_genie_chunk.sh.
# A separate file so that jobs running gen_genie_chunk.sh are never edited under them.
# (LambdaCCQE, official NuMI RHC overlay gen
# settings, no EXT) and convert the Sigma0-primary ones to hepmc (runs inside the SL7 container;
# submitted by submit_genie.sh).
# Usage: gen_genie_chunk_hypfix.sh <workdir> <nevents> [fcl (genie_lambdaccqe_numi_fhc_hypfix.fcl)]
# Output in <workdir>: genie_hyp_gen.root, gen.log, chunk.hepmcev, chunk.summary, DONE
WORK="$1"; NEV="$2"; FCL="${3:-genie_lambdaccqe_numi_fhc_hypfix.fcl}"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPTDIR/../setup_env.sh" >/dev/null 2>&1 || { echo "ERROR: uboonecode setup" >&2; exit 1; }
export FHICL_FILE_PATH="$SCRIPTDIR:$FHICL_FILE_PATH"
source "$SCRIPTDIR/hypfix/env_fix.sh" || { echo "ERROR: hypfix env" >&2; exit 1; }   # after setup_env.sh
mkdir -p "$WORK" && cd "$WORK" || exit 1
rm -f DONE
echo "$(date '+%F %T') GENIE LambdaCCQE ($FCL): $NEV events on $(hostname)"
# a rerun after a timeout in the conversion step reuses a complete GENIE output
if [ -s genie_hyp_gen.root ] && grep -q "TrigReport Events total = $NEV " gen.log 2>/dev/null \
   && grep -q "Art has completed and will exit with status 0" gen.log; then
  rc=0; echo "$(date '+%F %T') reusing complete genie_hyp_gen.root"
else
  lar -c "$FCL" -n "$NEV" -o genie_hyp_gen.root -T genie_hyp_gen_hist.root > gen.log 2>&1
  rc=$?; echo "$(date '+%F %T') lar rc=$rc"
fi
[ $rc -eq 0 ] && grep -q "TrigReport Events total = $NEV " gen.log || { tail -20 gen.log >&2; exit 2; }
root -l -b -q "$SCRIPTDIR/genie_to_hepmc.C(\"$WORK/genie_hyp_gen.root\", \"$WORK/chunk\")" > convert.log 2>&1
rc=$?; echo "$(date '+%F %T') convert rc=$rc  $(tr '\n' ' ' < chunk.summary 2>/dev/null)"
[ $rc -eq 0 ] && [ -s chunk.summary ] || { tail -20 convert.log >&2; exit 3; }
rm -f Pandora_Events.pndr mem.db time.db
touch DONE
