#!/bin/bash
# Run SigmaSelection.C in detector-variation mode on the run 3 RHC CV and its nine
# variations (NJOBS at a time), then compare them with DetVar.C.
#   ./run_detvar.sh [NJOBS]     outputs SigmaHists_dv_<name>.root, logs in dvlogs/
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=$HOME/lib/flexiblas-shim:$LD_LIBRARY_PATH
D=/data/uboone/temp/custom_pelee_ntuples/detvar/run3
NJOBS=${1:-3}
mkdir -p dvlogs
# compile once, so the parallel jobs only load the library
root -l -b -q -e '.L SigmaSelection.C+' > dvlogs/compile.log 2>&1 || { echo "compile failed"; exit 1; }
printf "%s\n" CV LYAttenuation LYDown LYRayleigh Recombination SCE WireModThetaXZ WireModThetaYZ WireModX WireModYZ |
  xargs -P "$NJOBS" -I{} sh -c "root -l -b -q 'SigmaSelection.C+(\"SigmaHists_dv_{}.root\",false,true,false,\"$D/numi_detvars_pion_ntuples_run3_rhc_{}_slim.root\")' > dvlogs/{}.log 2>&1"
root -l -b -q DetVar.C 2>&1 | tee dvlogs/DetVar.txt
