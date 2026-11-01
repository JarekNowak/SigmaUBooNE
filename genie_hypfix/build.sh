#!/bin/bash
# Build libGHypFix.so (+ dictionary and rootmap) against the production GENIE
# (uboonecode v08_00_00_85: GENIE v3_00_04_ub4, ROOT 6.12, gcc 7.3).
# Run inside the SL7 container:  bash build.sh
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../../setup_env.sh" >/dev/null 2>&1   # before set -e: the ups scripts fail under it
set -e
cd "$HERE"; mkdir -p lib
INC="-I$GENIE_INC/GENIE -I$(root-config --incdir) -I$LOG4CPP_INC -I$HERE"
rootcling -f lib/HypFixDict.cxx -rmf lib/libGHypFix.rootmap -rml libGHypFix.so \
          $INC HypQELTracePXSec.h LinkDef.h
g++ -O2 -fPIC -shared $(root-config --cflags) $INC HypQELTracePXSec.cxx lib/HypFixDict.cxx \
    -L$GENIE_LIB -lGPhQELXS -lGPhXSIg -lGFwEG -lGFwInt -lGFwAlg -lGFwParDat -lGFwUtl -lGFwMsg -lGFwReg \
    -L$LOG4CPP_LIB -llog4cpp -L$LIBXML2_FQ_DIR/lib -lxml2 $(root-config --libs) -o lib/libGHypFix.so
mv lib/HypFixDict_rdict.pcm lib/ 2>/dev/null || true
ls -la lib
