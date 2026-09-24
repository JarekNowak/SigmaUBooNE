// Compiled wrapper: ACLiC builds this together with SigmaSelection.C, which is the
// invocation path that works here (.L from inside a macro trips a TMVA/BLAS load).
#include "SigmaSelection.C"
void runNoPion2(){
  SigmaSelection("SigmaHists_nopion.root", true, false);
}
