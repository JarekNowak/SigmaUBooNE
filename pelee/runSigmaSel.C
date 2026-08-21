// runSigmaSel.C — compile and run the PeLEE-based Sigma0 selection.
//   root -l -b -q runSigmaSel.C
void runSigmaSel() {
  gROOT->ProcessLine(".L SigmaSelection.C+");
  gROOT->ProcessLine("SigmaSelection(\"SigmaHists.root\")");
}
