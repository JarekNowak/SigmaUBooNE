// Run the Sigma0 selection WITHOUT the pion-candidate requirement, so the BDT
// comparison can isolate the effect of that cut. Same feature set: events with
// no third track are kept and their pion inputs carry sentinels.
void runNoPion(){
  gROOT->ProcessLine(".L SigmaSelection.C+");
  gROOT->ProcessLine("SigmaSelection(\"SigmaHists_nopion.root\", true, false)");
}
