// TrainBDTSplit.C -- dedicated per-background-class BDTs, replacing the single
// discriminant trained against every background at once.
//
// Motivation (sigma0_note.tex, "Plan: dedicated BDTs"): the background is not one
// thing. It is pi0-producing RES and DIS, non-pi0 neutrino interactions, and cosmics,
// which look nothing like each other, so a single surface separating signal from all
// of them learns what they share and averages away what distinguishes each. gLEE runs
// five BDTs for its 1g1p channel for exactly this reason.
//
// The evidence that this is the binding constraint rather than variable quality: the
// Moliere angle moved the figure of merit by x1.013, and the six Second-Shower-Veto
// variables moved it from 0.0537 to 0.0540 -- nothing -- despite being aimed squarely
// at the 76% of background that is pi0-driven. Good variables in the wrong-shaped
// discriminant.
//
// Three classes, from truth branches already written to BDTInput.root:
//   pi0     HasPi0>0.5                        37.8% of background
//   cosmic  CosmicMu>0.5                       4.3%
//   nunopi0 HasPi0<0.5 && CosmicMu<0.5        59.6%
// Signal is the same TreeS in all three; only the background changes. The two-fold
// parity cross-validation is preserved per class, so six trainings in total.
//
//   usage: root -l -b -q TrainBDTSplit.C
//   output: SigmaBDT_<class>_<A|B>/weights/

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include <iostream>

static const char* kVars[] = {
  "MuonPID","MuonTrkLength","MuonTrkDistVtx","MuonTrkScore","MuonNHitsY",
  "PionPID","PionTrkLength","PionTrkDistVtx","PionTrkScore","PionNHitsY","PionCaloEnergyY",
  "ProtonPID","ProtonTrkLength","ProtonTrkDistVtx","ProtonTrkScore","ProtonNHitsY","ProtonCaloEnergyY","ProtonChi2",
  "Shower1NHitsY","Shower1VtxDist","Shower1DeDxY","Shower1Energy",
  "DistMuonPion","DistProtonPion","DistProtonMuon","DistShower1Muon","DistShower1Pion","DistShower1Proton",
  "AngleMuonPion","AngleProtonPion","AngleProtonMuon","AngleShower1Muon","AngleShower1Pion","AngleShower1Proton",
  "LambdaMass","SigmaMass","TopologicalScore","NTracks","NShowers",
  "LambdaPointAngle","LambdaImpact","LambdaDecayDist","PhotonPointAngle","PhotonImpact",
  "LambdaGammaAngle","SigmaMuonAngle",
  "Pi0Mass","Pi0Gammadot","Shower2Energy",
  // Shower Moliere angle -- must stay in the SAME ORDER as BDTV in
  // SigmaSelection.C, which keeps its own copy of this list. The reader there
  // registers variables positionally, so a mismatch silently feeds the BDT the
  // wrong numbers rather than failing.
  "Shower1MoliereAvg",
  // Second-Shower Veto block -- same order as BDTV in SigmaSelection.C.
  // SliceUnassocFrac/Hits excluded -- NoCR artefact, see SigmaSelection.C.
  "SliceClustFrac","SecondShowerNHit","SecondShowerDot","ShrSubClusters",
  "UnassocHitsCalc","UnassocFracCalc" };


// One (class, fold) training. `bkgCut` selects the background class; `foldCut`
// selects the events this fold trains on. They are ANDed for the background and the
// fold cut alone is used for the signal, which is common to every class.
static void trainOne(const char* dataset, const char* foldCut, const char* bkgCut,
                     TTree* sig, TTree* bkg)
{
  TFile* out = TFile::Open(Form("TMVA_%s.root", dataset), "RECREATE");
  TMVA::Factory factory("TMVASigma", out,
    "!V:!Silent:Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification");
  TMVA::DataLoader loader(dataset);

  loader.AddSignalTree(sig, 1.0);
  loader.AddBackgroundTree(bkg, 1.0);
  loader.SetSignalWeightExpression("Weight");
  loader.SetBackgroundWeightExpression("Weight");
  for (auto v : kVars) loader.AddVariable(v, 'F');
  loader.AddSpectator("Weight", 'F');
  loader.AddSpectator("IsOdd", 'F');

  // Signal keeps only the fold cut; background additionally keeps its class.
  TCut sigCut = foldCut;
  TCut bkgSel = TCut(foldCut) && TCut(bkgCut);
  loader.PrepareTrainingAndTestTree(sigCut, bkgSel,
    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

  // Same hyper-parameters as the single BDT, so any difference in performance is the
  // class split and not a retune. MinNodeSize is left at 3%: with the background now
  // as small as ~1500 events for the cosmic class, a looser value would overtrain.
  factory.BookMethod(&loader, TMVA::Types::kBDT, "BDTG",
    "!H:!V:NTrees=600:MinNodeSize=3%:MaxDepth=3:BoostType=Grad:Shrinkage=0.06:"
    "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:NegWeightTreatment=Pray");

  factory.TrainAllMethods();
  factory.TestAllMethods();
  factory.EvaluateAllMethods();
  out->Close();
  std::cout << "  " << dataset << " done -> " << dataset << "/weights/\n";
}

void TrainBDTSplit()
{
  TMVA::Tools::Instance();
  TFile* in = TFile::Open("BDTInput.root", "READ");
  if (!in || in->IsZombie()) { std::cerr << "cannot open BDTInput.root\n"; return; }
  TTree* sig = (TTree*)in->Get("TreeS");
  TTree* bkg = (TTree*)in->Get("TreeB");
  if (!sig || !bkg) { std::cerr << "TreeS/TreeB missing\n"; return; }

  struct Cls { const char* name; const char* cut; };
  const std::vector<Cls> classes = {
    { "pi0",     "HasPi0>0.5" },
    { "nunopi0", "HasPi0<0.5 && CosmicMu<0.5" },
    { "cosmic",  "CosmicMu>0.5" },
  };

  std::cout << "Signal=" << sig->GetEntries()
            << "  Background=" << bkg->GetEntries() << "\n";
  for (const auto& c : classes)
    std::cout << "  class " << c.name << " : " << bkg->GetEntries(c.cut) << " events\n";

  // Refuse to train a class with too little background to be meaningful. 500 is the
  // same floor used for the TKI binning in the CC1mu1pi analysis, and is well below
  // where TMVA becomes unstable, but it stops a silent near-empty training.
  for (const auto& c : classes) {
    Long64_t n = bkg->GetEntries(c.cut);
    if (n < 500) {
      std::cout << "  SKIP " << c.name << " (only " << n << " background events)\n";
      continue;
    }
    trainOne(Form("SigmaBDT_%s_A", c.name), "IsOdd<0.5", c.cut, sig, bkg);
    trainOne(Form("SigmaBDT_%s_B", c.name), "IsOdd>0.5", c.cut, sig, bkg);
  }
  std::cout << "Done. Per-class 2-fold weights in SigmaBDT_<class>_<A|B>/.\n";
}
