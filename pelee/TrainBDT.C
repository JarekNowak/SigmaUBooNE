// ============================================================================
//  TrainBDT.C — 2-fold cross-validation BDTG training for the PeLEE Σ⁰ selection.
//
//  Reads BDTInput.root (TreeS/TreeB at the NTracksCut preselection, dumped by
//  SigmaSelection.C) and trains TWO gradient BDTs on the full candidate feature
//  set (incl. Λ/Σ⁰ mass, shower energy/size, Λ/γ vertex-pointing):
//    fold A → trained on EVEN event numbers (IsOdd==0); applied to ODD events
//    fold B → trained on ODD  event numbers (IsOdd==1); applied to EVEN events
//  Cross-application (in SigmaSelection.C) scores each event with the fold that
//  did NOT see it, giving an unbiased integrated selection.
//
//  Output: SigmaBDT_A/weights/*.xml, SigmaBDT_B/weights/*.xml
//  Run:    root -l -b -q TrainBDT.C
// ============================================================================
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
  "Pi0Mass","Pi0Gammadot","Shower2Energy" };

// Train one fold.  `trainCut` selects the events used for THIS fold's training;
// the resulting weights are applied (elsewhere) to the complementary events.
void trainFold(const char* dataset, const char* trainCut, TTree* sig, TTree* bkg)
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

  // Preselection cut = the CV fold selector (applied to both classes).
  loader.PrepareTrainingAndTestTree(trainCut, trainCut,
    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

  factory.BookMethod(&loader, TMVA::Types::kBDT, "BDTG",
    "!H:!V:NTrees=600:MinNodeSize=3%:MaxDepth=3:BoostType=Grad:Shrinkage=0.06:"
    "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:NegWeightTreatment=Pray");

  factory.TrainAllMethods();
  factory.TestAllMethods();
  factory.EvaluateAllMethods();
  out->Close();
  std::cout << "  fold '" << dataset << "' done -> " << dataset << "/weights/\n";
}

void TrainBDT()
{
  TMVA::Tools::Instance();
  TFile* in = TFile::Open("BDTInput.root", "READ");
  if (!in || in->IsZombie()) { std::cerr << "cannot open BDTInput.root\n"; return; }
  TTree* sig = (TTree*)in->Get("TreeS");
  TTree* bkg = (TTree*)in->Get("TreeB");
  std::cout << "Signal=" << sig->GetEntries() << "  Background=" << bkg->GetEntries() << "\n";

  // Fold A: train on even events (IsOdd==0) -> applied to odd.
  trainFold("SigmaBDT_A", "IsOdd<0.5", sig, bkg);
  // Fold B: train on odd events (IsOdd==1) -> applied to even.
  trainFold("SigmaBDT_B", "IsOdd>0.5", sig, bkg);

  std::cout << "Done. 2-fold weights: SigmaBDT_A/ (even), SigmaBDT_B/ (odd).\n";
}
