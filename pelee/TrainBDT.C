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
//  Each fold's events are split again into TMVA's own train/test samples; testFrac is
//  the test share. That test sample only feeds TMVA's overtraining printout -- the scores
//  that matter are the cross-applied ones -- so it is kept at 5%. Until 2026-10-18 it was
//  TMVA's default half, so each fold trained on 25% of all events. Training on 95% of the
//  fold gave AUC 0.896 -> 0.901 and -8% background at 40% signal efficiency, against a
//  seed-to-seed spread of +-0.001 / +-2% (README, "BDT retrained with EXT").
//
//  Output: <prefix>_A/weights/*.xml, <prefix>_B/weights/*.xml
//  Run:    root -l -b -q TrainBDT.C
//          root -l -b -q 'TrainBDT.C("BDTInput.root", 0.5, "SigmaBDT_half", 1)'   # old split, seed 1
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
  "Pi0Mass","Pi0Gammadot","Shower2Energy",
  // Shower Moliere angle -- must stay in the SAME ORDER as BDTV in
  // SigmaSelection.C, which keeps its own copy of this list. The reader there
  // registers variables positionally, so a mismatch silently feeds the BDT the
  // wrong numbers rather than failing.
  "Shower1MoliereAvg",
  // Second-Shower Veto block -- same order as BDTV in SigmaSelection.C.
  // SliceUnassocFrac/Hits excluded -- NoCR artefact, see SigmaSelection.C.
  "SliceClustFrac","SecondShowerNHit","SecondShowerDot","ShrSubClusters",
  "UnassocHitsCalc","UnassocFracCalc",
  "HasPair"
};

// Train one fold.  `trainCut` selects the events used for THIS fold's training;
// the resulting weights are applied (elsewhere) to the complementary events.
void trainFold(const char* dataset, const char* trainCut, TTree* sig, TTree* bkg, double testFrac, int splitSeed)
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
  // nTrain=0 splits the fold in halves; otherwise keep testFrac of it for TMVA's test.
  TString split = "nTrain_Signal=0:nTrain_Background=0";
  if (testFrac < 0.5) {
    const Long64_t nS = sig->GetEntries(trainCut), nB = bkg->GetEntries(trainCut);
    const Long64_t tS = (Long64_t)(testFrac * nS), tB = (Long64_t)(testFrac * nB);
    split = Form("nTrain_Signal=%lld:nTest_Signal=%lld:nTrain_Background=%lld:nTest_Background=%lld",
                 nS - tS, tS, nB - tB, tB);
  }
  loader.PrepareTrainingAndTestTree(trainCut, trainCut, split + Form(":SplitMode=Random:SplitSeed=%d:NormMode=NumEvents:!V", splitSeed));

  factory.BookMethod(&loader, TMVA::Types::kBDT, "BDTG",
    "!H:!V:NTrees=600:MinNodeSize=3%:MaxDepth=3:BoostType=Grad:Shrinkage=0.06:"
    "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=30:NegWeightTreatment=Pray");

  factory.TrainAllMethods();
  factory.TestAllMethods();
  factory.EvaluateAllMethods();
  out->Close();
  std::cout << "  fold '" << dataset << "' done -> " << dataset << "/weights/\n";
}

void TrainBDT(const char* input = "BDTInput.root", double testFrac = 0.05, const char* prefix = "SigmaBDT",
              int splitSeed = 100)   // 100 = TMVA default; other seeds measure the training noise
{
  TMVA::Tools::Instance();
  TFile* in = TFile::Open(input, "READ");
  if (!in || in->IsZombie()) { std::cerr << "cannot open " << input << "\n"; return; }
  TTree* sig = (TTree*)in->Get("TreeS");
  TTree* bkg = (TTree*)in->Get("TreeB");
  std::cout << "Signal=" << sig->GetEntries() << "  Background=" << bkg->GetEntries() << "\n";

  // Fold A: train on even events (IsOdd==0) -> applied to odd.
  trainFold(Form("%s_A", prefix), "IsOdd<0.5", sig, bkg, testFrac, splitSeed);
  // Fold B: train on odd events (IsOdd==1) -> applied to even.
  trainFold(Form("%s_B", prefix), "IsOdd>0.5", sig, bkg, testFrac, splitSeed);

  std::cout << "Done. 2-fold weights: " << prefix << "_A/ (even), " << prefix << "_B/ (odd).\n";
}
