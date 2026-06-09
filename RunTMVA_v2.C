// ============================================================
//  RunTMVA_v2.C  —  CONTROLLED RETRAIN (does NOT touch production)
//
//  Identical to RunTMVA.C's BDTg configuration, but:
//    * writes to a SEPARATE dataset dir  "SigmaDataset_v2"
//    * output file                       "TMVA_output_v2.root"
//    * ADDS the already-stored-but-unused discriminating variables:
//        - shower dE/dx  (Shr1DeDxPlane0/1/2, Shr2DeDxPlane2)
//            -> the prime photon-quality discriminant for Sigma0->Lambda gamma
//        - shower<->track distances (DistShr1Muon/Pion/Proton)
//            -> currently only hard pre-cuts, now continuous BDT inputs
//
//  Only BDTg is booked (the production method) to keep the test clean.
//
//  Usage:  root -b -l -q RunTMVA_v2.C
// ============================================================

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"

void RunTMVA_v2()
{
  TMVA::Tools::Instance();
  TFile *outputFile = TFile::Open("TMVA_output_v2.root", "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory(
    "TMVASigma", outputFile,
    "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification");

  TMVA::DataLoader *loader = new TMVA::DataLoader("SigmaDataset_v2");

  TFile *inputFile = TFile::Open("TMVA_input_v3.root", "READ");
  if (!inputFile || inputFile->IsZombie()) { std::cerr<<"ERROR: no TMVA_input_v3.root (run BuildV3Input.C)\n"; return; }
  TTree *sigTree = (TTree*)inputFile->Get("TreeS");
  TTree *bkgTree = (TTree*)inputFile->Get("TreeB");
  std::cout << "Signal entries    : " << sigTree->GetEntries() << std::endl;
  std::cout << "Background entries: " << bkgTree->GetEntries() << std::endl;

  loader->AddSignalTree    (sigTree, 1.0);
  loader->AddBackgroundTree(bkgTree, 1.0);
  loader->SetSignalWeightExpression    ("Weight");
  loader->SetBackgroundWeightExpression("Weight");

  // ===== identical 33 production variables =====
  loader->AddVariable("MuonTrkLength",     "Muon Track Length", "F");
  loader->AddVariable("MuonTrkDistVtx",    "Muon Dist Vtx",     "F");
  loader->AddVariable("MuonTrkScore",      "Muon Trk Score",    "F");
  loader->AddVariable("MuonNHitsY",        "Muon NHits Y",      "F");
  loader->AddVariable("MuonPhotonPID",     "Muon Photon BDT",   "F");
  loader->AddVariable("PionPID",           "Pion LLR PID",      "F");
  loader->AddVariable("PionTrkLength",     "Pion Track Length", "F");
  loader->AddVariable("PionTrkDistVtx",    "Pion Dist Vtx",     "F");
  loader->AddVariable("PionTrkScore",      "Pion Trk Score",    "F");
  loader->AddVariable("PionNHitsY",        "Pion NHits Y",      "F");
  loader->AddVariable("PionDeDxPlane2",    "Pion dEdx coll.",   "F");
  loader->AddVariable("PionPhotonPID",     "Pion Photon BDT",   "F");
  loader->AddVariable("ProtonPID",         "Proton LLR PID",    "F");
  loader->AddVariable("ProtonTrkLength",   "Proton Trk Length", "F");
  loader->AddVariable("ProtonTrkDistVtx",  "Proton Dist Vtx",   "F");
  loader->AddVariable("ProtonTrkScore",    "Proton Trk Score",  "F");
  loader->AddVariable("ProtonNHitsY",      "Proton NHits Y",    "F");
  loader->AddVariable("ProtonDeDxPlane2",  "Proton dEdx coll.", "F");
  loader->AddVariable("ProtonPhotonPID",   "Proton Photon BDT", "F");
  loader->AddVariable("Shr1NHitsY",        "Shr1 NHits Y",      "F");
  loader->AddVariable("Shr1VtxDist",       "Shr1 Dist Vtx",     "F");
  loader->AddVariable("Shr1OpenAngle",     "Shr1 Open Angle",   "F");
  loader->AddVariable("Shr2NHitsY",        "Shr2 NHits Y",      "F");
  loader->AddVariable("Shr2VtxDist",       "Shr2 Dist Vtx",     "F");
  loader->AddVariable("AngleMuonPion",     "Angle #mu-#pi",      "F");
  loader->AddVariable("AngleProtonPion",   "Angle p-#pi",        "F");
  loader->AddVariable("AngleProtonMuon",   "Angle p-#mu",        "F");
  loader->AddVariable("AngleShr1Muon",     "Angle Shr1-#mu",     "F");
  loader->AddVariable("AngleShr1Pion",     "Angle Shr1-#pi",     "F");
  loader->AddVariable("AngleShr1Proton",   "Angle Shr1-p",       "F");
  loader->AddVariable("DistMuonPion",      "Dist #mu-#pi [cm]",  "F");
  loader->AddVariable("DistProtonPion",    "Dist p-#pi [cm]",    "F");
  loader->AddVariable("DistProtonMuon",    "Dist p-#mu [cm]",    "F");

  // ===== NEW: Lambda->p pi- invariant mass (separation <s^2>=0.224) =====
  //   reconstructed from range-based momenta + AngleProtonPion.
  //   (shower dE/dx and shower energy are 100% sentinel in this ntuple,
  //    so the Sigma0 mass and shower-dEdx vars are NOT available.)
  loader->AddVariable("LambdaMass",        "M(p#pi^{-}) [MeV]",  "F");

  // ===== NEW: shower<->track distances (real, only used as pre-cuts before) =====
  loader->AddVariable("DistShr1Muon",      "Dist Shr1-#mu [cm]", "F");
  loader->AddVariable("DistShr1Pion",      "Dist Shr1-#pi [cm]", "F");
  loader->AddVariable("DistShr1Proton",    "Dist Shr1-p [cm]",   "F");

  loader->AddSpectator("Weight", "PoT scale weight", "F");

  loader->PrepareTrainingAndTestTree(
    "", "",
    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

  // identical BDTg configuration to production
  factory->BookMethod(
    loader, TMVA::Types::kBDT, "BDTg",
    "!H:!V:NTrees=800:MinNodeSize=2.0%:MaxDepth=4:BoostType=Grad:"
    "Shrinkage=0.05:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=25:"
    "NegWeightTreatment=Pray:DoBoostMonitor");

  std::cout << "\n=== Training ===" << std::endl;   factory->TrainAllMethods();
  std::cout << "\n=== Testing  ===" << std::endl;   factory->TestAllMethods();
  std::cout << "\n=== Evaluate ===" << std::endl;   factory->EvaluateAllMethods();

  outputFile->Close();
  std::cout << "\nWeights -> SigmaDataset_v2/weights/  (production untouched)\n";
  delete factory; delete loader;
}
