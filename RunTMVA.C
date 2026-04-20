// ============================================================
//  RunTMVA.C
//  TMVA training script for Sigma0 signal/background separation
//  Methods: BDT, kNN, Likelihood
//
//  Usage (ROOT prompt):
//    root -l RunTMVA.C
//  or compiled:
//    root -l -b -q RunTMVA.C
// ============================================================

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void RunTMVA()
{
  // -----------------------------------------------------------
  // Initialise TMVA
  // -----------------------------------------------------------
  TMVA::Tools::Instance();

  // Output file for weights and evaluation results
  TFile *outputFile = TFile::Open("TMVA_output.root", "RECREATE");

  // Factory options:
  //   !V         : suppress verbose output
  //   !Silent    : print progress
  //   Color      : coloured terminal output
  //   DrawProgressBar : show progress bar
  //   Transformations=I : store input variable transformations
  //   AnalysisType=Classification : binary signal/background
  TMVA::Factory *factory = new TMVA::Factory(
    "TMVASigma", outputFile,
    "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification"
  );

  TMVA::DataLoader *loader = new TMVA::DataLoader("SigmaDataset");

  // -----------------------------------------------------------
  // Input file
  // -----------------------------------------------------------
  TFile *inputFile = TFile::Open("TMVA_input.root", "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: cannot open TMVA_input.root" << std::endl;
    return;
  }

  TTree *sigTree = (TTree*)inputFile->Get("TreeS");
  TTree *bkgTree = (TTree*)inputFile->Get("TreeB");

  if (!sigTree || !bkgTree) {
    std::cerr << "ERROR: TreeS or TreeB not found in TMVA_input.root" << std::endl;
    return;
  }

  std::cout << "Signal entries    : " << sigTree->GetEntries() << std::endl;
  std::cout << "Background entries: " << bkgTree->GetEntries() << std::endl;

  // -----------------------------------------------------------
  // Register signal and background trees
  // The "Weight" branch carries the PoT scale factor
  // -----------------------------------------------------------
  loader->AddSignalTree    (sigTree, 1.0);
  loader->AddBackgroundTree(bkgTree, 1.0);

  loader->SetSignalWeightExpression    ("Weight");
  loader->SetBackgroundWeightExpression("Weight");

  // -----------------------------------------------------------
  // Input variables
  // 'F' = float type; all variables saved in Sigma.C
  // -----------------------------------------------------------

  // --- Topological ---
 // loader->AddVariable("NTracks",           "N Tracks",          "F");
 // loader->AddVariable("NShowers",          "N Showers",         "F");

  // --- Muon candidate ---
  loader->AddVariable("MuonPID",           "Muon LLR PID",      "F");
  loader->AddVariable("MuonTrkLength",     "Muon Track Length", "F");
  loader->AddVariable("MuonTrkDistVtx",    "Muon Dist Vtx",     "F");
  loader->AddVariable("MuonTrkScore",      "Muon Trk Score",    "F");
 // loader->AddVariable("MuonDeDxPlane0",    "Muon dEdx Pl0",     "F");
 // loader->AddVariable("MuonDeDxPlane1",    "Muon dEdx Pl1",     "F");
 // loader->AddVariable("MuonDeDxPlane2",    "Muon dEdx Pl2",     "F");
  loader->AddVariable("MuonNHitsY",        "Muon NHits Y",      "F");
  loader->AddVariable("MuonPhotonPID",     "Muon Photon BDT",   "F");

  // --- Pion candidate ---
  loader->AddVariable("PionPID",           "Pion LLR PID",      "F");
  loader->AddVariable("PionTrkLength",     "Pion Track Length", "F");
  loader->AddVariable("PionTrkDistVtx",    "Pion Dist Vtx",     "F");
  loader->AddVariable("PionTrkScore",      "Pion Trk Score",    "F");
  //loader->AddVariable("PionDeDxPlane0",    "Pion dEdx Pl0",     "F");
  //loader->AddVariable("PionDeDxPlane1",    "Pion dEdx Pl1",     "F");
  //loader->AddVariable("PionDeDxPlane2",    "Pion dEdx Pl2",     "F");
  //loader->AddVariable("PionNHitsY",        "Pion NHits Y",      "F");
  loader->AddVariable("PionPhotonPID",     "Pion Photon BDT",   "F");

  // --- Proton candidate ---
  loader->AddVariable("ProtonPID",         "Proton LLR PID",    "F");
  loader->AddVariable("ProtonTrkLength",   "Proton Trk Length", "F");
  loader->AddVariable("ProtonTrkDistVtx",  "Proton Dist Vtx",   "F");
  loader->AddVariable("ProtonTrkScore",    "Proton Trk Score",  "F");
 // loader->AddVariable("ProtonDeDxPlane0",  "Proton dEdx Pl0",   "F");
 // loader->AddVariable("ProtonDeDxPlane1",  "Proton dEdx Pl1",   "F");
 // loader->AddVariable("ProtonDeDxPlane2",  "Proton dEdx Pl2",   "F");
  loader->AddVariable("ProtonNHitsY",      "Proton NHits Y",    "F");
  loader->AddVariable("ProtonPhotonPID",   "Proton Photon BDT", "F");

  // --- Shower 1 ---
  loader->AddVariable("Shr1NHitsY",        "Shr1 NHits Y",      "F");
  loader->AddVariable("Shr1VtxDist",       "Shr1 Dist Vtx",     "F");
  //loader->AddVariable("Shr1DeDxPlane0",    "Shr1 dEdx Pl0",     "F");
 // loader->AddVariable("Shr1DeDxPlane1",    "Shr1 dEdx Pl1",     "F");
 // loader->AddVariable("Shr1DeDxPlane2",    "Shr1 dEdx Pl2",     "F");
  loader->AddVariable("Shr1OpenAngle",     "Shr1 Open Angle",   "F");

  // --- Shower 2 ---
  //loader->AddVariable("Shr2NHitsY",        "Shr2 NHits Y",      "F");
 // loader->AddVariable("Shr2VtxDist",       "Shr2 Dist Vtx",     "F");
 // loader->AddVariable("Shr2DeDxPlane0",    "Shr2 dEdx Pl0",     "F");
 // loader->AddVariable("Shr2DeDxPlane1",    "Shr2 dEdx Pl1",     "F");
 // loader->AddVariable("Shr2DeDxPlane2",    "Shr2 dEdx Pl2",     "F");
 // loader->AddVariable("Shr2OpenAngle",     "Shr2 Open Angle",   "F");

  // --- Opening angles ---
  loader->AddVariable("AngleMuonPion",     "Angle #mu-#pi",      "F");
  loader->AddVariable("AngleProtonPion",   "Angle p-#pi",        "F");
  loader->AddVariable("AngleProtonMuon",   "Angle p-#mu",        "F");
  loader->AddVariable("AngleShr1Muon",     "Angle Shr1-#mu",     "F");
  loader->AddVariable("AngleShr1Pion",     "Angle Shr1-#pi",     "F");
  loader->AddVariable("AngleShr1Proton",   "Angle Shr1-p",       "F");
 // loader->AddVariable("AngleShr2Muon",     "Angle Shr2-#mu",     "F");
 // loader->AddVariable("AngleShr2Pion",     "Angle Shr2-#pi",     "F");
 // loader->AddVariable("AngleShr2Proton",   "Angle Shr2-p",       "F");

  // --- Inter-candidate distances ---
  loader->AddVariable("DistMuonPion",      "Dist #mu-#pi [cm]",  "F");
  loader->AddVariable("DistProtonPion",    "Dist p-#pi [cm]",    "F");
  loader->AddVariable("DistProtonMuon",    "Dist p-#mu [cm]",    "F");
  //loader->AddVariable("DistShr1Muon",      "Dist Shr1-#mu [cm]", "F");
  //loader->AddVariable("DistShr1Pion",      "Dist Shr1-#pi [cm]", "F");
  //loader->AddVariable("DistShr1Proton",    "Dist Shr1-p [cm]",   "F");
  //loader->AddVariable("DistShr2Muon",      "Dist Shr2-#mu [cm]", "F");
  //loader->AddVariable("DistShr2Pion",      "Dist Shr2-#pi [cm]", "F");
 // loader->AddVariable("DistShr2Proton",    "Dist Shr2-p [cm]",   "F");

  // -----------------------------------------------------------
  // Spectator variables (stored but NOT used in training)
  // Useful for cross-checks and cut-and-count after TMVA
  // -----------------------------------------------------------
  loader->AddSpectator("Weight", "PoT scale weight", "F");

  // -----------------------------------------------------------
  // Data preparation
  //   nTrain_Signal      : events used for training (0 = use all)
  //   nTrain_Background  : same
  //   SplitMode=Random   : random train/test split
  //   NormMode=NumEvents : normalise signal and background to
  //                        equal number of events during training
  //   !V                 : suppress verbose
  // -----------------------------------------------------------
  loader->PrepareTrainingAndTestTree(
    "",   // signal cut     (empty = use all)
    "",   // background cut (empty = use all)
    "nTrain_Signal=0:nTrain_Background=0:"
    "SplitMode=Random:NormMode=NumEvents:!V"
  );

  // -----------------------------------------------------------
  // Book methods
  // -----------------------------------------------------------

  // ---- 1. BDT (Gradient Boosted Decision Trees) ----
  // NTrees         : number of trees in the forest
  // MinNodeSize    : minimum fraction of training events in a leaf
  // MaxDepth       : maximum tree depth
  // BoostType=Grad : use gradient boosting (better than AdaBoost for HEP)
  // Shrinkage      : learning rate for gradient boosting
  // UseBaggedBoost : bagging improves generalisation
  // BaggedSampleFraction : fraction of events used per tree
  // nCuts          : grid points for finding optimal split
  // DoBoostMonitor : write boost-monitor histogram
  factory->BookMethod(
    loader,
    TMVA::Types::kBDT,
    "BDT",
    "!H:!V:"
    "NTrees=500:"
    "MinNodeSize=2.5%:"
    "MaxDepth=3:"
    "BoostType=Grad:"
    "Shrinkage=0.10:"
    "UseBaggedBoost:"
    "BaggedSampleFraction=0.5:"
    "nCuts=20:"
    "DoBoostMonitor"
  );

  // ---- 2. BDT with AdaBoost (complementary to Grad) ----
  factory->BookMethod(
    loader,
    TMVA::Types::kBDT,
    "BDTAda",
    "!H:!V:"
    "NTrees=400:"
    "MinNodeSize=2.5%:"
    "MaxDepth=3:"
    "BoostType=AdaBoost:"
    "AdaBoostBeta=0.5:"
    "nCuts=20"
  );

  // ---- 3. k-Nearest Neighbours ----
  // nkNN           : number of neighbours
  // ScaleFrac      : fraction of events used to compute variable scale
  // SigmaFact      : kernel width scaling factor
  // Kernel=Gaus    : Gaussian kernel weighting
  // UseKernel      : enable kernel weighting (smoother boundaries)
  // UseWeight      : use event weights
  // !Trim          : do not trim outlier events
/*  factory->BookMethod(
    loader,
    TMVA::Types::kKNN,
    "KNN",
    "!H:!V:"
    "nkNN=20:"
    "ScaleFrac=0.8:"
    "SigmaFact=1.0:"
    "Kernel=Gaus:"
    "UseKernel=T:"
    "UseWeight=T:"
    "!Trim"
  );
*/
  // ---- 4. Likelihood (Naive Bayes with PDF estimation) ----
  // TransformOutput : apply logit transform to output (maps to [0,1])
  // PDFInterpol=Spline2 : use spline interpolation for PDFs
  // NSmooth         : number of smoothing iterations on PDFs
  // NAvEvtPerBin    : average number of events per PDF bin
  // !H              : suppress help text
  factory->BookMethod(
    loader,
    TMVA::Types::kLikelihood,
    "Likelihood",
    "!H:!V:"
    "TransformOutput:"
    "PDFInterpol=Spline2:"
    "NSmooth=5:"
    "NAvEvtPerBin=50"
  );

  // ---- 5. Likelihood with decorrelation pre-processing ----
  // Decorrelation removes linear correlations between variables
  // before PDF estimation — often improves performance
  factory->BookMethod(
    loader,
    TMVA::Types::kLikelihood,
    "LikelihoodD",
    "!H:!V:"
    "TransformOutput:"
    "PDFInterpol=Spline2:"
    "NSmooth=5:"
    "NAvEvtPerBin=50:"
    "VarTransform=Decorrelate"
  );

  // -----------------------------------------------------------
  // Train, test and evaluate all methods
  // -----------------------------------------------------------
  std::cout << "\n=== Training all methods ===" << std::endl;
  factory->TrainAllMethods();

  std::cout << "\n=== Testing all methods ===" << std::endl;
  factory->TestAllMethods();

  std::cout << "\n=== Evaluating all methods ===" << std::endl;
  factory->EvaluateAllMethods();

  // -----------------------------------------------------------
  // Save and close
  // -----------------------------------------------------------
  outputFile->Close();
  std::cout << "\nResults written to TMVA_output.root" << std::endl;
  std::cout << "Weights written to SigmaDataset/weights/" << std::endl;

  // -----------------------------------------------------------
  // Launch the TMVA GUI (comment out for batch mode)
  // -----------------------------------------------------------
  // TMVA::TMVAGui("TMVA_output.root");

  delete factory;
  delete loader;
}
