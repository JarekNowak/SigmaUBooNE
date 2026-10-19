// RescoreBDT.C -- re-score BDTInput.root with a pair of fold weights (cross-applied by IsOdd:
// odd -> fold A (trained on even), even -> fold B), writing Weight/IsOdd/Sample/BdtScore,
// which CompareBDT.C reads. kVars must stay in the order of TrainBDT.C.
#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"
#include <vector>
#include <string>
#include <cstdio>
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
  "Shower1MoliereAvg",
  "SliceClustFrac","SecondShowerNHit","SecondShowerDot","ShrSubClusters",
  "UnassocHitsCalc","UnassocFracCalc",
  "HasPair"
};
static const int NV = sizeof(kVars) / sizeof(kVars[0]);
void RescoreBDT(const char* in, const char* out, const char* dirA, const char* dirB){
  float v[NV], w = 0, odd = 0, sc = 0, smp = 0;
  TMVA::Reader reader("!Color:Silent");
  for (int i = 0; i < NV; i++) reader.AddVariable(kVars[i], &v[i]);
  reader.AddSpectator("Weight", &w); reader.AddSpectator("IsOdd", &odd);
  reader.BookMVA("A", Form("%s/weights/TMVASigma_BDTG.weights.xml", dirA));
  reader.BookMVA("B", Form("%s/weights/TMVASigma_BDTG.weights.xml", dirB));
  TFile fi(in); TFile fo(out, "RECREATE");
  for (const char* tn : {"TreeS", "TreeB"}) {
    TTree* t = (TTree*)fi.Get(tn);
    for (int i = 0; i < NV; i++) t->SetBranchAddress(kVars[i], &v[i]);
    t->SetBranchAddress("Weight", &w); t->SetBranchAddress("IsOdd", &odd); t->SetBranchAddress("Sample", &smp);
    fo.cd(); TTree o(tn, tn); float nsc = 0;
    o.Branch("Weight", &w); o.Branch("IsOdd", &odd); o.Branch("Sample", &smp); o.Branch("BdtScore", &nsc);
    for (Long64_t e = 0; e < t->GetEntries(); e++) {
      t->GetEntry(e);
      nsc = reader.EvaluateMVA(odd > 0.5 ? "A" : "B");
      o.Fill();
    }
    o.Write();
    printf("  %s: %lld rescored\n", tn, t->GetEntries());
  }
  fo.Close();
}
