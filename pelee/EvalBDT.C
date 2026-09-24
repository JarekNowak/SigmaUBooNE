// EvalBDT.C -- full efficiency/rejection curve for the trained 49-variable BDT.
// Reuses the two-fold weights from AblateBDT.C; every event is scored by the fold it
// was NOT trained on, so the curve is out-of-sample throughout.
R__LOAD_LIBRARY(libTMVA)
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
static const char* ALL[] = {
 "MuonPID","MuonTrkLength","MuonTrkDistVtx","MuonTrkScore","MuonNHitsY",
 "PionPID","PionTrkLength","PionTrkDistVtx","PionTrkScore","PionNHitsY",
 "ProtonPID","ProtonTrkLength","ProtonTrkDistVtx","ProtonTrkScore","ProtonNHitsY",
 "Shower1NHitsY","Shower1VtxDist","Shower1Energy",
 "DistMuonPion","DistProtonPion","DistProtonMuon",
 "DistShower1Muon","DistShower1Pion","DistShower1Proton",
 "AngleMuonPion","AngleProtonPion","AngleProtonMuon",
 "AngleShower1Muon","AngleShower1Pion","AngleShower1Proton",
 "LambdaMass","SigmaMass","TopologicalScore","NTracks","NShowers",
 "CosPoint","LenSum","PPion","MomAsym","CosLamMu",
 "ct_hits_U","ct_hits_V","ct_hits_Y","ct_vtxoff_U","ct_vtxoff_V","ct_vtxoff_Y",
 "ct_nisl_U","ct_nisl_V","ct_nisl_Y","HasPair",
 "IslNPlanes","IslAny","IslSeeded"};
static const int NV = sizeof(ALL)/sizeof(ALL[0]);
void EvalBDT(){
  TMVA::Tools::Instance();
  TFile* f=TFile::Open("BDTHyp.root");
  TTree* sig=(TTree*)f->Get("TreeS"); TTree* bkg=(TTree*)f->Get("TreeB");
  std::vector<double> sS,bS,sW,bW;
  for(int fold=0; fold<2; ++fold){
    TMVA::Reader reader("!Color:Silent");
    std::vector<float> buf(NV); float spW=0,spO=0;
    for(int i=0;i<NV;++i) reader.AddVariable(ALL[i],&buf[i]);
    reader.AddSpectator("Weight",&spW); reader.AddSpectator("IsOdd",&spO);
    reader.BookMVA("BDTG",TString::Format("abl_full_f%d/weights/A_BDTG.weights.xml",fold));
    for(int isSig=0;isSig<2;++isSig){
      TTree* t=isSig?sig:bkg;
      std::vector<float> val(NV); float w=0,od=0;
      for(int i=0;i<NV;++i) t->SetBranchAddress(ALL[i],&val[i]);
      t->SetBranchAddress("Weight",&w); t->SetBranchAddress("IsOdd",&od);
      for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
        bool inTest = fold==0 ? (od>0.5) : (od<0.5);
        if(!inTest) continue;
        for(int i=0;i<NV;++i) buf[i]=val[i];
        double r=reader.EvaluateMVA("BDTG");
        if(isSig){ sS.push_back(r); sW.push_back(w);} else { bS.push_back(r); bW.push_back(w);} }
      t->ResetBranchAddresses(); }
  }
  double S0=0,B0=0; for(double w:sW)S0+=w; for(double w:bW)B0+=w;
  double base=S0/sqrt(S0+B0);
  printf("\n  ===== BDT working-point curve (49 vars, out-of-sample) =====\n");
  printf("  baseline at NTracksCut: S=%.2f  B=%.1f  pur=%.3f%%  FoM=%.4f\n\n",
         S0,B0,100*S0/(S0+B0),base);
  printf("  %7s %8s %10s %9s %9s %9s %9s %8s\n",
         "cut","S","B","eff_S[%]","eff_B[%]","pur[%]","FoM","xbase");
  // quote at fixed signal-retention points, which is the constraint that matters
  std::vector<double> ss=sS; std::sort(ss.begin(),ss.end());
  for(double keep : {1.00,0.95,0.90,0.85,0.80,0.70,0.60,0.50,0.40,0.30}){
    size_t idx = (size_t)((1.0-keep)*ss.size());
    if(idx>=ss.size()) idx=ss.size()-1;
    double c = keep>=0.999 ? -1e9 : ss[idx];
    double s=0,b=0;
    for(size_t i=0;i<sS.size();++i) if(sS[i]>c) s+=sW[i];
    for(size_t i=0;i<bS.size();++i) if(bS[i]>c) b+=bW[i];
    if(s<=0||b<=0) continue;
    double fom=s/sqrt(s+b);
    printf("  %7.3f %8.2f %10.1f %9.1f %9.1f %9.3f %9.4f %8.2f\n",
           c,s,b,100*s/S0,100*b/B0,100*s/(s+b),fom,fom/base);
  }
}
