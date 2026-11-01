// PeakBDT.C -- fine scan of the hypana BDTG cut to locate the figure-of-merit
// peak, out-of-sample (each event scored by the fold it was not trained on).
// Prints the peak and a neighbourhood so the plateau width is visible: a peak
// sitting on a knife edge is not a working point you want to adopt.
//
// Optional per-event weights (2026-10-31): BDTHyp.root carries no event identifiers, but its
// TreeS / TreeB are filled in the same order as the island caches isl_sig_genie_qe.root /
// isl_bkg.root (checked: identical (nTracks, nShowers) at all 47,033 / 123,000 indices), which
// do. sigW / bkgW are "run sub evt w" tables (HypSignalWeights.C: GENIE hyperon cross section x
// Sigma0 polarisation; for the background, the GENIE hyperon-QE events with their cross-section
// weight); events not listed keep w = 1.
//   root -l -b -q 'PeakBDT.C("hypw_genie_qe.tsv","bkg_hypqe.tsv")'
R__LOAD_LIBRARY(libTMVA)
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "HypConfig.h"
#include "FluxWeights.h"
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
static const int NV=sizeof(ALL)/sizeof(ALL[0]);
// per-index extra weight from the island cache's (run, sub, evt) and a "run sub evt w" table
static std::vector<double> indexWeights(const char* cache, const char* table, Long64_t n){
  std::vector<double> w(n, 1.0);
  if(!table||!*table) return w;
  const auto m = loadFluxWeights(table);
  TFile* c=TFile::Open(cache); TTree* t=c?(TTree*)c->Get("isl"):nullptr;
  if(!t||t->GetEntries()!=n){ printf("  %s does not match the BDT tree (%lld entries): weights ignored\n",cache,n); return w; }
  t->SetEstimate(n+1); t->Draw("run:sub:evt","","goff para");
  long hit=0; double sum=0;
  for(Long64_t i=0;i<n;++i){ auto it=m.find(rseKey((unsigned)t->GetVal(0)[i],(unsigned)t->GetVal(1)[i],(unsigned)t->GetVal(2)[i]));
    if(it!=m.end()){ w[i]=it->second; hit++; sum+=it->second; } }
  printf("  %s: %ld of %lld events weighted (<w> %.3f)\n",table,hit,n,hit?sum/hit:0.);
  delete c; return w;
}
void PeakBDT(const char* sigW="", const char* bkgW=""){
  TMVA::Tools::Instance();
  TFile* f=TFile::Open("BDTHyp.root");
  TTree* sig=(TTree*)f->Get("TreeS"); TTree* bkg=(TTree*)f->Get("TreeB");
  const auto xS = indexWeights("isl_sig_genie_qe.root", sigW, sig->GetEntries());
  const auto xB = indexWeights("isl_bkg.root",          bkgW, bkg->GetEntries());
  std::vector<double> sS,bS,sW,bW;
  for(int fold=0;fold<2;++fold){
    TMVA::Reader rd("!Color:Silent");
    std::vector<float> buf(NV); float spW=0,spO=0;
    for(int i=0;i<NV;++i) rd.AddVariable(ALL[i],&buf[i]);
    rd.AddSpectator("Weight",&spW); rd.AddSpectator("IsOdd",&spO);
    rd.BookMVA("BDTG",TString::Format("abl_full_f%d/weights/A_BDTG.weights.xml",fold));
    for(int isSig=0;isSig<2;++isSig){
      TTree* t=isSig?sig:bkg; std::vector<float> val(NV); float w=0,od=0;
      for(int i=0;i<NV;++i) t->SetBranchAddress(ALL[i],&val[i]);
      t->SetBranchAddress("Weight",&w); t->SetBranchAddress("IsOdd",&od);
      for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
        bool inTest = fold==0 ? (od>0.5):(od<0.5); if(!inTest) continue;
        for(int i=0;i<NV;++i) buf[i]=val[i];
        double r=rd.EvaluateMVA("BDTG");
        if(isSig){sS.push_back(r);sW.push_back(w*xS[e]);} else {bS.push_back(r);bW.push_back(w*xB[e]);} }
      t->ResetBranchAddresses(); }
  }
  double S0=0,B0=0; for(double w:sW)S0+=w; for(double w:bW)B0+=w;
  const double base=S0/std::sqrt(S0+B0);
  printf("\n  baseline (no BDT cut): S=%.2f B=%.1f pur=%.3f%% FoM=%.4f\n",S0,B0,100*S0/(S0+B0),base);
  auto at=[&](double c,double&s,double&b){ s=0;b=0;
    for(size_t i=0;i<sS.size();++i) if(sS[i]>c) s+=sW[i];
    for(size_t i=0;i<bS.size();++i) if(bS[i]>c) b+=bW[i]; };
  double best=-1,bc=0;
  for(double c=-1.0;c<=1.0;c+=0.002){ double s,b; at(c,s,b);
    if(s<=0||b<=0) continue; double fm=s/std::sqrt(s+b); if(fm>best){best=fm;bc=c;} }
  printf("\n  === FoM peak at BDTG > %+.3f ===\n",bc);
  printf("  %8s %9s %10s %9s %9s %9s %8s\n","cut","S","B","eff_S[%]","eff_B[%]","pur[%]","FoM");
  for(double d=-0.30; d<=0.3001; d+=0.05){
    double c=bc+d, s,b; at(c,s,b); if(s<=0||b<=0) continue;
    printf("  %8.3f %9.2f %10.1f %9.1f %9.2f %9.3f %8.4f%s\n",
           c,s,b,100*s/S0,100*b/B0,100*s/(s+b),s/std::sqrt(s+b),
           std::fabs(d)<1e-6?"   <-- peak":"");
  }
  double s,b; at(bc,s,b);
  printf("\n  peak: BDTG > %+.3f   S=%.2f  B=%.1f  eff=%.1f%%  pur=%.3f%%  FoM=%.4f  (x%.2f baseline)\n",
         bc,s,b,100*s/S0,100*s/(s+b),best,best/base);
  double sc,bc2; at(HypCfg::BDTG_CUT,sc,bc2);
  printf("  configured operating point (HypConfig.h) BDTG > %+.3f:\n",HypCfg::BDTG_CUT);
  printf("    S=%.2f  B=%.1f  eff=%.1f%%  pur=%.3f%%  FoM=%.4f  (x%.2f baseline)\n",
         sc,bc2,100*sc/S0,100*sc/(sc+bc2),sc/std::sqrt(sc+bc2),(sc/std::sqrt(sc+bc2))/base);
}
