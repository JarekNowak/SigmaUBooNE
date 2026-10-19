// CurveBDT.C -- fine-grained efficiency curve + ROC for the 53-variable BDT.
// Out-of-sample throughout: each event scored by the fold it was not trained on.
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
static const int NV=sizeof(ALL)/sizeof(ALL[0]);
void CurveBDT(){
  TMVA::Tools::Instance();
  TFile* f=TFile::Open("BDTHyp.root");
  TTree* sig=(TTree*)f->Get("TreeS"); TTree* bkg=(TTree*)f->Get("TreeB");
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
        if(isSig){sS.push_back(r);sW.push_back(w);} else {bS.push_back(r);bW.push_back(w);} }
      t->ResetBranchAddresses(); }
  }
  double S0=0,B0=0; for(double w:sW)S0+=w; for(double w:bW)B0+=w;
  double base=S0/sqrt(S0+B0);
  std::vector<double> ss=sS; std::sort(ss.begin(),ss.end());
  printf("\n  ===== BDT efficiency curve, 5%% steps (out-of-sample) =====\n");
  printf("  baseline S=%.2f  B=%.1f  pur=%.3f%%  FoM=%.4f\n\n",S0,B0,100*S0/(S0+B0),base);
  printf("  %8s %8s %9s %9s %9s %8s %8s %8s\n",
         "eff_S[%]","cut","S","B","eff_B[%]","pur[%]","FoM","xbase");
  const int NP=20; std::vector<double> ex,ey,fx,fy,px,py;
  for(int i=NP;i>=1;--i){
    double keep=i/(double)NP;
    size_t idx=(size_t)((1.0-keep)*ss.size()); if(idx>=ss.size()) idx=ss.size()-1;
    double c = (i==NP)? -1e9 : ss[idx];
    double s=0,b=0;
    for(size_t k=0;k<sS.size();++k) if(sS[k]>c) s+=sW[k];
    for(size_t k=0;k<bS.size();++k) if(bS[k]>c) b+=bW[k];
    if(s<=0||b<=0) continue;
    double es=100*s/S0, eb=100*b/B0, pur=100*s/(s+b), fom=s/sqrt(s+b);
    printf("  %8.0f %8.3f %9.2f %9.1f %9.2f %8.3f %8.4f %8.2f\n",es,c,s,b,eb,pur,fom,fom/base);
    ex.push_back(es); ey.push_back(100-eb);
    fx.push_back(es); fy.push_back(fom/base);
    px.push_back(es); py.push_back(pur);
  }
  gStyle->SetOptStat(0);
  TCanvas cv("cv","",1500,500); cv.Divide(3,1);
  cv.cd(1); gPad->SetGrid();
  TGraph gr(ex.size(),&ex[0],&ey[0]); gr.SetTitle("ROC;signal efficiency [%];background rejection [%]");
  gr.SetLineWidth(3); gr.SetLineColor(kAzure+2); gr.SetMarkerStyle(20); gr.Draw("ALP");
  cv.cd(2); gPad->SetGrid();
  TGraph gf(fx.size(),&fx[0],&fy[0]); gf.SetTitle("figure of merit;signal efficiency [%];FoM / baseline");
  gf.SetLineWidth(3); gf.SetLineColor(kRed+1); gf.SetMarkerStyle(20); gf.Draw("ALP");
  TLine l1(0,1,100,1); l1.SetLineStyle(2); l1.Draw();
  cv.cd(3); gPad->SetGrid();
  TGraph gp(px.size(),&px[0],&py[0]); gp.SetTitle("purity;signal efficiency [%];purity [%]");
  gp.SetLineWidth(3); gp.SetLineColor(kGreen+2); gp.SetMarkerStyle(20); gp.Draw("ALP");
  cv.SaveAs("bdt_efficiency_curve.png");
  printf("\n  wrote bdt_efficiency_curve.png\n");
}
