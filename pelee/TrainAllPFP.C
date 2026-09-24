// ============================================================================
//  TrainAllPFP.C -- does a BDT given EVERY reconstructed quantity for EVERY
//  particle in the selected events beat one given the four hand-picked candidate
//  roles?
//
//  Three configurations, all trained and evaluated identically (two-fold CV by
//  event parity, every quoted number out-of-sample):
//    curated  the 57 role-based inputs the selection already uses
//    allpfp   the per-PFP block alone: 8 slots x 15 quantities, ordered by
//             collection-plane hit count, sentinel-padded, plus the multiplicity
//    both     the union
//
//  The BDT is deliberately shallow: with ~3.2k signal events, 121 inputs is
//  ~27 events per variable and overtraining is the main risk.
//
//  Run:  root -l -b -q TrainAllPFP.C
// ============================================================================
R__LOAD_LIBRARY(libTMVA)
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

static const char* CURATED[] = {
  "MuonPID","MuonTrkLength","MuonTrkDistVtx","MuonTrkScore","MuonNHitsY",
  "PionPID","PionTrkLength","PionTrkDistVtx","PionTrkScore","PionNHitsY","PionCaloEnergyY",
  "ProtonPID","ProtonTrkLength","ProtonTrkDistVtx","ProtonTrkScore","ProtonNHitsY","ProtonCaloEnergyY","ProtonChi2",
  "Shower1NHitsY","Shower1VtxDist","Shower1DeDxY","Shower1Energy",
  "DistMuonPion","DistProtonPion","DistProtonMuon","DistShower1Muon","DistShower1Pion","DistShower1Proton",
  "AngleMuonPion","AngleProtonPion","AngleProtonMuon","AngleShower1Muon","AngleShower1Pion","AngleShower1Proton",
  "LambdaMass","SigmaMass","TopologicalScore","NTracks","NShowers",
  "LambdaPointAngle","LambdaImpact","LambdaDecayDist","PhotonPointAngle","PhotonImpact",
  "LambdaGammaAngle","SigmaMuonAngle","Pi0Mass","Pi0Gammadot","Shower2Energy",
  "Shower1MoliereAvg","SliceClustFrac","SecondShowerNHit","SecondShowerDot",
  "ShrSubclusters","UnassocHitsCalc","UnassocFracCalc","HasPair" };
static const int NCUR = sizeof(CURATED)/sizeof(CURATED[0]);
static const char* PFPV[] = {"score","len","llr","chipr","chipi","dvtx","nhity",
                             "caloy","eproton","rmumom","dedxy","shrE","moliere",
                             "cosmu","costh"};
static const int NSLOT=8, NPV=15;

static double aucOf(std::vector<double> a, std::vector<double> b){
  if(a.size()<10||b.size()<10) return 0.5;
  std::sort(b.begin(),b.end()); double t=0;
  for(double x:a) t += (std::lower_bound(b.begin(),b.end(),x)-b.begin())
     + 0.5*(std::upper_bound(b.begin(),b.end(),x)-std::lower_bound(b.begin(),b.end(),x));
  double v=t/(double(a.size())*b.size()); return v<0.5?1-v:v;
}

static void runConfig(const char* tag, const std::vector<std::string>& vars,
                      TTree* sig, TTree* bkg){
  std::vector<double> sS,bS,sW,bW;
  for(int fold=0; fold<2; ++fold){
    TString ds = TString::Format("pfp_%s_f%d",tag,fold);
    gSystem->Exec("rm -rf "+ds);
    TString trainCut = fold==0 ? "IsOdd<0.5" : "IsOdd>0.5";
    TFile* out = TFile::Open(TString::Format("tmvapfp_%s_f%d.root",tag,fold),"RECREATE");
    TMVA::Factory factory("A",out,"!V:Silent:!DrawProgressBar:AnalysisType=Classification");
    TMVA::DataLoader loader(ds.Data());
    for(auto&v:vars) loader.AddVariable(v.c_str(),'F');
    loader.AddSpectator("Weight",'F'); loader.AddSpectator("IsOdd",'F');
    loader.AddSignalTree(sig,1.0); loader.AddBackgroundTree(bkg,1.0);
    loader.SetSignalWeightExpression("Weight");
    loader.SetBackgroundWeightExpression("Weight");
    loader.PrepareTrainingAndTestTree(TCut(trainCut),TCut(trainCut),
      "SplitMode=Random:NormMode=NumEvents:!V");
    factory.BookMethod(&loader,TMVA::Types::kBDT,"BDTG",
      "!H:!V:NTrees=300:MinNodeSize=5%:MaxDepth=3:BoostType=Grad:Shrinkage=0.05:"
      "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:NegWeightTreatment=Pray");
    factory.TrainAllMethods();
    out->Close(); delete out;

    TMVA::Reader reader("!Color:Silent");
    std::vector<float> buf(vars.size()); float spW=0,spO=0;
    for(size_t i=0;i<vars.size();++i) reader.AddVariable(vars[i].c_str(),&buf[i]);
    reader.AddSpectator("Weight",&spW); reader.AddSpectator("IsOdd",&spO);
    reader.BookMVA("BDTG", ds+"/weights/A_BDTG.weights.xml");
    for(int isSig=0;isSig<2;++isSig){
      TTree* t = isSig ? sig : bkg;
      std::vector<float> val(vars.size()); float w=0,od=0;
      for(size_t i=0;i<vars.size();++i) t->SetBranchAddress(vars[i].c_str(),&val[i]);
      t->SetBranchAddress("Weight",&w); t->SetBranchAddress("IsOdd",&od);
      for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
        bool inTest = fold==0 ? (od>0.5) : (od<0.5);
        if(!inTest) continue;
        for(size_t i=0;i<vars.size();++i) buf[i]=val[i];
        double r=reader.EvaluateMVA("BDTG");
        if(isSig){ sS.push_back(r); sW.push_back(w);} else { bS.push_back(r); bW.push_back(w);} }
      t->ResetBranchAddresses();
    }
  }
  double S0=0,B0=0; for(double w:sW)S0+=w; for(double w:bW)B0+=w;
  const double base=S0/std::sqrt(S0+B0);
  auto at=[&](double c,double&s,double&b){ s=0;b=0;
    for(size_t i=0;i<sS.size();++i) if(sS[i]>c) s+=sW[i];
    for(size_t i=0;i<bS.size();++i) if(bS[i]>c) b+=bW[i]; };
  double best=-1,bc=0;
  for(double c=-1;c<=1;c+=0.005){ double s,b; at(c,s,b);
    if(s<=0||b<=0) continue; double f=s/std::sqrt(s+b); if(f>best){best=f;bc=c;} }
  double sp,bp; at(bc,sp,bp);
  printf("  %-9s nvar=%3zu  AUC=%.3f | presel FoM %.4f | peak %+.3f  S=%6.3f B=%8.1f pur=%.3f%%  FoM=%.4f (x%.2f)\n",
         tag,vars.size(),aucOf(sS,bS),base,bc,sp,bp,100*sp/(sp+bp),best,best/base);
  fflush(stdout);
}

void TrainAllPFP(){
  TMVA::Tools::Instance();
  TFile* f=TFile::Open("BDTInput.root");
  TTree* sig=(TTree*)f->Get("TreeS"); TTree* bkg=(TTree*)f->Get("TreeB");
  if(!sig||!bkg){ printf("  no trees\n"); return; }
  if(!sig->GetBranch("p0_score")){ printf("  BDTInput.root has no all-PFP block\n"); return; }
  printf("\n  signal %lld  background %lld   (two-fold CV, all numbers out-of-sample)\n\n",
         sig->GetEntries(),bkg->GetEntries());
  std::vector<std::string> cur, pfp, both;
  for(int i=0;i<NCUR;++i) if(sig->GetBranch(CURATED[i])) cur.push_back(CURATED[i]);
  pfp.push_back("nPFPtot");
  // slot 7 is never occupied (the preselection caps the event at 7 PFPs), so it
  // would feed TMVA 15 identically-constant inputs; use 7 slots.
  for(int k=0;k<NSLOT-1;++k) for(int j=0;j<NPV;++j)
    pfp.push_back(Form("p%d_%s",k,PFPV[j]));
  both=cur; both.insert(both.end(),pfp.begin(),pfp.end());
  runConfig("curated", cur,  sig,bkg);
  runConfig("allpfp",  pfp,  sig,bkg);
  runConfig("both",    both, sig,bkg);
}
