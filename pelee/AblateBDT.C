// AblateBDT.C -- variable-block ablation for the Sigma0 BDT.
//
// Two-fold cross-validation by event parity: train on one fold, evaluate ONLY on the
// other, so every quoted number is out-of-sample. With 3043 signal events against 53
// variables (~1450 training events per fold) overtraining is still the main risk, so
// the BDT is deliberately shallow and regularised.
//
// Reported per configuration: out-of-sample AUC, and the POT-weighted figure of merit
// at a working point holding 90% of the signal -- the retention constraint, rather
// than the FoM-optimal cut which tends to sit at much lower efficiency.
//
// Each fold is split again into TMVA's own train/test samples. That test sample only
// feeds TMVA's printout -- the numbers here come from the held-out fold -- so testFrac
// keeps 5% of the fold for it. Until 2026-10-19 it was TMVA's default half, so each fold
// trained on 25% of all events (the same fix as TrainBDT.C on 2026-10-18).
//
// Run:  root -l -b -q AblateBDT.C
//       root -l -b -q 'AblateBDT.C("BDTHyp.root", 0.5)'      # old split
//       root -l -b -q 'AblateBDT.C("BDTHyp.root", 0.05, 7)'  # another split seed
R__LOAD_LIBRARY(libTMVA)
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
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
static const int NALL = sizeof(ALL)/sizeof(ALL[0]);

static bool isCT(const std::string& v){ return v.rfind("ct_",0)==0; }
static bool isShower(const std::string& v){ return v.find("Shower")!=std::string::npos; }
static bool isCore(const std::string& v){
  return v=="PionPID"||v=="MuonPID"||v=="CosPoint"||v=="LenSum"||v=="PPion"
       ||v=="MomAsym"||v=="CosLamMu"||v=="LambdaMass"; }

static double aucOf(std::vector<double> a,std::vector<double> b){
  if(a.size()<10||b.size()<10) return 0.5;
  std::sort(b.begin(),b.end()); double t=0;
  for(double x:a) t += (std::lower_bound(b.begin(),b.end(),x)-b.begin())
     + 0.5*(std::upper_bound(b.begin(),b.end(),x)-std::lower_bound(b.begin(),b.end(),x));
  double v=t/(double(a.size())*b.size()); return v<0.5?1-v:v;
}

// TMVA split options for one fold: nTrain=0 halves it; otherwise keep testFrac of it.
static TString splitOpt(TTree* sig, TTree* bkg, const char* cut, double testFrac, int seed){
  TString s="nTrain_Signal=0:nTrain_Background=0";
  if(testFrac<0.5){
    const Long64_t nS=sig->GetEntries(cut), nB=bkg->GetEntries(cut);
    const Long64_t tS=(Long64_t)(testFrac*nS), tB=(Long64_t)(testFrac*nB);
    s=Form("nTrain_Signal=%lld:nTest_Signal=%lld:nTrain_Background=%lld:nTest_Background=%lld",
           nS-tS,tS,nB-tB,tB); }
  return s+Form(":SplitMode=Random:SplitSeed=%d:NormMode=NumEvents:!V",seed);
}

static void runConfig(const char* tag, const std::vector<std::string>& vars,
                      TTree* sig, TTree* bkg, double testFrac, int seed){
  std::vector<double> sScore,bScore,sW,bW;
  for(int fold=0; fold<2; ++fold){
    TString ds = TString::Format("abl_%s_f%d",tag,fold);
    gSystem->Exec("rm -rf "+ds);
    TString trainCut = fold==0 ? "IsOdd<0.5" : "IsOdd>0.5";
    TString testCut  = fold==0 ? "IsOdd>0.5" : "IsOdd<0.5";
    TFile* out = TFile::Open(TString::Format("tmva_%s_f%d.root",tag,fold),"RECREATE");
    TMVA::Factory factory("A",out,"!V:Silent:!DrawProgressBar:AnalysisType=Classification");
    TMVA::DataLoader loader(ds.Data());
    for(auto&v:vars) loader.AddVariable(v.c_str(),'F');
    loader.AddSpectator("Weight",'F');
    loader.AddSpectator("IsOdd",'F');
    loader.AddSignalTree(sig,1.0);
    loader.AddBackgroundTree(bkg,1.0);
    loader.SetSignalWeightExpression("Weight");
    loader.SetBackgroundWeightExpression("Weight");
    loader.PrepareTrainingAndTestTree(TCut(trainCut),TCut(trainCut),
      splitOpt(sig,bkg,trainCut,testFrac,seed));
    // shallow + regularised: ~1.1k signal events cannot support a deep forest
    factory.BookMethod(&loader,TMVA::Types::kBDT,"BDTG",
      "!H:!V:NTrees=300:MinNodeSize=5%:MaxDepth=3:BoostType=Grad:Shrinkage=0.05:"
      "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:NegWeightTreatment=Pray");
    factory.TrainAllMethods();
    out->Close(); delete out;

    // evaluate on the HELD-OUT fold
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
        if(isSig){ sScore.push_back(r); sW.push_back(w); }
        else     { bScore.push_back(r); bW.push_back(w); } }
      t->ResetBranchAddresses();
    }
  }
  double S0=0,B0=0; for(double w:sW)S0+=w; for(double w:bW)B0+=w;
  double base=S0/sqrt(S0+B0);
  // working point at 90% signal retention
  std::vector<double> ss=sScore; std::sort(ss.begin(),ss.end());
  double c90=ss[(size_t)(0.10*ss.size())];
  double s90=0,b90=0;
  for(size_t i=0;i<sScore.size();++i) if(sScore[i]>c90) s90+=sW[i];
  for(size_t i=0;i<bScore.size();++i) if(bScore[i]>c90) b90+=bW[i];
  // best FoM over a scan
  double best=base,bestc=-99;
  for(double c=-1.0;c<=1.0001;c+=0.05){
    double s=0,b=0;
    for(size_t i=0;i<sScore.size();++i) if(sScore[i]>c) s+=sW[i];
    for(size_t i=0;i<bScore.size();++i) if(bScore[i]>c) b+=bW[i];
    if(s<=0||b<=0) continue;
    double f=s/sqrt(s+b); if(f>best){best=f;bestc=c;} }
  printf("  %-22s nvar=%2zu  AUC=%.3f | base FoM %.4f | best %.4f (x%.2f) | @90%%: B kept %.1f%%  FoM %.4f (x%.2f)\n",
    tag,vars.size(),aucOf(sScore,bScore),base,best,best/base,
    B0>0?100.*b90/B0:0., s90/sqrt(s90+b90), (s90/sqrt(s90+b90))/base);
  fflush(stdout);
}

void AblateBDT(const char* input="BDTHyp.root", double testFrac=0.05,
               int seed=100){   // 100 = TMVA default; other seeds measure the training noise
  TMVA::Tools::Instance();
  TFile* f=TFile::Open(input);
  TTree* sig=(TTree*)f->Get("TreeS"); TTree* bkg=(TTree*)f->Get("TreeB");
  printf("\n  signal %lld  background %lld\n",sig->GetEntries(),bkg->GetEntries());
  printf("  two-fold CV, all numbers OUT-OF-SAMPLE; TMVA test share %.2f, seed %d\n\n",testFrac,seed);
  std::vector<std::string> all,noCT,noCTnoShw,ctOnly,core;
  for(int i=0;i<NALL;++i){ std::string v=ALL[i];
    all.push_back(v);
    if(!isCT(v)) noCT.push_back(v);
    if(!isCT(v)&&!isShower(v)) noCTnoShw.push_back(v);
    if(isCT(v)) ctOnly.push_back(v);
    if(isCore(v)) core.push_back(v); }
  runConfig("full",        all,       sig,bkg,testFrac,seed);
  runConfig("noConnected", noCT,      sig,bkg,testFrac,seed);
  runConfig("noCT_noShow", noCTnoShw, sig,bkg,testFrac,seed);
  runConfig("connectedOnly",ctOnly,   sig,bkg,testFrac,seed);
  runConfig("core8",       core,      sig,bkg,testFrac,seed);
}
