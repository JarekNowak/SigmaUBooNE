// pi0_characteristics.C -- what distinguishes the pi0-bearing background from the
// Sigma0 signal, and from the other background classes?
//
// This is the input to a dedicated pi0-rejection BDT: which of the available
// variables actually carry the pi0 signature, as opposed to separating background
// in general. Reported per variable as the TMVA-style separation <S^2> against
// signal, alongside the same quantity for the non-pi0 neutrino class, so a variable
// that is genuinely pi0-specific can be told from one that merely separates any
// background.
#include <vector>
#include <string>
#include <algorithm>

void pi0_characteristics(){
  TFile* f=TFile::Open("BDTInput.root");
  TTree* S=(TTree*)f->Get("TreeS");
  TTree* B=(TTree*)f->Get("TreeB");
  const char* CUT_PI0 = "HasPi0>0.5";
  const char* CUT_NU  = "HasPi0<0.5 && CosmicMu<0.5";
  printf("\n  signal %lld   pi0 background %lld   non-pi0 nu background %lld\n",
         S->GetEntries(), B->GetEntries(CUT_PI0), B->GetEntries(CUT_NU));

  struct V { const char* br; int nb; double lo, hi; };
  std::vector<V> vs = {
    {"SigmaMass",80,1.0,2.0},        {"LambdaMass",80,1.0,1.8},
    {"Pi0Mass",80,-0.05,0.5},        {"Pi0Gammadot",80,-1.1,1.1},
    {"Shower1Energy",80,0,0.8},      {"Shower1NHitsY",80,0,300},
    {"Shower1VtxDist",80,0,80},      {"Shower1DeDxY",80,0,12},
    {"Shower1MoliereAvg",80,0,90},   {"Shower2Energy",80,-0.05,0.5},
    {"SecondShowerNHit",80,0,200},   {"SecondShowerDot",80,-1.1,1.1},
    {"NShowers",8,0,8},              {"NTracks",10,0,10},
    {"PhotonPointAngle",80,0,3.15},  {"PhotonImpact",80,0,40},
    {"LambdaPointAngle",80,0,3.15},  {"LambdaDecayDist",80,0,60},
    {"LambdaImpact",80,0,40},        {"LambdaGammaAngle",80,0,3.15},
    {"SigmaMuonAngle",80,0,3.15},    {"ProtonPID",80,-1,1},
    {"PionPID",80,-1,1},             {"ProtonChi2",80,0,400},
    {"TopologicalScore",80,0,1},     {"UnassocFracCalc",80,0,1.05},
    {"SliceClustFrac",80,0,1.05},    {"ShrSubClusters",80,0,40},
  };
  auto sep=[&](TTree* t1,const char* c1,TTree* t2,const char* c2,const V& v)->double{
    TH1D a("a","",v.nb,v.lo,v.hi), b("b","",v.nb,v.lo,v.hi);
    t1->Draw(Form("%s>>a",v.br),c1,"goff"); t2->Draw(Form("%s>>b",v.br),c2,"goff");
    if(a.Integral()<=0||b.Integral()<=0) return -1;
    a.Scale(1./a.Integral()); b.Scale(1./b.Integral());
    double s=0; for(int i=1;i<=v.nb;++i){ double x=a.GetBinContent(i),y=b.GetBinContent(i);
      if(x+y>0) s += 0.5*(x-y)*(x-y)/(x+y); }
    return s;
  };
  struct R { double sp, sn, d; std::string nm; };
  std::vector<R> res;
  for(auto& v:vs){
    double sp=sep(S,"",B,CUT_PI0,v);   // signal vs pi0 background
    double sn=sep(S,"",B,CUT_NU ,v);   // signal vs non-pi0 nu background
    if(sp<0||sn<0) continue;
    res.push_back({sp,sn,sp-sn,v.br});
  }
  std::sort(res.begin(),res.end(),[](const R&a,const R&b){return a.sp>b.sp;});
  printf("\n  separation from SIGNAL (higher = more useful to a pi0 BDT)\n");
  printf("  %-22s %10s %10s %10s\n","variable","vs pi0","vs nu-no-pi0","difference");
  printf("  %s\n", std::string(56,'-').c_str());
  for(auto& r:res)
    printf("  %-22s %10.4f %10.4f %+10.4f%s\n", r.nm.c_str(), r.sp, r.sn, r.d,
           (r.d>0.02?"   <- pi0-specific":""));
}
