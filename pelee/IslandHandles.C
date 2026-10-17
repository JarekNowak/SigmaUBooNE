// ============================================================================
//  IslandHandles.C -- the Michel tag and the photon conversion gap (PhotonMichel.C)
//  at the Lambda-style operating points, on the full-statistics island caches.
//
//  Both handles were measured on the Run 3 list alone (2026-09-12), where the
//  two-plane point had 8-10 background events and showed nothing. Here the
//  background is split into that list and the RHC extension, which played no part
//  in defining the handles, so the extension column is out of sample:
//    B_ext = (B_full * POT_full - B_run3 * POT_run3) / POT_ext,  raw_ext = raw_full - raw_run3.
//  Signal: GENIE (normalisation without hyperon FSI) and NuWro (with), see the README.
//
//  Run:  root -l -b -q 'IslandHandles.C("isl_sig_genie.root","isl_sig_nuwro.root","isl_bkg.root","isl_bkg_run3.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>

struct HE { int hasPair, islN, muContained; double w, mass, len, cosA, phDist, michelDist; };

static std::vector<HE> loadH(const char* fn){
  std::vector<HE> v; TFile f(fn); TTree* t=(TTree*)f.Get("isl"); HE e;
  t->SetBranchAddress("hasPair",&e.hasPair); t->SetBranchAddress("islNPlanes",&e.islN);
  t->SetBranchAddress("muContained",&e.muContained); t->SetBranchAddress("weight",&e.w);
  t->SetBranchAddress("mass",&e.mass); t->SetBranchAddress("decayLen",&e.len);
  t->SetBranchAddress("cosAlpha",&e.cosA); t->SetBranchAddress("phDist",&e.phDist);
  t->SetBranchAddress("michelDist",&e.michelDist);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); v.push_back(e); }
  return v;
}

static double fomH(double s,double b,double f=0.15){ return (s+b>0)? s/std::sqrt(s+b+f*f*b*b) : 0; }

void IslandHandles(const char* sigGenie="isl_sig_genie.root", const char* sigNuwro="isl_sig_nuwro.root",
                   const char* bkgFull="isl_bkg.root", const char* bkgRun3="isl_bkg_run3.root",
                   double potFull=7.0488e21, double potRun3=2.7795e21)
{
  const auto SG=loadH(sigGenie), SN=loadH(sigNuwro), BF=loadH(bkgFull), B3=loadH(bkgRun3);
  const double potExt=potFull-potRun3;
  typedef std::function<bool(const HE&)> Cut;
  const Cut none =[](const HE&){ return true; };
  const Cut michel=[](const HE& e){ return !e.muContained || (e.michelDist>=0 && e.michelDist<5); };
  const Cut gap  =[](const HE& e){ return e.phDist>2; };
  const Cut both =[&](const HE& e){ return michel(e) && gap(e); };
  const char* name[4]={"no handle","Michel tag","photon gap > 2 cm","Michel + gap"};
  const Cut* cut[4]={&none,&michel,&gap,&both};

  for(int planes=2; planes>=1; --planes){
    auto op=[&](const HE& e){ return e.hasPair && std::fabs(e.mass-1.115683)<0.020 && e.len>0.5 && e.cosA>0.95 && e.islN>=planes; };
    auto sum=[&](const std::vector<HE>& v, const Cut& c, long& n){ double s=0; n=0; for(auto&e:v) if(op(e)&&c(e)){ s+=e.w; n++; } return s; };
    printf("\n  >= %d plane%s   %-17s %14s %13s | %13s %13s %13s | FoM15 run3 / ext / full\n", planes, planes>1?"s":" ",
           "", "S GENIE", "S NuWro", "B run3", "B ext", "B full");
    for(int k=0;k<4;++k){
      long nG,nN,nF,n3; const double sG=sum(SG,*cut[k],nG), sN=sum(SN,*cut[k],nN), bF=sum(BF,*cut[k],nF), b3=sum(B3,*cut[k],n3);
      const double bX=(bF*potFull-b3*potRun3)/potExt; const long nX=nF-n3;
      printf("  %-30s %6.3f (%3ld) %6.3f (%3ld) | %6.2f (%3ld) %6.2f (%3ld) %6.2f (%3ld) | %.3f / %.3f / %.3f\n",
             name[k], sG,nG, sN,nN, b3,n3, bX,nX, bF,nF, fomH(sG,b3), fomH(sG,bX), fomH(sG,bF));
    }
  }
  fflush(stdout);
}
