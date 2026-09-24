// ============================================================================
//  IslandSideband.C -- sidebands for the Lambda-style selection, on the cache
//  written by IslandProcess.C. Their recipe (Hyperon technote): cut on one of
//  {invariant mass, alpha} and look at the other, signal bin blinded. Added here:
//  the island-fail region, which is the largest background-only sample.
//
//  For each region: POT-scaled signal contamination, background, raw counts, and
//  whether the background COMPOSITION (mode, candidate truth) matches the signal
//  region -- the property that makes a sideband usable for a normalisation
//  constraint.
//
//  Run:  root -l -b -q 'IslandSideband.C("isl_sig_genie.root","isl_bkg.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <vector>

static const double ML=1.115683;
struct E { int isSig, mode, hasPair, islN, prTrue, piTrue, muTrue; double w, mass, len, cosA; };
static std::vector<E> load(const char* fn){ std::vector<E> v; TFile f(fn); TTree* t=(TTree*)f.Get("isl"); E e;
  t->SetBranchAddress("isSig",&e.isSig); t->SetBranchAddress("mode",&e.mode); t->SetBranchAddress("hasPair",&e.hasPair);
  t->SetBranchAddress("islNPlanes",&e.islN); t->SetBranchAddress("prTrue",&e.prTrue); t->SetBranchAddress("piTrue",&e.piTrue);
  t->SetBranchAddress("muTrue",&e.muTrue); t->SetBranchAddress("weight",&e.w); t->SetBranchAddress("mass",&e.mass);
  t->SetBranchAddress("decayLen",&e.len); t->SetBranchAddress("cosAlpha",&e.cosA);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); v.push_back(e); } return v; }

struct Region { const char* name; bool (*f)(const E&); };
static bool base(const E& e){ return e.hasPair && e.len>0.5; }
static bool inMass(const E& e){ return std::fabs(e.mass-ML)<0.020; }
static bool sigReg (const E& e){ return base(e) && inMass(e) && e.cosA>0.95 && e.islN>=1; }
static bool massSB (const E& e){ return base(e) && e.mass>0 && std::fabs(e.mass-ML)>0.040 && std::fabs(e.mass-ML)<0.200 && e.cosA>0.95 && e.islN>=1; }
static bool alphaSB(const E& e){ return base(e) && inMass(e) && e.cosA>0.0 && e.cosA<0.80 && e.islN>=1; }
static bool islSB  (const E& e){ return base(e) && inMass(e) && e.cosA>0.95 && e.islN==0; }
static bool loose  (const E& e){ return base(e) && inMass(e) && e.cosA>0.95; }          // island not applied
static bool massSBnoIsl(const E& e){ return base(e) && e.mass>0 && std::fabs(e.mass-ML)>0.040 && std::fabs(e.mass-ML)<0.200 && e.cosA>0.95; }

void IslandSideband(const char* sigFile="isl_sig_genie.root", const char* bkgFile="isl_bkg.root")
{
  auto S=load(sigFile), B=load(bkgFile);
  Region R[]={{"SIGNAL REGION (blind)",sigReg},{"mass sideband 40-200 MeV",massSB},{"alpha sideband 0<cos<0.8",alphaSB},
              {"island-fail region",islSB},{"mass sideband, no island",massSBnoIsl},{"signal region, no island",loose}};
  const char* mn[]={"QE","RES","DIS","COH","MEC","Oth"}; auto mi=[](int m){ return m==0?0:m==1?1:m==2?2:m==3?3:m==10?4:5; };
  printf("\n  %-28s %6s %7s %5s %5s %7s | %-28s | %s\n","region","S","B","rawS","rawB","S/B[%]","bkg mode QE/RES/DIS/MEC [%]","bkg: p true / pi- true / mu- [%]");
  for(auto& r:R){ double s=0,b=0, bm[6]={0}, tp=0,tpi=0,tmu=0; long ns=0,nb=0;
    for(auto&e:S) if(r.f(e)){ s+=e.w; ns++; }
    for(auto&e:B) if(r.f(e)){ b+=e.w; nb++; bm[mi(e.mode)]+=e.w; if(e.prTrue==2212)tp+=e.w; if(e.piTrue==-211)tpi+=e.w; if(e.muTrue==13)tmu+=e.w; }
    printf("  %-28s %6.2f %7.1f %5ld %5ld %7.1f | %5.0f %5.0f %5.0f %5.0f        | %5.0f %5.0f %5.0f\n", r.name, s, b, ns, nb, b>0?100*s/b:0,
           b>0?100*bm[0]/b:0, b>0?100*bm[1]/b:0, b>0?100*bm[2]/b:0, b>0?100*bm[4]/b:0, b>0?100*tp/b:0, b>0?100*tpi/b:0, b>0?100*tmu/b:0); }
  // mass spectrum in the alpha-passing, island-passing sample: signal vs background per 20 MeV bin (blind bin marked)
  printf("\n  mass spectrum with cos alpha>0.95 & island>=1 (20 MeV bins; * = blinded signal bin)\n  %-14s %7s %7s %6s\n","m-mL [MeV]","S","B","rawB");
  for(int k=-8;k<10;++k){ const double lo=0.020*k, hi=lo+0.020; double s=0,b=0; long nb=0;
    for(auto&e:S) if(base(e)&&e.mass>0&&e.mass-ML>=lo&&e.mass-ML<hi&&e.cosA>0.95&&e.islN>=1) s+=e.w;
    for(auto&e:B) if(base(e)&&e.mass>0&&e.mass-ML>=lo&&e.mass-ML<hi&&e.cosA>0.95&&e.islN>=1){ b+=e.w; nb++; }
    printf("  %+4.0f to %+4.0f %s %7.2f %7.1f %6ld\n", 1000*lo, 1000*hi, (k==-1||k==0)?"*":" ", s, b, nb); }
  // alpha spectrum in the mass window, island passing
  printf("\n  cos alpha spectrum in the mass window & island>=1 (* = signal bin)\n  %-14s %7s %7s %6s\n","cos alpha","S","B","rawB");
  const double ce[]={-1,0,0.5,0.8,0.9,0.95,0.98,1.0001};
  for(int k=0;k<7;++k){ double s=0,b=0; long nb=0;
    for(auto&e:S) if(base(e)&&inMass(e)&&e.cosA>=ce[k]&&e.cosA<ce[k+1]&&e.islN>=1) s+=e.w;
    for(auto&e:B) if(base(e)&&inMass(e)&&e.cosA>=ce[k]&&e.cosA<ce[k+1]&&e.islN>=1){ b+=e.w; nb++; }
    printf("  %5.2f-%-5.2f %s %7.2f %7.1f %6ld\n", ce[k], ce[k+1]>1?1.0:ce[k+1], ce[k]>=0.95?"*":" ", s, b, nb); }
  fflush(stdout);
}
