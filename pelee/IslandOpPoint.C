// ============================================================================
//  IslandOpPoint.C -- S and B, with their MC-statistical errors, at the
//  Lambda-style operating points, read from IslandProcess.C caches.
//
//  Operating point as in SigmaLambdaSel.C: best pair (no precut veto), Lambda mass
//  within +-20 MeV, decay length > 0.5 cm, cos(alpha) > 0.95, island on >= N planes.
//  The error on each yield is sqrt(sum w^2), i.e. the MC statistics of that sample.
//  The FoM is S/sqrt(S+B+(fB)^2) for the background systematics f quoted in the README.
//
//  Run:  root -l -b -q 'IslandOpPoint.C("isl_sig_nuwro.root","isl_bkg.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>

struct Yield { long raw=0; double w=0, w2=0; };

static void count(const char* fn, int planes, double win, double cmin, Yield& y){
  TFile f(fn); TTree* t=(TTree*)f.Get("isl");
  if(!t){ printf("  no isl tree in %s\n",fn); return; }
  int hasPair=0, islN=0; double w=0, mass=-1, len=-1, cosA=-2;
  t->SetBranchAddress("hasPair",&hasPair); t->SetBranchAddress("islNPlanes",&islN);
  t->SetBranchAddress("weight",&w);        t->SetBranchAddress("mass",&mass);
  t->SetBranchAddress("decayLen",&len);    t->SetBranchAddress("cosAlpha",&cosA);
  for(Long64_t i=0;i<t->GetEntries();++i){
    t->GetEntry(i);
    if(!hasPair || std::fabs(mass-1.115683)>=win || len<=0.5 || cosA<=cmin || islN<planes) continue;
    y.raw++; y.w+=w; y.w2+=w*w;
  }
}

static double fom(double s,double b,double f){ return (s+b>0)? s/std::sqrt(s+b+f*f*b*b) : 0; }

void IslandOpPoint(const char* sigFile="isl_sig_nuwro.root", const char* bkgFile="isl_bkg.root",
                   double win=0.020, double cmin=0.95)
{
  printf("\n  signal %s\n  background %s\n  mass +-%.0f MeV, cos(alpha) > %.2f, decay length > 0.5 cm\n",
         sigFile, bkgFile, 1e3*win, cmin);
  printf("\n  %-8s %16s %18s %8s %8s   FoM f=0 / 15%% / 30%%\n", "planes","S (raw)","B (raw)","purity","dB/B");
  for(int planes=1; planes<=2; ++planes){
    Yield s, b; count(sigFile,planes,win,cmin,s); count(bkgFile,planes,win,cmin,b);
    printf("  >= %-5d %6.3f+-%.3f (%3ld) %7.2f+-%.2f (%4ld) %7.1f%% %7.1f%%   %.3f / %.3f / %.3f\n",
           planes, s.w, std::sqrt(s.w2), s.raw, b.w, std::sqrt(b.w2), b.raw,
           (s.w+b.w>0)?100*s.w/(s.w+b.w):0, b.w>0?100*std::sqrt(b.w2)/b.w:0,
           fom(s.w,b.w,0), fom(s.w,b.w,0.15), fom(s.w,b.w,0.30));
  }
}
