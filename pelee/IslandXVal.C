// ============================================================================
//  IslandXVal.C -- odd/even cross-validation of the Lambda-style operating point.
//
//  The thresholds in SigmaLambdaSel.C were chosen by scanning S and B on the full
//  cached samples (IslandAna.C), so the FoM quoted there is optimistically biased.
//  Here the scan is repeated on even-numbered events only, the chosen point is
//  evaluated on the odd events, and vice versa; the two out-of-sample halves add
//  to an unbiased full-sample estimate. Both folds must also pick compatible
//  thresholds, otherwise the "optimum" is noise.
//
//  Run:  root -l -b -q 'IslandXVal.C("isl_sig_genie.root","isl_bkg.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <vector>

struct E { int evt, hasPair, precutOK, islN; double w, mass, decayLen, cosA; };
static std::vector<E> load(const char* fn){
  std::vector<E> v; TFile f(fn); TTree* t=(TTree*)f.Get("isl"); E e;
  t->SetBranchAddress("evt",&e.evt); t->SetBranchAddress("hasPair",&e.hasPair); t->SetBranchAddress("precutOK",&e.precutOK);
  t->SetBranchAddress("islNPlanes",&e.islN); t->SetBranchAddress("weight",&e.w); t->SetBranchAddress("mass",&e.mass);
  t->SetBranchAddress("decayLen",&e.decayLen); t->SetBranchAddress("cosAlpha",&e.cosA);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); v.push_back(e); } return v; }

struct P { double win, cmin; int planes, precut; };
static bool pass(const E& e, const P& p){
  if(p.precut ? !e.precutOK : !e.hasPair) return false;
  return std::fabs(e.mass-1.115683)<p.win && e.decayLen>0.5 && e.cosA>p.cmin && e.islN>=p.planes; }
static double fom(double s,double b,double f){ return (s+b>0)? s/std::sqrt(s+b+f*f*b*b) : 0; }

void IslandXVal(const char* sigFile="isl_sig_genie.root", const char* bkgFile="isl_bkg.root", double fsys=0.15)
{
  auto S=load(sigFile), B=load(bkgFile);
  std::vector<P> grid;
  for(double w : {0.020,0.030,0.040,0.060}) for(double c : {0.80,0.90,0.95,0.98})
    for(int pl : {1,2}) for(int pc : {0,1}) grid.push_back({w,c,pl,pc});

  auto count=[&](const P& p, int parity, double& s, double& b, long& ns, long& nb){
    s=b=0; ns=nb=0;
    for(auto&e:S) if((e.evt%2)==parity && pass(e,p)){ s+=e.w; ns++; }
    for(auto&e:B) if((e.evt%2)==parity && pass(e,p)){ b+=e.w; nb++; } };

  printf("\n  cross-validation, objective S/sqrt(S+B+(%.0f%% B)^2), grid of %zu points\n", 100*fsys, grid.size());
  double Soos=0, Boos=0; long nSoos=0, nBoos=0;
  for(int train=0; train<2; ++train){
    const int test=1-train; P best{}; double bestF=-1;
    for(auto& p:grid){ double s,b; long ns,nb; count(p,train,s,b,ns,nb); const double f=fom(s,b,fsys); if(f>bestF){ bestF=f; best=p; } }
    double s,b; long ns,nb; count(best,test,s,b,ns,nb);
    printf("  train on %s: best = win %.0f MeV, cos>%.2f, planes>=%d, precuts %d  (train FoM %.3f)\n",
           train?"odd ":"even", 1000*best.win, best.cmin, best.planes, best.precut, bestF);
    printf("     evaluated on %s: S %.2f  B %.1f  (raw %ld / %ld)  FoM %.3f\n", test?"odd ":"even", s, b, ns, nb, fom(s,b,fsys));
    Soos+=s; Boos+=b; nSoos+=ns; nBoos+=nb;
  }
  printf("  OUT-OF-SAMPLE total: S %.2f  B %.1f  (raw %ld / %ld)  purity %.1f%%  FoM(f=0) %.3f  FoM(f=%.0f%%) %.3f\n",
         Soos, Boos, nSoos, nBoos, 100*Soos/(Soos+Boos), fom(Soos,Boos,0), 100*fsys, fom(Soos,Boos,fsys));
  // in-sample reference
  P ref{0.020,0.95,1,0}; double s,b; long ns,nb; double s2,b2; long n2,m2;
  count(ref,0,s,b,ns,nb); count(ref,1,s2,b2,n2,m2);
  printf("  IN-SAMPLE chosen point (20 MeV, cos>0.95, >=1, no precuts): S %.2f  B %.1f  FoM(f=0) %.3f  FoM(f=%.0f%%) %.3f\n",
         s+s2, b+b2, fom(s+s2,b+b2,0), 100*fsys, fom(s+s2,b+b2,fsys));
  // fold stability of the chosen point
  printf("  chosen point by fold: even S %.2f B %.1f (%ld/%ld) | odd S %.2f B %.1f (%ld/%ld)\n", s,b,ns,nb,s2,b2,n2,m2);
  fflush(stdout);
}
