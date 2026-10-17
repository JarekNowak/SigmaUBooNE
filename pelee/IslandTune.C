// ============================================================================
//  IslandTune.C -- tune the island-test parameters on the cache written by
//  IslandCache.C. The seeds land on activity 98% of the time, so the signal
//  losses at 2-5 cm flight are bridging: vertex activity connecting the decay V
//  to the muon under the default connectivity (3 wires, 40 ticks). Nothing here
//  was ever tuned, so a grid is scanned over connectivity, minimum cluster size,
//  seed radius and an ADC occupancy threshold.
//
//  Objective: S/sqrt(S+B+(0.15 B)^2) at the operating point of SigmaLambdaSel.C
//  (no precuts, |m-m_Lambda|<20 MeV, cos alpha>0.95, L>0.5), island on >= 1 plane,
//  with the >= 2 plane variant printed alongside. The default parameters must
//  reproduce IslandAna.C (S=1.58, B=12.0) -- that line is printed first as a check.
//
//  Run:  root -l -b -q 'IslandTune.C("islcache_sig.root","islcache_bkg.root")'
// ============================================================================
#include "Island.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

struct Hit { int w; float t, a; };
struct Ev { int isSig, evt; double w, mass, cosA, len, flight; int seedW[3][3]; float seedT[3]; std::vector<Hit> h[3]; };

static std::vector<Ev> load(const char* fn){
  std::vector<Ev> v; TFile f(fn); TTree* t=(TTree*)f.Get("cache");
  int isSig,evt,firstWire[3],seedW[3][3]; float seedT[3]; double w,mass,len,cosA,fl;
  std::vector<int>* hw[3]={0,0,0}; std::vector<float>* ht[3]={0,0,0}; std::vector<float>* ha[3]={0,0,0};
  t->SetBranchAddress("isSig",&isSig); t->SetBranchAddress("evt",&evt); t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("mass",&mass); t->SetBranchAddress("decayLen",&len); t->SetBranchAddress("cosAlpha",&cosA);
  t->SetBranchAddress("trueFlight",&fl); t->SetBranchAddress("firstWire",firstWire); t->SetBranchAddress("seedW",seedW); t->SetBranchAddress("seedT",seedT);
  for(int p=0;p<3;++p){ t->SetBranchAddress(Form("hw%d",p),&hw[p]); t->SetBranchAddress(Form("ht%d",p),&ht[p]); t->SetBranchAddress(Form("ha%d",p),&ha[p]); }
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); Ev e; e.isSig=isSig; e.evt=evt; e.w=w; e.mass=mass; e.len=len; e.cosA=cosA; e.flight=fl;
    for(int p=0;p<3;++p){ for(int k=0;k<3;++k) e.seedW[p][k]=seedW[p][k]; 
      if(firstWire[p]<0) continue;
      for(size_t k=0;k<hw[p]->size();++k) e.h[p].push_back({(*hw[p])[k],(*ht[p])[k],(*ha[p])[k]});
      std::sort(e.h[p].begin(),e.h[p].end(),[](const Hit&a,const Hit&b){ return a.w<b.w || (a.w==b.w && a.t<b.t); }); }
    for(int k=0;k<3;++k) e.seedT[k]=seedT[k];
    v.push_back(e); }
  return v; }

// same logic as Island::planeTestP, but O(n k) with wire-sorted hits
static bool planeTest(const std::vector<Hit>& H, const int seedW[3], const float seedT[3], const Island::Params& P, bool& ok){
  ok=false;
  std::vector<int> idx; idx.reserve(H.size());
  for(int i=0;i<(int)H.size();++i) if(P.adcMin<0 || H[i].a>=P.adcMin) idx.push_back(i);
  const int n=idx.size(); if(n<P.minclust) return false;
  std::vector<int> par(n); std::iota(par.begin(),par.end(),0);
  std::function<int(int)> find=[&](int x){ while(par[x]!=x){ par[x]=par[par[x]]; x=par[x]; } return x; };
  for(int i=0;i<n;++i){ const Hit& a=H[idx[i]];
    for(int j=i+1;j<n;++j){ const Hit& b=H[idx[j]]; if(b.w-a.w>P.cdw) break;
      if(std::fabs(a.t-b.t)<=P.cdt){ int ra=find(i), rb=find(j); if(ra!=rb) par[ra]=rb; } } }
  std::vector<int> lab(n); for(int i=0;i<n;++i) lab[i]=find(i);
  std::vector<int> csz(n,0); for(int i=0;i<n;++i) csz[lab[i]]++;
  int sl[3]={-1,-1,-1};
  for(int s=0;s<3;++s){ double best=1e18; int bi=-1;
    for(int i=0;i<n;++i){ const Hit& h=H[idx[i]]; const double dw=h.w-seedW[s], dt=(h.t-seedT[s])*0.055/0.3, d=dw*dw+dt*dt; if(d<best){ best=d; bi=i; } }
    if(bi>=0 && best<P.seedR2) sl[s]=lab[bi]; }
  if(sl[0]<0||sl[1]<0||sl[2]<0) return false;
  ok=true;
  if(csz[sl[0]]<P.minclust||csz[sl[1]]<P.minclust||csz[sl[2]]<P.minclust) return false;
  return sl[1]==sl[2] && sl[0]!=sl[1] && sl[0]!=sl[2];
}

struct Res { double S1=0,B1=0,S2=0,B2=0; long nS1=0,nB1=0,nS2=0,nB2=0; };
static Res eval(const std::vector<Ev>& S, const std::vector<Ev>& B, const Island::Params& P, double win, double cmin, std::vector<int>* sigPlanes=nullptr){
  Res r; auto run=[&](const std::vector<Ev>& v, bool sig){
    for(size_t ie=0; ie<v.size(); ++ie){ const Ev& e=v[ie];
      const bool op = std::fabs(e.mass-1.115683)<win && e.cosA>cmin && e.len>0.5;
      int np=0; for(int p=0;p<3;++p){ if(e.h[p].empty()) continue; bool ok; if(planeTest(e.h[p],e.seedW[p],e.seedT,P,ok)) np++; }
      if(sig && sigPlanes) (*sigPlanes)[ie]=np;
      if(!op) continue;
      if(np>=1){ if(sig){ r.S1+=e.w; r.nS1++; } else { r.B1+=e.w; r.nB1++; } }
      if(np>=2){ if(sig){ r.S2+=e.w; r.nS2++; } else { r.B2+=e.w; r.nB2++; } } } };
  run(S,true); run(B,false); return r; }
static double fom(double s,double b,double f=0.15){ return (s+b>0)? s/std::sqrt(s+b+f*f*b*b):0; }

void IslandTune(const char* sigFile="islcache_sig.root", const char* bkgFile="islcache_bkg.root",
                double win=0.020, double cmin=0.95)
{
  auto S=load(sigFile), B=load(bkgFile);
  printf("\n  cache: %zu signal / %zu background events (pre-island superset)\n", S.size(), B.size());
  { double s=0,b=0; long ns=0,nb=0; for(auto&e:S) if(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5){ s+=e.w; ns++; }
    for(auto&e:B) if(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5){ b+=e.w; nb++; }
    printf("  at the operating point before the island: S %.2f (%ld)  B %.1f (%ld)\n", s, ns, b, nb); }

  auto line=[&](const char* tag, const Island::Params& P, const Res& r){
    printf("  %-34s  >=1: S %5.2f B %6.1f (%3ld/%3ld) pur %5.1f%% FoM %.3f   >=2: S %5.2f B %6.1f (%3ld/%3ld) pur %5.1f%% FoM %.3f\n",
      tag, r.S1, r.B1, r.nS1, r.nB1, 100*r.S1/(r.S1+r.B1+1e-9), fom(r.S1,r.B1), r.S2, r.B2, r.nS2, r.nB2, 100*r.S2/(r.S2+r.B2+1e-9), fom(r.S2,r.B2)); };

  Island::Params def; Res r0=eval(S,B,def,win,cmin);
  printf("\n  --- default parameters (must reproduce IslandAna: S 1.58, B 12.0 at >=1) ---\n"); line("cdw 3 cdt 40 min 3 R 20 adc none", def, r0);

  // ADC scale, for the threshold grid
  { std::vector<float> a; for(auto&e:S) for(int p=0;p<3;++p) for(auto&h:e.h[p]) a.push_back(h.a); std::sort(a.begin(),a.end());
    if(!a.empty()) printf("  ADC quantiles (signal hits): 10%% %.1f  25%% %.1f  50%% %.1f  75%% %.1f  90%% %.1f\n",
      a[a.size()/10], a[a.size()/4], a[a.size()/2], a[3*a.size()/4], a[9*a.size()/10]); }

  printf("\n  --- stage 1: connectivity (min 3, R 20, no adc) ---\n");
  Island::Params best=def; double bestF=fom(r0.S1,r0.B1);
  for(int cdw : {1,2,3,4,6}) for(double cdt : {10.,20.,40.,80.,160.}){
    Island::Params P=def; P.cdw=cdw; P.cdt=cdt; Res r=eval(S,B,P,win,cmin);
    line(Form("cdw %d cdt %.0f", cdw, cdt), P, r); if(fom(r.S1,r.B1)>bestF){ bestF=fom(r.S1,r.B1); best=P; } }
  printf("\n  --- stage 2: minimum cluster size and seed radius, at cdw %d cdt %.0f ---\n", best.cdw, best.cdt);
  Island::Params b1=best;
  for(int mc : {1,3,5,10,20}) for(double R : {5.,10.,20.,40.}){
    Island::Params P=b1; P.minclust=mc; P.seedR2=R*R; Res r=eval(S,B,P,win,cmin);
    line(Form("min %d R %.0f", mc, R), P, r); if(fom(r.S1,r.B1)>bestF){ bestF=fom(r.S1,r.B1); best=P; } }
  printf("\n  --- stage 3: ADC occupancy threshold, at cdw %d cdt %.0f min %d R %.0f ---\n", best.cdw, best.cdt, best.minclust, std::sqrt(best.seedR2));
  Island::Params b2=best;
  { std::vector<float> a; for(auto&e:S) for(int p=0;p<3;++p) for(auto&h:e.h[p]) a.push_back(h.a); std::sort(a.begin(),a.end());
    for(double q : {0.0,0.10,0.25,0.50,0.75}){ Island::Params P=b2; P.adcMin = q>0 ? a[(size_t)(q*a.size())] : -1.f; Res r=eval(S,B,P,win,cmin);
      line(Form("adc >= %.1f (q%.0f)", P.adcMin, 100*q), P, r); if(fom(r.S1,r.B1)>bestF){ bestF=fom(r.S1,r.B1); best=P; } } }

  printf("\n  === BEST: cdw %d cdt %.0f min %d R %.0f adc %.1f  FoM(>=1) %.3f ===\n", best.cdw, best.cdt, best.minclust, std::sqrt(best.seedR2), best.adcMin, bestF);
  // flight dependence, default vs best (signal at the operating point)
  std::vector<int> npD(S.size(),0), npB(S.size(),0); eval(S,B,def,win,cmin,&npD); eval(S,B,best,win,cmin,&npB);
  const double edges[]={0,1,2,3,5,10,1e9};
  printf("  flight[cm]     n   default>=1  best>=1   default>=2  best>=2\n");
  for(int b=0;b<6;++b){ long n=0,d1=0,b1c=0,d2=0,b2c=0; for(size_t i=0;i<S.size();++i){ const Ev&e=S[i];
      if(!(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5)) continue; if(!(e.flight>=edges[b]&&e.flight<edges[b+1])) continue;
      n++; if(npD[i]>=1)d1++; if(npB[i]>=1)b1c++; if(npD[i]>=2)d2++; if(npB[i]>=2)b2c++; }
    printf("  %4.0f-%-5.0f %5ld  %8.1f%% %8.1f%%   %8.1f%% %8.1f%%\n", edges[b], edges[b+1]>1e8?99:edges[b+1], n, n?100.0*d1/n:0, n?100.0*b1c/n:0, n?100.0*d2/n:0, n?100.0*b2c/n:0); }
  // odd/even check of the best point
  { Res re, ro; auto cnt=[&](const Island::Params& P, int par, Res& r){ for(auto&e:S){ if(e.evt%2!=par) continue; if(!(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5)) continue; int np=0; for(int p=0;p<3;++p){ if(e.h[p].empty()) continue; bool ok; if(planeTest(e.h[p],e.seedW[p],e.seedT,P,ok)) np++; } if(np>=1){ r.S1+=e.w; r.nS1++; } }
      for(auto&e:B){ if(e.evt%2!=par) continue; if(!(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5)) continue; int np=0; for(int p=0;p<3;++p){ if(e.h[p].empty()) continue; bool ok; if(planeTest(e.h[p],e.seedW[p],e.seedT,P,ok)) np++; } if(np>=1){ r.B1+=e.w; r.nB1++; } } };
    cnt(best,0,re); cnt(best,1,ro);
    printf("  best point by parity: even S %.2f B %.1f (%ld/%ld) FoM %.3f | odd S %.2f B %.1f (%ld/%ld) FoM %.3f\n", re.S1,re.B1,re.nS1,re.nB1,fom(re.S1,re.B1), ro.S1,ro.B1,ro.nS1,ro.nB1,fom(ro.S1,ro.B1)); }
  fflush(stdout);
}

// Evaluate one parameter set fully: both plane requirements, parity split, flight curve.
void IslandPoint(const char* sigFile, const char* bkgFile, int cdw, double cdt, int minclust=3, double R=20.0,
                 double win=0.020, double cmin=0.95)
{
  static std::vector<Ev> S, B; static bool loaded=false;
  if(!loaded){ S=load(sigFile); B=load(bkgFile); loaded=true; }
  Island::Params P; P.cdw=cdw; P.cdt=cdt; P.minclust=minclust; P.seedR2=R*R;
  std::vector<int> np(S.size(),0); Res r=eval(S,B,P,win,cmin,&np);
  printf("\n  === cdw %d cdt %.0f min %d R %.0f ===\n", cdw, cdt, minclust, R);
  printf("  >=1: S %.2f B %.1f (%ld/%ld) pur %.1f%% FoM %.3f | >=2: S %.2f B %.1f (%ld/%ld) pur %.1f%% FoM %.3f\n",
    r.S1,r.B1,r.nS1,r.nB1,100*r.S1/(r.S1+r.B1),fom(r.S1,r.B1), r.S2,r.B2,r.nS2,r.nB2,100*r.S2/(r.S2+r.B2),fom(r.S2,r.B2));
  for(int planes=1; planes<=2; ++planes){
    double se=0,be=0,so=0,bo=0; long nbe=0,nbo=0;
    for(size_t i=0;i<S.size();++i){ const Ev&e=S[i]; if(!(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5)) continue; if(np[i]<planes) continue; (e.evt%2?so:se)+=e.w; }
    for(auto&e:B){ if(!(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5)) continue; int n=0; for(int p=0;p<3;++p){ if(e.h[p].empty()) continue; bool ok; if(planeTest(e.h[p],e.seedW[p],e.seedT,P,ok)) n++; } if(n<planes) continue; if(e.evt%2){ bo+=e.w; nbo++; } else { be+=e.w; nbe++; } }
    printf("  >=%d by parity: even S %.2f B %.1f (%ld) FoM %.3f | odd S %.2f B %.1f (%ld) FoM %.3f\n", planes, se,be,nbe,fom(se,be), so,bo,nbo,fom(so,bo));
  }
  const double edges[]={0,1,2,3,5,10,1e9}; printf("  flight[cm]   n   >=1     >=2\n");
  for(int b=0;b<6;++b){ long n=0,n1=0,n2=0; for(size_t i=0;i<S.size();++i){ const Ev&e=S[i]; if(!(std::fabs(e.mass-1.115683)<win&&e.cosA>cmin&&e.len>0.5)) continue; if(!(e.flight>=edges[b]&&e.flight<edges[b+1])) continue; n++; if(np[i]>=1)n1++; if(np[i]>=2)n2++; }
    printf("  %4.0f-%-4.0f %4ld %6.1f%% %6.1f%%\n", edges[b], edges[b+1]>1e8?99:edges[b+1], n, n?100.0*n1/n:0, n?100.0*n2/n:0); }
  fflush(stdout);
}

// Connectivity scan at the N-plane point with an honest split. The adopted 2x20 was
// chosen on the Run 3 list (bkgTrain); the RHC extension (bkgTest) played no part in
// that choice, so its column is out of sample. Each background cache must be weighted
// to data POT on its own; the full-sample B is their POT-weighted mean. An optional
// second signal cache (NuWro) is shown as a pass rate only -- its normalisation is a
// generator choice, see the README.
//  Run:  root -l -b -q -e '.L IslandTune.C' -e 'IslandScan2("islcache_sig.root","islcache_bkg.root","islcache_bkg_ext.root","islcache_sig_nuwro.root")'
void IslandScan2(const char* sigFile, const char* bkgTrain, const char* bkgTest, const char* sig2File="",
                 double potTrain=2.7795e21, double potTest=4.2693e21, int planes=2,
                 double win=0.020, double cmin=0.95)
{
  auto S=load(sigFile), BT=load(bkgTrain), BX=load(bkgTest);
  std::vector<Ev> S2; if(sig2File && *sig2File) S2=load(sig2File);
  const double fT=potTrain/(potTrain+potTest), fX=potTest/(potTrain+potTest);
  auto atOp=[&](const Ev& e){ return std::fabs(e.mass-1.115683)<win && e.cosA>cmin && e.len>0.5; };
  long nS2op=0; for(auto&e:S2) if(atOp(e)) nS2op++;
  printf("\n  >= %d planes; cache sizes: sig %zu, train bkg %zu, test bkg %zu, sig2 %zu (%ld at the op point)\n",
         planes, S.size(), BT.size(), BX.size(), S2.size(), nS2op);
  printf("  %-12s %12s | %16s %16s %16s | %21s | %s\n", "connectivity", "S (raw)", "B train (raw)", "B test (raw)",
         "B full (raw)", "FoM15 train/test/full", "sig2 pass");
  struct Row { int cdw; double cdt, s, bt, bx, ft, fx; };
  std::vector<Row> rows;
  for(int cdw : {1,2,3,4,6}) for(double cdt : {10.,20.,30.,40.,60.,80.}){
    Island::Params P; P.cdw=cdw; P.cdt=cdt;
    Res rT=eval(S,BT,P,win,cmin), rX=eval(S,BX,P,win,cmin);
    const double s = planes>=2 ? rT.S2 : rT.S1;  const long ns = planes>=2 ? rT.nS2 : rT.nS1;
    const double bt= planes>=2 ? rT.B2 : rT.B1;  const long nbt= planes>=2 ? rT.nB2 : rT.nB1;
    const double bx= planes>=2 ? rX.B2 : rX.B1;  const long nbx= planes>=2 ? rX.nB2 : rX.nB1;
    const double bf= fT*bt + fX*bx;
    const double ebt= nbt? bt/std::sqrt((double)nbt):0, ebx= nbx? bx/std::sqrt((double)nbx):0;
    const double ebf= std::sqrt(fT*fT*ebt*ebt + fX*fX*ebx*ebx);
    long n2=0; for(auto&e:S2){ if(!atOp(e)) continue; int np=0;
      for(int p=0;p<3;++p){ if(e.h[p].empty()) continue; bool ok; if(planeTest(e.h[p],e.seedW[p],e.seedT,P,ok)) np++; }
      if(np>=planes) n2++; }
    rows.push_back({cdw,cdt,s,bt,bx,fom(s,bt),fom(s,bx)});
    printf("  cdw %d cdt %2.0f %6.2f (%3ld) | %6.2f+-%4.2f (%3ld) %6.2f+-%4.2f (%3ld) %6.2f+-%4.2f (%3ld) | %.3f / %.3f / %.3f | %5.1f%%\n",
           cdw, cdt, s, ns, bt, ebt, nbt, bx, ebx, nbx, bf, ebf, nbt+nbx, fom(s,bt), fom(s,bx), fom(s,bf),
           nS2op? 100.0*n2/nS2op : 0.0);
  }
  auto bestT=std::max_element(rows.begin(),rows.end(),[](const Row&a,const Row&b){ return a.ft<b.ft; });
  auto bestX=std::max_element(rows.begin(),rows.end(),[](const Row&a,const Row&b){ return a.fx<b.fx; });
  printf("\n  chosen on train: cdw %d cdt %.0f  (train FoM %.3f -> test FoM %.3f)\n", bestT->cdw, bestT->cdt, bestT->ft, bestT->fx);
  printf("  best on test   : cdw %d cdt %.0f  (test FoM %.3f)\n", bestX->cdw, bestX->cdt, bestX->fx);
  fflush(stdout);
}
