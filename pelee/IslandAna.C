// ============================================================================
//  IslandAna.C -- analysis of the cached island-finding pass (IslandProcess.C).
//  Reads isl_<sig>.root and isl_bkg.root; every number is POT-scaled to data
//  (1.1526e21) through the `weight` branch. Blind: signal + overlay only.
//
//  Answers, in order:
//    1. the Lambda-style cut flow with the island test at the end, S/B/purity,
//       and S/sqrt(S+B+(fB)^2) for background systematics f = 0, 7, 15, 30%;
//    2. where the island test belongs -- straight after the pair, or last --
//       and whether one plane or two should be demanded;
//    3. what the test does as a function of true Lambda flight (signal only);
//    4. what the surviving background actually is (mode, true candidate pdgs);
//    5. a scan of the mass window and pointing cut with the island fixed, since
//       the thresholds in SigmaLambdaSel.C were set by intuition, not optimised.
//
//  Run:  root -l -b -q 'IslandAna.C("isl_sig_genie.root","isl_bkg.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <string>

static const double M_LAM = 1.115683;

struct Ev {
  int isSig, mode, stage, nTr, nSh, hasPair, precutOK, hasPhoton;
  double w, mass, decayLen, cosAlpha, sep, llrP, llrPi, shrE, trueFlight;
  int islPass[3], islOK[3], islNPlanes, islAny, islSeeded, muTrue, prTrue, piTrue, s1True;
};
static std::vector<Ev> load(const char* fn){
  std::vector<Ev> v; TFile f(fn); TTree* t=(TTree*)f.Get("isl"); if(!t){ printf("  no tree in %s\n",fn); return v; }
  Ev e;
  t->SetBranchAddress("isSig",&e.isSig); t->SetBranchAddress("mode",&e.mode); t->SetBranchAddress("stage",&e.stage);
  t->SetBranchAddress("nTr",&e.nTr); t->SetBranchAddress("nSh",&e.nSh); t->SetBranchAddress("hasPair",&e.hasPair);
  t->SetBranchAddress("precutOK",&e.precutOK); t->SetBranchAddress("hasPhoton",&e.hasPhoton);
  t->SetBranchAddress("weight",&e.w); t->SetBranchAddress("mass",&e.mass); t->SetBranchAddress("decayLen",&e.decayLen);
  t->SetBranchAddress("cosAlpha",&e.cosAlpha); t->SetBranchAddress("sep",&e.sep); t->SetBranchAddress("llrP",&e.llrP);
  t->SetBranchAddress("llrPi",&e.llrPi); t->SetBranchAddress("shrE",&e.shrE); t->SetBranchAddress("trueFlight",&e.trueFlight);
  t->SetBranchAddress("islPass",e.islPass); t->SetBranchAddress("islOK",e.islOK); t->SetBranchAddress("islNPlanes",&e.islNPlanes);
  t->SetBranchAddress("islAny",&e.islAny); t->SetBranchAddress("islSeeded",&e.islSeeded);
  t->SetBranchAddress("muTrue",&e.muTrue); t->SetBranchAddress("prTrue",&e.prTrue); t->SetBranchAddress("piTrue",&e.piTrue);
  t->SetBranchAddress("s1True",&e.s1True);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); v.push_back(e); }
  return v;
}

// selection pieces, so cut orders can be composed freely
static bool cMass (const Ev& e, double win){ return e.mass>0 && std::fabs(e.mass-M_LAM)<win; }
static bool cAlpha(const Ev& e, double cmin, double lmin){ return e.decayLen>lmin && e.cosAlpha>cmin; }
static bool cIsl  (const Ev& e, int nplanes){ return e.islNPlanes>=nplanes; }

struct SB { double S=0,B=0; long nS=0,nB=0; };
static void row(const char* name, const SB& r){
  const double pur = (r.S+r.B>0)?100*r.S/(r.S+r.B):0;
  auto fom=[&](double f){ return (r.S+r.B>0)?r.S/std::sqrt(r.S+r.B+(f*r.B)*(f*r.B)):0; };
  printf("  %-34s %8.2f %10.1f %6ld %7ld %7.2f %7.3f %7.3f %7.3f %7.3f\n", name, r.S, r.B, r.nS, r.nB, pur,
         fom(0), fom(0.07), fom(0.15), fom(0.30));
}
static void head(){ printf("\n  %-34s %8s %10s %6s %7s %7s %7s %7s %7s %7s\n", "selection","S","B","rawS","rawB","pur[%]",
                           "FoM","f=7%","f=15%","f=30%"); }
template<class F> static SB count(const std::vector<Ev>& S, const std::vector<Ev>& B, F sel){
  SB r; for(auto& e:S) if(sel(e)){ r.S+=e.w; r.nS++; } for(auto& e:B) if(sel(e)){ r.B+=e.w; r.nB++; } return r; }

void IslandAna(const char* sigFile="isl_sig_genie.root", const char* bkgFile="isl_bkg.root")
{
  auto S = load(sigFile), B = load(bkgFile);
  printf("\n  loaded %zu signal / %zu background events at MuonCandidate  (%s, %s)\n", S.size(), B.size(), sigFile, bkgFile);

  // ── 1. the Lambda-style chain, cumulative ───────────────────────────────
  head();
  row("MuonCandidate",                  count(S,B,[](const Ev&e){ return true; }));
  row("+ any pair",                     count(S,B,[](const Ev&e){ return e.hasPair; }));
  row("+ precuts (PairCandidate)",      count(S,B,[](const Ev&e){ return e.precutOK; }));
  row("+ Lambda mass +-40 MeV",         count(S,B,[](const Ev&e){ return e.precutOK&&cMass(e,0.040); }));
  row("+ cos alpha>0.90, L>0.5",        count(S,B,[](const Ev&e){ return e.precutOK&&cMass(e,0.040)&&cAlpha(e,0.90,0.5); }));
  row("+ photon",                       count(S,B,[](const Ev&e){ return e.precutOK&&cMass(e,0.040)&&cAlpha(e,0.90,0.5)&&e.hasPhoton; }));
  row("+ island >=1 plane  [FINAL]",    count(S,B,[](const Ev&e){ return e.precutOK&&cMass(e,0.040)&&cAlpha(e,0.90,0.5)&&e.hasPhoton&&cIsl(e,1); }));
  row("+ island >=2 planes",            count(S,B,[](const Ev&e){ return e.precutOK&&cMass(e,0.040)&&cAlpha(e,0.90,0.5)&&e.hasPhoton&&cIsl(e,2); }));

  // ── 2. placement of the island test ─────────────────────────────────────
  printf("\n  --- island placement: what it removes at each position ---");
  head();
  row("any pair + island>=1",           count(S,B,[](const Ev&e){ return e.hasPair&&cIsl(e,1); }));
  row("any pair + island>=2",           count(S,B,[](const Ev&e){ return e.hasPair&&cIsl(e,2); }));
  row("precuts + island>=1",            count(S,B,[](const Ev&e){ return e.precutOK&&cIsl(e,1); }));
  row("precuts + island>=2",            count(S,B,[](const Ev&e){ return e.precutOK&&cIsl(e,2); }));
  row("precuts + island>=1 + mass",     count(S,B,[](const Ev&e){ return e.precutOK&&cIsl(e,1)&&cMass(e,0.040); }));
  row("precuts + island>=1 + mass + ph",count(S,B,[](const Ev&e){ return e.precutOK&&cIsl(e,1)&&cMass(e,0.040)&&e.hasPhoton; }));
  row("any pair + island>=1 + mass + ph",count(S,B,[](const Ev&e){ return e.hasPair&&cIsl(e,1)&&cMass(e,0.040)&&e.hasPhoton; }));
  row("any pair + island>=2 + mass + ph",count(S,B,[](const Ev&e){ return e.hasPair&&cIsl(e,2)&&cMass(e,0.040)&&e.hasPhoton; }));

  // per-plane behaviour, on precut pairs
  printf("\n  --- per plane, precut-passing pairs (pass = p+pi merged & mu separate; seeded = all 3 seeds on activity) ---\n");
  printf("  %-6s %10s %10s %10s %10s\n","plane","sig pass","sig seeded","bkg pass","bkg seeded");
  for(int p=0;p<3;++p){
    double sp=0,ss=0,st=0,bp=0,bs=0,bt=0;
    for(auto&e:S) if(e.precutOK){ st+=e.w; if(e.islPass[p]) sp+=e.w; if(e.islOK[p]) ss+=e.w; }
    for(auto&e:B) if(e.precutOK){ bt+=e.w; if(e.islPass[p]) bp+=e.w; if(e.islOK[p]) bs+=e.w; }
    printf("  %-6d %9.1f%% %9.1f%% %9.1f%% %9.1f%%\n", p, st?100*sp/st:0, st?100*ss/st:0, bt?100*bp/bt:0, bt?100*bs/bt:0);
  }

  // ── 3. flight dependence (signal) ───────────────────────────────────────
  printf("\n  --- signal island pass rate vs true Lambda flight, precut-passing pairs ---\n");
  const double edges[]={0,1,2,3,5,10,1e9}; const int nb=6;
  printf("  %-10s %6s %8s %8s %8s\n","flight[cm]","n","pass>=1","pass>=2","seeded3");
  for(int b=0;b<nb;++b){ long n=0,p1=0,p2=0,s3=0;
    for(auto&e:S) if(e.precutOK && e.trueFlight>=edges[b] && e.trueFlight<edges[b+1]){ n++; if(e.islNPlanes>=1)p1++; if(e.islNPlanes>=2)p2++; if(e.islSeeded==3)s3++; }
    printf("  %4.0f-%-5.0f %6ld %7.1f%% %7.1f%% %7.1f%%\n", edges[b], edges[b+1]>1e8?99:edges[b+1], n, n?100.0*p1/n:0, n?100.0*p2/n:0, n?100.0*s3/n:0); }
  { long nAll=0,nLong=0,nPre=0,nPreLong=0;
    for(auto&e:S){ nAll++; if(e.trueFlight>2)nLong++; if(e.precutOK){ nPre++; if(e.trueFlight>2)nPreLong++; } }
    printf("  flight > 2 cm: %.1f%% at MuonCandidate, %.1f%% after precuts\n", nAll?100.0*nLong/nAll:0, nPre?100.0*nPreLong/nPre:0); }

  // ── 4. what the surviving background is ─────────────────────────────────
  auto finalSel=[](const Ev&e){ return e.precutOK&&cMass(e,0.040)&&cAlpha(e,0.90,0.5)&&e.hasPhoton&&cIsl(e,1); };
  auto looseSel=[](const Ev&e){ return e.hasPair&&cIsl(e,1)&&cMass(e,0.040)&&e.hasPhoton; };
  const char* modeN[]={"QE","RES","DIS","COH","MEC","Other"};
  auto mi=[](int m){ return m==0?0:m==1?1:m==2?2:m==3?3:m==10?4:5; };
  for(int which=0; which<2; ++which){
    printf("\n  --- surviving BACKGROUND composition, %s ---\n", which==0?"FINAL selection":"loose (any pair + island + mass + photon)");
    double tot=0, byMode[6]={0}, truePair=0, trueP=0, truePi=0, trueMu=0, trueG=0, trueV=0; long raw=0;
    for(auto&e:B){ if(!(which==0?finalSel(e):looseSel(e))) continue; raw++; tot+=e.w; byMode[mi(e.mode)]+=e.w;
      if(e.prTrue==2212) trueP+=e.w; if(e.piTrue==-211) truePi+=e.w; if(e.prTrue==2212&&e.piTrue==-211) truePair+=e.w;
      if(std::abs(e.muTrue)==13) trueMu+=e.w; if(e.s1True==22) trueG+=e.w; }
    printf("  raw %ld, weighted %.1f\n", raw, tot);
    if(tot>0){ printf("  mode:"); for(int m=0;m<6;++m) printf("  %s %.0f%%", modeN[m], 100*byMode[m]/tot); printf("\n");
      printf("  candidate truth: muon is mu %.0f%%, proton is p %.0f%%, pion is pi- %.0f%%, both p AND pi- %.0f%%, shower is gamma %.0f%%\n",
             100*trueMu/tot, 100*trueP/tot, 100*truePi/tot, 100*truePair/tot, 100*trueG/tot); }
    // the same for signal, for reference
    double stot=0, sPair=0; for(auto&e:S){ if(!(which==0?finalSel(e):looseSel(e))) continue; stot+=e.w; if(e.prTrue==2212&&e.piTrue==-211) sPair+=e.w; }
    if(stot>0) printf("  (signal: both candidates correct %.0f%%)\n", 100*sPair/stot);
  }

  // ── 5. threshold scan with the island fixed ─────────────────────────────
  printf("\n  --- scan: precuts + photon + island>=1, then mass window x cos(alpha) (L>0.5) ; cells show S / B / FoM(f=15%%) ---\n");
  const double wins[]={0.020,0.040,0.060,0.100,9.0}; const double cas[]={-2.0,0.0,0.5,0.8,0.90,0.95,0.98};
  printf("  %-8s","win\\cos"); for(double c:cas) printf(" %19s", c<-1?"none":Form("cos>%.2f",c)); printf("\n");
  for(double w:wins){ printf("  %-8s", w>1?"none":Form("%.0fMeV",w*1000));
    for(double c:cas){ SB r=count(S,B,[&](const Ev&e){ return e.precutOK&&e.hasPhoton&&cIsl(e,1)&&(w>1||cMass(e,w))&&(c<-1||cAlpha(e,c,0.5)); });
      const double f=(r.S+r.B>0)?r.S/std::sqrt(r.S+r.B+0.0225*r.B*r.B):0; printf(" %5.2f/%6.1f/%5.3f", r.S, r.B, f); }
    printf("\n"); }
  printf("\n  --- same scan WITHOUT the precuts (any pair) ---\n");
  printf("  %-8s","win\\cos"); for(double c:cas) printf(" %19s", c<-1?"none":Form("cos>%.2f",c)); printf("\n");
  for(double w:wins){ printf("  %-8s", w>1?"none":Form("%.0fMeV",w*1000));
    for(double c:cas){ SB r=count(S,B,[&](const Ev&e){ return e.hasPair&&e.hasPhoton&&cIsl(e,1)&&(w>1||cMass(e,w))&&(c<-1||cAlpha(e,c,0.5)); });
      const double f=(r.S+r.B>0)?r.S/std::sqrt(r.S+r.B+0.0225*r.B*r.B):0; printf(" %5.2f/%6.1f/%5.3f", r.S, r.B, f); }
    printf("\n"); }
  fflush(stdout);
}
