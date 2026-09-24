// ============================================================================
//  PhotonMichel.C -- two unexploited handles, measured on the v2 island cache
//  (IslandProcess.C): (3) the Sigma0 photon, which the Lambda-style chain has
//  so far not used at all (the photon requirement is a no-op after NShowers), and
//  (4) a Michel-electron tag on the muon candidate, since 68% of the surviving
//  background muons are true mu- (wrong-sign numu CC) while the signal muon is
//  a mu+: mu+ always decays, mu- captures 75% of the time in argon.
//
//  Two samples: the pre-island stage (any pair, |m-mL|<20 MeV, cos alpha>0.95,
//  L>0.5: 144 / 267 raw) for shapes, and the operating point (+ island >= 1
//  plane, old 3x40 connectivity as cached: 82 / 29 raw) for the numbers.
//
//  Run:  root -l -b -q 'PhotonMichel.C("isl_sig_genie.root","isl_bkg.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>

static const double ML=1.115683, MS=1.192642;
struct E { int isSig, hasPair, islN, s1True, lepTrue, muContained, michelN5, mode; double w, mass, len, cosA, sigmaMass, cosPhLam, shrE, phDist, phDedx, phMoliere, truePhE, michelDist, michelE; };

static double fom(double s,double b,double f=0.15){ return (s+b>0)? s/std::sqrt(s+b+f*f*b*b):0; }

void PhotonMichel(const char* sigFile="isl_sig_genie.root", const char* bkgFile="isl_bkg.root")
{
  // load with explicit names (islNPlanes->islN, weight->w, decayLen->len, cosAlpha->cosA)
  auto loadX=[](const char* fn){ std::vector<E> v; TFile f(fn); TTree* t=(TTree*)f.Get("isl"); E e;
    t->SetBranchAddress("isSig",&e.isSig); t->SetBranchAddress("hasPair",&e.hasPair); t->SetBranchAddress("islNPlanes",&e.islN);
    t->SetBranchAddress("s1True",&e.s1True); t->SetBranchAddress("lepTrue",&e.lepTrue); t->SetBranchAddress("muContained",&e.muContained);
    t->SetBranchAddress("michelN5",&e.michelN5); t->SetBranchAddress("mode",&e.mode); t->SetBranchAddress("weight",&e.w);
    t->SetBranchAddress("mass",&e.mass); t->SetBranchAddress("decayLen",&e.len); t->SetBranchAddress("cosAlpha",&e.cosA);
    t->SetBranchAddress("sigmaMass",&e.sigmaMass); t->SetBranchAddress("cosPhLam",&e.cosPhLam); t->SetBranchAddress("shrE",&e.shrE);
    t->SetBranchAddress("phDist",&e.phDist); t->SetBranchAddress("phDedx",&e.phDedx); t->SetBranchAddress("phMoliere",&e.phMoliere);
    t->SetBranchAddress("truePhE",&e.truePhE); t->SetBranchAddress("michelDist",&e.michelDist); t->SetBranchAddress("michelE",&e.michelE);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); v.push_back(e); } return v; };
  auto S=loadX(sigFile), B=loadX(bkgFile);
  auto pre=[](const E& e){ return e.hasPair && std::fabs(e.mass-ML)<0.020 && e.cosA>0.95 && e.len>0.5; };
  auto op =[&](const E& e){ return pre(e) && e.islN>=1; };
  auto tab=[&](const char* title, std::function<bool(const E&)> sel, std::function<bool(const E&)> cut){
    double s=0,b=0,sc=0,bc=0; long ns=0,nb=0,nsc=0,nbc=0;
    for(auto&e:S) if(sel(e)){ s+=e.w; ns++; if(cut(e)){ sc+=e.w; nsc++; } }
    for(auto&e:B) if(sel(e)){ b+=e.w; nb++; if(cut(e)){ bc+=e.w; nbc++; } }
    printf("  %-44s keeps sig %5.1f%% (%3ld/%3ld)  bkg %5.1f%% (%3ld/%3ld)   S %.2f B %5.1f  FoM15 %.3f -> %.3f\n",
           title, s?100*sc/s:0, nsc, ns, b?100*bc/b:0, nbc, nb, sc, bc, fom(s,b), fom(sc,bc)); };
  auto hist=[&](const char* title, std::function<bool(const E&)> sel, std::function<double(const E&)> var, const std::vector<double>& edges, const char* unit){
    printf("\n  %s [%s]\n  %-16s %8s %8s %6s\n", title, unit, "bin", "sig[%]", "bkg[%]", "rawB");
    double st=0,bt=0; for(auto&e:S) if(sel(e)) st+=e.w; for(auto&e:B) if(sel(e)) bt+=e.w;
    for(size_t k=0;k+1<edges.size();++k){ double s=0,b=0; long nb=0;
      for(auto&e:S) if(sel(e)&&var(e)>=edges[k]&&var(e)<edges[k+1]) s+=e.w;
      for(auto&e:B) if(sel(e)&&var(e)>=edges[k]&&var(e)<edges[k+1]){ b+=e.w; nb++; }
      printf("  %7.3g-%-7.3g %7.1f%% %7.1f%% %6ld\n", edges[k], edges[k+1], st?100*s/st:0, bt?100*b/bt:0, nb); } };

  printf("\n================ STEP 3: the photon ================\n");
  { double n=0,g=0; for(auto&e:S) if(pre(e)){ n+=e.w; if(e.s1True==22) g+=e.w; } printf("  signal at pre-island: photon candidate is a true photon %.0f%%\n", n?100*g/n:0);
    double nb=0,gb=0; for(auto&e:B) if(pre(e)){ nb+=e.w; if(e.s1True==22) gb+=e.w; } printf("  background            : photon candidate is a true photon %.0f%%\n", nb?100*gb/nb:0); }
  hist("Sigma mass - Lambda mass, pre-island (true 77 MeV)", pre, [](const E&e){ return e.sigmaMass>0 ? 1000*(e.sigmaMass-e.mass) : -1; }, {-1,0,20,40,60,80,100,130,170,250,400,1000}, "MeV");
  hist("photon energy, pre-island", pre, [](const E&e){ return e.shrE>0? 1000*e.shrE : -1; }, {-1,0,20,40,60,80,120,160,250,400,10000}, "MeV");
  hist("cos(photon, Lambda), pre-island", pre, [](const E&e){ return e.cosPhLam; }, {-2,-1,-0.5,0,0.3,0.6,0.8,0.9,1.001}, "");
  hist("photon start distance to vertex, pre-island", pre, [](const E&e){ return e.phDist; }, {0,2,5,10,15,20,30,50,1000}, "cm");
  hist("true photon momentum, signal only", [&](const E&e){ return pre(e)&&e.isSig; }, [](const E&e){ return 1000*e.truePhE; }, {0,40,60,80,100,120,150,200,400}, "MeV");
  printf("\n  --- photon cuts at the OPERATING POINT (island>=1) ---\n");
  tab("Sigma-Lambda mass diff in [40,130] MeV", op, [](const E&e){ const double d=1000*(e.sigmaMass-e.mass); return e.sigmaMass>0 && d>40 && d<130; });
  tab("Sigma-Lambda mass diff in [20,170] MeV", op, [](const E&e){ const double d=1000*(e.sigmaMass-e.mass); return e.sigmaMass>0 && d>20 && d<170; });
  tab("Sigma mass diff < 170 MeV (veto high only)", op, [](const E&e){ const double d=1000*(e.sigmaMass-e.mass); return !(e.sigmaMass>0) || d<170; });
  tab("photon E < 250 MeV", op, [](const E&e){ return !(e.shrE>0) || e.shrE<0.25; });
  tab("photon E < 160 MeV", op, [](const E&e){ return !(e.shrE>0) || e.shrE<0.16; });
  tab("photon E in [20,250] MeV", op, [](const E&e){ return e.shrE>0.02 && e.shrE<0.25; });
  tab("photon dist to vertex > 2 cm (conversion gap)", op, [](const E&e){ return e.phDist>2; });
  tab("photon dist to vertex > 5 cm", op, [](const E&e){ return e.phDist>5; });
  printf("  --- same, at the pre-island stage (more statistics) ---\n");
  tab("Sigma-Lambda mass diff in [40,130] MeV", pre, [](const E&e){ const double d=1000*(e.sigmaMass-e.mass); return e.sigmaMass>0 && d>40 && d<130; });
  tab("photon E < 250 MeV", pre, [](const E&e){ return !(e.shrE>0) || e.shrE<0.25; });
  tab("photon dist to vertex > 2 cm", pre, [](const E&e){ return e.phDist>2; });

  printf("\n================ STEP 4: the Michel tag ================\n");
  { const char* nm[2]={"signal","background"}; for(int w=0;w<2;++w){ auto& V = w?B:S; double n=0,mp=0,mm=0,cont=0;
      for(auto&e:V) if(pre(e)){ n+=e.w; if(e.lepTrue==-13) mp+=e.w; if(e.lepTrue==13) mm+=e.w; if(e.muContained) cont+=e.w; }
      printf("  %-10s pre-island: true lepton mu+ %.0f%%  mu- %.0f%%  ; reco muon end contained %.0f%%\n", nm[w], n?100*mp/n:0, n?100*mm/n:0, n?100*cont/n:0); } }
  hist("distance muon end -> nearest other PFP start, contained muons, pre-island", [&](const E&e){ return pre(e)&&e.muContained; }, [](const E&e){ return e.michelDist; }, {0,1,2,3,5,8,12,20,50,1e4}, "cm");
  // physics validation: Michel rate for true mu+ vs true mu- in BACKGROUND (both present there)
  { double np=0,mp=0,nm=0,mm=0; for(auto&e:B) if(pre(e)&&e.muContained){ if(e.lepTrue==-13){ np+=e.w; if(e.michelDist>=0&&e.michelDist<5) mp+=e.w; } if(e.lepTrue==13){ nm+=e.w; if(e.michelDist>=0&&e.michelDist<5) mm+=e.w; } }
    printf("  background, contained: PFP within 5 cm of muon end for true mu+ %.0f%% vs true mu- %.0f%%\n", np?100*mp/np:0, nm?100*mm/nm:0);
    double ns=0,ms=0; for(auto&e:S) if(pre(e)&&e.muContained){ ns+=e.w; if(e.michelDist>=0&&e.michelDist<5) ms+=e.w; } printf("  signal, contained    : %.0f%%\n", ns?100*ms/ns:0); }
  printf("\n  --- Michel cuts at the OPERATING POINT (island>=1) ---\n");
  tab("require PFP within 5 cm of muon end (contained only)", op, [](const E&e){ return !e.muContained || (e.michelDist>=0 && e.michelDist<5); });
  tab("require PFP within 3 cm of muon end (contained only)", op, [](const E&e){ return !e.muContained || (e.michelDist>=0 && e.michelDist<3); });
  tab("require Michel-like SHOWER within 5 cm (contained)", op, [](const E&e){ return !e.muContained || (e.michelDist>=0 && e.michelDist<5 && e.michelE>0 && e.michelE<0.08); });
  tab("require contained muon", op, [](const E&e){ return e.muContained; });
  printf("  --- same, pre-island ---\n");
  tab("require PFP within 5 cm of muon end (contained only)", pre, [](const E&e){ return !e.muContained || (e.michelDist>=0 && e.michelDist<5); });
  tab("require contained muon", pre, [](const E&e){ return e.muContained; });
  fflush(stdout);
}
