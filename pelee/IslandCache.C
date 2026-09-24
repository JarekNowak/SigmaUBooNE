// ============================================================================
//  IslandCache.C -- cache the raw wire windows (hits, ADC, per-hit truth pdg) and
//  the three candidate seeds for every event at the PRE-island stage of the
//  Lambda-style chain, so the island parameters can be tuned offline
//  (IslandTune.C) without re-reading 8000 files per trial.
//
//  Pre-island stage here is a SUPERSET of the operating point: any pair (no
//  precuts), |m - m_Lambda| < 60 MeV, decay length > 0.5 cm, cos alpha > 0.80.
//
//  Run:  root -l -b -q 'IslandCache.C("lists/sig_genie.list",0,true,"islcache_sig.root",1.9284e-2)'
//        root -l -b -q 'IslandCache.C("lists/bkg_overlay.list",0,false,"islcache_bkg.root",1.6536e-1)'
// ============================================================================
#include "HypInput.h"
#include "SigmaFV.h"
#include "Island.h"
#include "LambdaStyle.h"
#include "TChain.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <fstream>
#include <vector>
#include <cstdio>
using namespace LambdaStyle;

static const double MUON_TRK_SCORE=0.8, MUON_TRK_DIST=20.0, MUON_TRK_LEN=10.0, MUON_LLR=0.7;
static const int NSHOWERS_MAX=4, NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5, TOPO_MIN=0.20;
static const double PRE_WIN=0.060, PRE_COS=0.80, PRE_LEN=0.5;

void IslandCache(const char* list="lists/sig_genie.list", int maxFiles=0, bool signalOnly=true,
                 const char* outName="islcache_sig.root", double potScale=1.0)
{
  // vector<vector<int>> has no compiled dictionary in this ROOT; without it the
  // per-hit truth pdg branch reads back corrupted.
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf=0;
  while(std::getline(in,l)) if(l.size() && (maxFiles<=0 || nf<maxFiles)){ ch.Add(l.c_str()); nf++; }
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){ printf("  not a hypana tree\n"); return; }
  int fmSlice=-1; const bool haveSlice=ch.GetBranch("flash_match_nu_slice_ID")!=nullptr;
  if(haveSlice) ch.SetBranchAddress("flash_match_nu_slice_ID",&fmSlice);
  std::vector<double>* d3=nullptr; const bool haveDedx=ch.GetBranch("trk_three_plane_dedx")!=nullptr;
  if(haveDedx) ch.SetBranchAddress("trk_three_plane_dedx",&d3);
  std::vector<std::vector<float>>* win[3]={0,0,0}; std::vector<std::vector<float>>* adc[3]={0,0,0};
  std::vector<std::vector<int>>* pdg[3]={0,0,0}; std::vector<int>* vw=nullptr;
  for(int p=0;p<3;++p){ ch.SetBranchAddress(Form("ct_test_window_plane%d",p),&win[p]);
    ch.SetBranchAddress(Form("ct_test_adc_plane%d",p),&adc[p]); ch.SetBranchAddress(Form("ct_test_pdg_plane%d",p),&pdg[p]); }
  ch.SetBranchAddress("ct_test_primary_vtx_wires",&vw);
  double lsx=0,lsy=0,lsz=0,lex=0,ley=0,lez=0;
  auto B=[&](const char* n, void* q){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,q); };
  B("mc_sigmazero_lambda_start_x",&lsx); B("mc_sigmazero_lambda_start_y",&lsy); B("mc_sigmazero_lambda_start_z",&lsz);
  B("mc_sigmazero_lambda_end_x",&lex);   B("mc_sigmazero_lambda_end_y",&ley);   B("mc_sigmazero_lambda_end_z",&lez);

  TFile fout(outName,"RECREATE"); TTree T("cache","pre-island events with wire windows");
  int isSig=0, evt=0, precutOK=0; double weight=potScale, mass=-1, decayLen=-1, cosAlpha=-2, trueFlight=-1;
  int firstWire[3]={0,0,0}; int seedW[3][3]; float seedT[3];   // seedW[plane][mu,p,pi], seedT[mu,p,pi]
  std::vector<int> hw[3]; std::vector<float> ht[3], ha[3]; std::vector<int> hp[3];
  T.Branch("isSig",&isSig,"isSig/I"); T.Branch("evt",&evt,"evt/I"); T.Branch("precutOK",&precutOK,"precutOK/I");
  T.Branch("weight",&weight,"weight/D"); T.Branch("mass",&mass,"mass/D"); T.Branch("decayLen",&decayLen,"decayLen/D");
  T.Branch("cosAlpha",&cosAlpha,"cosAlpha/D"); T.Branch("trueFlight",&trueFlight,"trueFlight/D");
  T.Branch("firstWire",firstWire,"firstWire[3]/I"); T.Branch("seedW",seedW,"seedW[3][3]/I"); T.Branch("seedT",seedT,"seedT[3]/F");
  for(int p=0;p<3;++p){ T.Branch(Form("hw%d",p),&hw[p]); T.Branch(Form("ht%d",p),&ht[p]); T.Branch(Form("ha%d",p),&ha[p]); T.Branch(Form("hp%d",p),&hp[p]); }

  const PairCuts pc; long nWritten=0; const Long64_t N=ch.GetEntries();
  printf("\n  %s : %d files, %lld events\n", signalOnly?"SIGNAL":"BACKGROUND", nf, N); fflush(stdout);
  for(Long64_t ie=0; ie<N; ++ie){
    ch.GetEntry(ie); H.convert();
    isSig=(H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
    if(signalOnly!=(bool)isSig) continue;
    if(haveSlice && fmSlice<0) continue;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int n=(int)H.trk_score.size(); if(n<1) continue;
    const TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int i){ return H.trk_score[i]>TRK_SHR_SPLIT; };
    std::vector<TVector3> st(n),dir(n); std::vector<double> dv(n,1e9);
    for(int i=0;i<n;++i){ if(isTrk(i)) st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]); else st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]);
      dir[i].SetXYZ(H.trk_dx[i],H.trk_dy[i],H.trk_dz[i]); dv[i]=(st[i]-vtx).Mag(); }
    int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue; if(!(ntr<NTRACKS_MAX)) continue;
    int mu=-1; for(int i=0;i<n;++i) if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dv[i]<MUON_TRK_DIST&&H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR) if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
    if(mu<0) continue;
    int pr=-1,pi=-1,prA=-1,piA=-1; double best=-1e30,bestA=-1e30;
    for(int a=0;a<n;++a){ if(a==mu||!isTrk(a)) continue; for(int b=0;b<n;++b){ if(b==mu||b==a||!isTrk(b)) continue;
      PairFeat f; f.sep=(st[a]-st[b]).Mag(); f.scoreP=H.trk_score[a]; f.scorePi=H.trk_score[b]; f.llrP=H.trk_llr[a]; f.llrPi=H.trk_llr[b];
      if(haveDedx&&d3&&a<(int)d3->size()&&std::isfinite(d3->at(a))) f.dedxP=d3->at(a);
      if(haveDedx&&d3&&b<(int)d3->size()&&std::isfinite(d3->at(b))) f.dedxPi=d3->at(b);
      const double s=provisionalPairScore(f);
      if(s>bestA){ bestA=s; prA=a; piA=b; }
      if(!passPrecuts(f,pc)) continue;
      if(s>best){ best=s; pr=a; pi=b; } } }
    precutOK = pr>=0; if(pr<0){ pr=prA; pi=piA; } if(pr<0) continue;
    const DecayObs obs=decayObservables(vtx, st[pr],dir[pr],H.trk_len[pr], st[pi],dir[pi],H.trk_len[pi], protonMomFromRange(H.trk_len[pr]), H.trk_rmumom[pi]);
    if(!obs.ok) continue;
    mass=obs.mass; decayLen=obs.decayLength; cosAlpha=obs.cosAlpha;
    if(!(std::fabs(mass-M_LAMBDA)<PRE_WIN && decayLen>PRE_LEN && cosAlpha>PRE_COS)) continue;
    if(!(vw && vw->size()>2)) continue;
    evt=H.evt; trueFlight = isSig ? std::sqrt((lex-lsx)*(lex-lsx)+(ley-lsy)*(ley-lsy)+(lez-lsz)*(lez-lsz)) : -1;
    const int cand[3]={mu,pr,pi};
    for(int k=0;k<3;++k) seedT[k]=Island::tickOf(st[cand[k]].X());
    for(int p=0;p<3;++p){
      hw[p].clear(); ht[p].clear(); ha[p].clear(); hp[p].clear();
      firstWire[p] = (win[p]&&!win[p]->empty()&&vw->at(p)>=0) ? vw->at(p)-(int)win[p]->size()/2 : -1;
      for(int k=0;k<3;++k) seedW[p][k]=Island::wireOf(p, st[cand[k]].Y(), st[cand[k]].Z());
      if(firstWire[p]<0) continue;
      for(int i=0;i<(int)win[p]->size();++i) for(size_t k=0;k<(*win[p])[i].size();++k){
        hw[p].push_back(firstWire[p]+i); ht[p].push_back((*win[p])[i][k]);
        ha[p].push_back((adc[p]&&i<(int)adc[p]->size()&&k<(*adc[p])[i].size())?(*adc[p])[i][k]:-1.f);
        hp[p].push_back((pdg[p]&&i<(int)pdg[p]->size()&&k<(*pdg[p])[i].size())?(*pdg[p])[i][k]:0); }
    }
    T.Fill(); nWritten++;
  }
  fout.cd(); T.Write(); printf("  wrote %ld events to %s\n", nWritten, outName); fflush(stdout);
}
