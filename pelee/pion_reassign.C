// pion_reassign.C -- the pion candidate's biggest contaminant is photons (25%).
// A track-score CUT trades signal for purity (FoM falls monotonically). Choosing
// a BETTER candidate among the PFPs already present costs nothing: the event is
// kept either way, only the index changes. Test adding photon-rejection terms --
// wiggliness (true pi- 0.0253 vs true gamma 0.0442) and dE/dx -- to the ranking.
#include "HypInput.h"
#include "SigmaFV.h"
#include "TVector3.h"
#include <fstream>
static const double MUON_TRK_SCORE=0.8,MUON_TRK_DIST=20.0,MUON_TRK_LEN=10.0,MUON_LLR=0.7;
static const int NSHOWERS_MAX=4,NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5,TOPO_MIN=0.20,LLR_SCALE=0.591,DCHI_SCALE=124.5;

void pion_reassign(const char* list="/home/t2k/nowak/.claude/jobs/a2cc9f4f/tmp/nuwro_hypana.list"){
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; while(std::getline(in,l)) if(l.size()) ch.Add(l.c_str());
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){printf(" not hypana\n");return;}
  // must bind AFTER attach(): on a TChain, attach() re-binds the branch list and
  // silently drops addresses set before it, leaving these vectors empty.
  std::vector<double> *rms=nullptr,*d3=nullptr;
  ch.SetBranchAddress("trk_rms_deflection",&rms);
  ch.SetBranchAddress("trk_three_plane_dedx",&d3);

  struct Cfg{ const char* name; double wW, wD; };
  std::vector<Cfg> cfg = {
    {"current D only",              0.0,  0.0},
    {"D - 2*wiggle",                2.0,  0.0},
    {"D - 5*wiggle",                5.0,  0.0},
    {"D - 10*wiggle",              10.0,  0.0},
    {"D - 20*wiggle",              20.0,  0.0},
    {"D - 10*wiggle - 0.1*dedx",   10.0,  0.1},
    {"D - 10*wiggle - 0.3*dedx",   10.0,  0.3},
  };
  const int NC=cfg.size();
  std::vector<long> sig(NC,0),pok(NC,0),pgam(NC,0),bkg(NC,0);

  for(Long64_t ie=0;ie<ch.GetEntries();++ie){
    ch.GetEntry(ie); H.convert();
    const bool isSig=(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0);
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z))continue;
    if(topo<TOPO_MIN)continue;
    const int n=(int)H.trk_score.size();
    std::vector<TVector3> st(n); std::vector<double> dv(n,1e9);
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int i){return H.trk_score[i]>TRK_SHR_SPLIT;};
    for(int i=0;i<n;++i){ if(isTrk(i))st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]); dv[i]=(st[i]-vtx).Mag(); }
    int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX))continue;
    if(!(ntr<NTRACKS_MAX))continue;
    int mu=-1; for(int i=0;i<n;++i)
      if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dv[i]<MUON_TRK_DIST&&H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
        if(mu==-1||H.trk_len[i]>H.trk_len[mu])mu=i;
    if(mu<0)continue;
    auto D=[&](int i){ double d=H.trk_llr[i]/LLR_SCALE; double dc=H.trk_chipr[i]-H.trk_chipi[i];
      if(std::isfinite(dc)&&H.trk_chipr[i]>-900&&H.trk_chipi[i]>-900)d+=dc/DCHI_SCALE; return d; };
    int pr=-1; for(int i=0;i<n;++i){ if(i==mu||!isTrk(i))continue; if(pr==-1||D(i)<D(pr))pr=i; }
    if(pr<0)continue;
    for(int c=0;c<NC;++c){
      auto Dp=[&](int i){ double d=D(i);
        if(cfg[c].wW>0&&rms&&i<(int)rms->size()&&std::isfinite(rms->at(i))) d-=cfg[c].wW*rms->at(i);
        if(cfg[c].wD>0&&d3 &&i<(int)d3->size() &&std::isfinite(d3->at(i)) &&d3->at(i)>0) d-=cfg[c].wD*d3->at(i);
        return d; };
      int pi=-1; for(int i=0;i<n;++i){ if(i==mu||i==pr||!isTrk(i))continue;
        if(pi==-1||Dp(i)>Dp(pi))pi=i; }
      if(pi<0)continue;
      if(isSig){ sig[c]++; if(H.bkt_pdg[pi]==-211)pok[c]++; if(H.bkt_pdg[pi]==22)pgam[c]++; }
      else bkg[c]++;
    }
  }
  printf("\n  === pion-candidate RANKING (no efficiency cost: same events kept) ===\n\n");
  printf("  %-28s %7s %7s %10s %9s %8s\n","ranking","signal","bkg","pi- purity","gamma","FoM rel");
  double base=sig[0]/std::sqrt((double)sig[0]+bkg[0]);
  for(int c=0;c<NC;++c){
    if(sig[c]+bkg[c]<=0)continue;
    double fom=sig[c]/std::sqrt((double)sig[c]+bkg[c]);
    printf("  %-28s %7ld %7ld %9.1f%% %8.1f%% %8.3f\n",cfg[c].name,sig[c],bkg[c],
      sig[c]?100.*pok[c]/sig[c]:0., sig[c]?100.*pgam[c]/sig[c]:0., fom/base);
  }
}
