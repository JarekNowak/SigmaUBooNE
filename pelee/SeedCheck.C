// SeedCheck.C -- do the muon/proton/pion seeds land on the right particle? Uses the
// per-hit truth label of the occupied cell nearest each seed. Separates a seeding
// (transform) failure from a genuine geometry limit.
#include "HypInput.h"
#include "PidAssign.h"
#include "SigmaFV.h"
#include "Island.h"
#include <fstream>
#include <map>
static const double MUON_TRK_SCORE=0.8,MUON_TRK_DIST=20.0,MUON_TRK_LEN=10.0,MUON_LLR=0.7;
static const int NSHOWERS_MAX=4,NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5,TOPO_MIN=0.20;
void SeedCheck(int maxFiles=700){
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int nf=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||nf<maxFiles)){ch.Add(l.c_str());nf++;}
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  int fm=-1; if(ch.GetBranch("flash_match_nu_slice_ID")) ch.SetBranchAddress("flash_match_nu_slice_ID",&fm);
  std::vector<std::vector<float>>* win[3]={0,0,0}; std::vector<std::vector<int>>* pdg[3]={0,0,0};
  std::vector<int>* vw=0;
  for(int p=0;p<3;++p){ ch.SetBranchAddress(Form("ct_test_window_plane%d",p),&win[p]);
                        ch.SetBranchAddress(Form("ct_test_pdg_plane%d",p),&pdg[p]); }
  ch.SetBranchAddress("ct_test_primary_vtx_wires",&vw);
  std::vector<int>* tpdg=0; if(ch.GetBranch("pfp_true_pdg")) ch.SetBranchAddress("pfp_true_pdg",&tpdg);
  // per role: seeds tried, seed found a cell, nearest cell's truth == candidate's truth
  long tried[3]={0,0,0}, found[3]={0,0,0}, match[3]={0,0,0}; std::vector<double> dist[3];
  const char* ROLE[3]={"muon","proton","pion"};
  for(Long64_t ie=0;ie<ch.GetEntries();++ie){ ch.GetEntry(ie); H.convert();
    if(!(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0)) continue;
    if(ch.GetBranch("flash_match_nu_slice_ID")&&fm<0) continue;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int nn=(int)H.trk_score.size(); if(nn<1||!tpdg) continue;
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int k){return H.trk_score[k]>TRK_SHR_SPLIT;};
    std::vector<double> dv(nn,1e9); std::vector<TVector3> st(nn);
    for(int k=0;k<nn;++k){ if(isTrk(k)) st[k].SetXYZ(H.trk_sx[k],H.trk_sy[k],H.trk_sz[k]);
      else st[k].SetXYZ(H.shr_sx[k],H.shr_sy[k],H.shr_sz[k]); dv[k]=(st[k]-vtx).Mag(); }
    int ntr=0,nsh=0; for(int k=0;k<nn;++k)(isTrk(k)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;
    int mu=-1; for(int k=0;k<nn;++k)
      if(isTrk(k)&&H.trk_score[k]>MUON_TRK_SCORE&&dv[k]<MUON_TRK_DIST&&H.trk_len[k]>MUON_TRK_LEN&&H.trk_llr[k]>MUON_LLR)
        if(mu==-1||H.trk_len[k]>H.trk_len[mu]) mu=k;
    if(mu<0) continue; if(!(ntr<NTRACKS_MAX)) continue;
    std::vector<int> cand; for(int k=0;k<nn;++k) if(k!=mu&&isTrk(k)) cand.push_back(k);
    std::vector<double> dvv(nn),nhv(nn);
    for(int k=0;k<nn;++k){ dvv[k]=dv[k]; nhv[k]=(k<(int)H.pfnplanehits_Y.size())?H.pfnplanehits_Y[k]:0; }
    int pr=-1,pi=-1;
    if(!PidAssign::assign(cand,H.trk_score,H.trk_llr,H.trk_chipr,H.trk_chipi,H.trk_len,dvv,nhv,pr,pi)) continue;
    if(!vw||vw->size()<3) continue;
    const int idx[3]={mu,pr,pi};
    for(int p=0;p<3;++p){
      if(!win[p]||win[p]->empty()||!pdg[p]||vw->at(p)<0) continue;
      const int nw=win[p]->size(), w0=vw->at(p)-nw/2;
      for(int r=0;r<3;++r){
        const int k=idx[r]; if(k<0||k>=(int)tpdg->size()) continue;
        const int truthOfCand=tpdg->at(k);
        Island::Cell sd{Island::wireOf(p,st[k].Y(),st[k].Z()),Island::tickOf(st[k].X())};
        double best=1e18; int bw=-1,bh=-1;
        for(int w=0;w<nw;++w) for(size_t h=0;h<(*win[p])[w].size();++h){
          double dw=(w0+w)-sd.w, dt=((*win[p])[w][h]-sd.t)*0.055/0.3, d=dw*dw+dt*dt;
          if(d<best){best=d;bw=w;bh=h;} }
        tried[r]++;
        if(bw<0||best>=400.0) continue;
        found[r]++; dist[r].push_back(std::sqrt(best)*0.3);   // ~cm
        const int hitPdg = (bw<(int)pdg[p]->size()&&bh<(int)(*pdg[p])[bw].size()) ? (*pdg[p])[bw][bh] : 0;
        if(hitPdg==truthOfCand) match[r]++;
      }
    }
  }
  printf("\n  seed placement, signal, all planes pooled\n");
  printf("  %-7s %8s %10s %14s %12s\n","role","seeds","found cell","truth agrees","median dist");
  for(int r=0;r<3;++r){ std::sort(dist[r].begin(),dist[r].end());
    printf("  %-7s %8ld %9.1f%% %13.1f%% %9.2f cm\n",ROLE[r],tried[r],
      tried[r]?100.*found[r]/tried[r]:0., found[r]?100.*match[r]/found[r]:0.,
      dist[r].empty()?0.:dist[r][dist[r].size()/2]); }
  printf("  ('truth agrees' = nearest hit's truth label equals the candidate track's truth label)\n");
}
