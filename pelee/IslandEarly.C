// IslandEarly.C -- can the island test act as a CUT earlier in the chain?
// The full test needs proton and pion seeds, so it cannot precede the candidate
// stages. The only form that can is MUON-SEEDED: cluster the activity, take the
// muon's cluster, and ask whether a second cluster of real size exists that is
// NOT connected to it -- i.e. is there displaced activity at all. That is
// applicable immediately after MuonCandidate, with no pair required.
#include "HypInput.h"
#include "Island.h"
#include "SigmaFV.h"
#include <fstream>
static const double MUON_TRK_SCORE=0.8,MUON_TRK_DIST=20.0,MUON_TRK_LEN=10.0,MUON_LLR=0.7;
static const int NSHOWERS_MAX=4,NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5,TOPO_MIN=0.20;

void IslandEarly(const char* list, int maxFiles, bool signalOnly){
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||nf<maxFiles)){ch.Add(l.c_str());nf++;}
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  int fm=-1; if(ch.GetBranch("flash_match_nu_slice_ID")) ch.SetBranchAddress("flash_match_nu_slice_ID",&fm);
  std::vector<std::vector<float>>* win[3]={0,0,0}; std::vector<int>* vw=0;
  ch.SetBranchAddress("ct_test_window_plane0",&win[0]);
  ch.SetBranchAddress("ct_test_window_plane1",&win[1]);
  ch.SetBranchAddress("ct_test_window_plane2",&win[2]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires",&vw);
  long nMu=0; long pass[4]={0,0,0,0};   // >=1 plane, for min sizes 3/5/10/20
  const int MS[4]={3,5,10,20};
  for(Long64_t ie=0;ie<ch.GetEntries();++ie){ ch.GetEntry(ie); H.convert();
    const bool isSig=(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0);
    if(signalOnly!=isSig) continue;
    if(ch.GetBranch("flash_match_nu_slice_ID")&&fm<0) continue;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int n=(int)H.trk_score.size(); if(n<1) continue;
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int i){return H.trk_score[i]>TRK_SHR_SPLIT;};
    std::vector<double> dv(n,1e9); std::vector<TVector3> st(n);
    for(int i=0;i<n;++i){ if(isTrk(i)) st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]); dv[i]=(st[i]-vtx).Mag(); }
    int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;
    int mu=-1; for(int i=0;i<n;++i)
      if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dv[i]<MUON_TRK_DIST&&
         H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
        if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
    if(mu<0) continue;                         // <- stage: MuonCandidate
    if(!vw||vw->size()<3) continue;
    nMu++;
    bool ok[4]={false,false,false,false};
    for(int p=0;p<3;++p){
      if(!win[p]||win[p]->empty()||vw->at(p)<0) continue;
      const int nw=win[p]->size(), w0=vw->at(p)-nw/2;
      std::vector<Island::Cell> pts;
      for(int i=0;i<nw;++i) for(float t:(*win[p])[i]) pts.push_back({w0+i,t});
      const int N=pts.size(); if(N<3) continue;
      auto conn=[&](const Island::Cell&a,const Island::Cell&b){
        return std::fabs(a.t-b.t)<=Island::CDT && std::abs(a.w-b.w)<=Island::CDW; };
      std::vector<int> lab(N,-1); int nl=0;
      for(int i=0;i<N;++i){ if(lab[i]>=0) continue; std::vector<int> s{i}; lab[i]=nl;
        while(!s.empty()){ int q=s.back(); s.pop_back();
          for(int j=0;j<N;++j) if(lab[j]<0&&conn(pts[j],pts[q])){lab[j]=nl;s.push_back(j);} } ++nl; }
      std::vector<int> csz(nl,0); for(int i=0;i<N;++i) csz[lab[i]]++;
      // muon cluster
      Island::Cell sd{Island::wireOf(p,st[mu].Y(),st[mu].Z()),Island::tickOf(st[mu].X())};
      double best=1e18; int bi=-1;
      for(int i=0;i<N;++i){ double dw=pts[i].w-sd.w, dt=(pts[i].t-sd.t)*0.055/0.3;
        double d=dw*dw+dt*dt; if(d<best){best=d;bi=i;} }
      if(bi<0||best>=400.0) continue;
      const int muLab=lab[bi];
      for(int m=0;m<4;++m)
        for(int c=0;c<nl;++c) if(c!=muLab && csz[c]>=MS[m]) { ok[m]=true; break; }
    }
    for(int m=0;m<4;++m) if(ok[m]) pass[m]++;
  }
  printf("\n  %s : %ld events at MuonCandidate\n", signalOnly?"SIGNAL":"BACKGROUND", nMu);
  for(int m=0;m<4;++m)
    printf("    island off the muon, size>=%2d : %6ld  (%.1f%%)\n",MS[m],pass[m],nMu?100.*pass[m]/nMu:0.);
}
