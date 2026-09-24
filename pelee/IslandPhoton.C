// IslandPhoton.C -- is the Sigma0 photon bridging the muon and Lambda islands?
// Run the seeded island test on signal three ways: with every hit in the window,
// with hits truth-labelled as photons removed, and with all EM hits (22, +-11)
// removed. Bin by true Lambda flight distance. If the photon bridges, removing it
// lifts the pass rate in the long-flight bins where the Lambda genuinely separates.
#include "HypInput.h"
#include "PidAssign.h"
#include "SigmaFV.h"
#include "Island.h"
#include <fstream>
#include <map>
static const double MUON_TRK_SCORE=0.8,MUON_TRK_DIST=20.0,MUON_TRK_LEN=10.0,MUON_LLR=0.7;
static const int NSHOWERS_MAX=4,NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5,TOPO_MIN=0.20;

static bool testGrid(const std::vector<std::vector<float>>& win,
                     const std::vector<std::vector<int>>* pdg, int mask, // 0 none,1 photon,2 EM
                     int w0,int plane,const TVector3& m,const TVector3& p,const TVector3& i,bool& ok){
  std::vector<std::vector<float>> filt(win.size());
  for(size_t w=0;w<win.size();++w) for(size_t h=0;h<win[w].size();++h){
    int pd = (pdg && w<pdg->size() && h<(*pdg)[w].size()) ? (*pdg)[w][h] : 0;
    bool drop = (mask==1 && pd==22) || (mask==2 && (pd==22||pd==11||pd==-11));
    if(!drop) filt[w].push_back(win[w][h]); }
  return Island::planeTest(filt,w0,plane,m.X(),m.Y(),m.Z(),p.X(),p.Y(),p.Z(),i.X(),i.Y(),i.Z(),ok);
}

void IslandPhoton(int maxFiles=700){
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
  double lsx=0,lsy=0,lsz=0,lex=0,ley=0,lez=0;
  auto B=[&](const char*n,void*q){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,q); };
  B("mc_sigmazero_lambda_start_x",&lsx);B("mc_sigmazero_lambda_start_y",&lsy);B("mc_sigmazero_lambda_start_z",&lsz);
  B("mc_sigmazero_lambda_end_x",&lex);  B("mc_sigmazero_lambda_end_y",&ley);  B("mc_sigmazero_lambda_end_z",&lez);

  const char* BIN[5]={"<1 cm","1-2 cm","2-5 cm","5-10 cm",">10 cm"};
  long n[5]={0}, pass[3][5]={{0}};
  std::map<int,long> pdgHist;
  for(Long64_t ie=0;ie<ch.GetEntries();++ie){ ch.GetEntry(ie); H.convert();
    if(!(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0)) continue;
    if(ch.GetBranch("flash_match_nu_slice_ID")&&fm<0) continue;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int nn=(int)H.trk_score.size(); if(nn<1) continue;
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
    const double fl=std::sqrt((lex-lsx)*(lex-lsx)+(ley-lsy)*(ley-lsy)+(lez-lsz)*(lez-lsz));
    int b = fl<1?0:(fl<2?1:(fl<5?2:(fl<10?3:4)));
    n[b]++;
    // what truth labels sit in the window? (plane 2 only, for the histogram)
    if(pdg[2]) for(auto& col:*pdg[2]) for(int pd:col) pdgHist[pd]++;
    for(int mask=0;mask<3;++mask){
      bool any=false;
      for(int p=0;p<3;++p){
        if(!win[p]||win[p]->empty()||vw->at(p)<0) continue;
        bool ok=false;
        if(testGrid(*win[p],pdg[p],mask,vw->at(p)-(int)win[p]->size()/2,p,st[mu],st[pr],st[pi],ok)) any=true; }
      if(any) pass[mask][b]++; }
  }
  printf("\n  truth labels of hits in the plane-2 window (top 8):\n");
  { std::vector<std::pair<long,int>> v; for(auto&kv:pdgHist) v.push_back({kv.second,kv.first});
    std::sort(v.rbegin(),v.rend()); long tot=0; for(auto&x:v) tot+=x.first;
    for(size_t k=0;k<v.size()&&k<8;++k) printf("    pdg %6d : %5.1f%%\n",v[k].second,100.*v[k].first/tot); }
  printf("\n  island pass rate (>=1 plane) by true Lambda flight, signal\n");
  printf("  %-9s %6s %12s %14s %14s\n","flight","n","all hits","no photon(22)","no EM(22,+-11)");
  long N=0,P[3]={0,0,0};
  for(int b=0;b<5;++b){ N+=n[b]; for(int m=0;m<3;++m) P[m]+=pass[m][b];
    printf("  %-9s %6ld %11.1f%% %13.1f%% %13.1f%%\n",BIN[b],n[b],
      n[b]?100.*pass[0][b]/n[b]:0.,n[b]?100.*pass[1][b]/n[b]:0.,n[b]?100.*pass[2][b]/n[b]:0.); }
  printf("  %-9s %6ld %11.1f%% %13.1f%% %13.1f%%\n","ALL",N,100.*P[0]/N,100.*P[1]/N,100.*P[2]/N);
}
