// ============================================================================
//  IslandVeto.C -- does removing the photon's activity recover the Lambda CCQE
//  island performance?  Strip the photon and a Sigma0 event IS a Lambda CCQE
//  event, so the seeded test should behave like theirs.
//
//  Their technote (sec 7.3) says the test is "highly efficient for true Lambda
//  ranges greater than 2 cm". 55.9% of our Lambdas exceed 2 cm, so matching them
//  means ~50-56% pass; we measure 27.2%. This asks whether the photon explains
//  the factor of two.
//
//  Masking is RECO-ONLY -- no truth branches -- both because it is what we would
//  actually deploy and because ct_test_pdg_plane* is so large that truth-based
//  masking is I/O-bound to the point of not finishing (4 h elapsed for 9 min of
//  CPU before it was killed). Each occupied cell is assigned to whichever of
//  {muon, proton, pion, shower1} seed is nearest in wire-time; mode 1 drops the
//  cells assigned to the shower, mode 2 instead drops cells within a fixed
//  radius of the shower seed.
//
//  Run:  root -l -b -q 'IslandVeto.C("lists/sig_genie.list",700,true)'
// ============================================================================
#include "HypInput.h"
#include "PidAssign.h"
#include "SigmaFV.h"
#include "Island.h"
#include <fstream>

static const double MUON_TRK_SCORE=0.8, MUON_TRK_DIST=20.0, MUON_TRK_LEN=10.0, MUON_LLR=0.7;
static const int    NSHOWERS_MAX=4, NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5, TOPO_MIN=0.20;

void IslandVeto(const char* list="lists/sig_genie.list", int maxFiles=700, bool signalOnly=true)
{
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf=0;
  while(std::getline(in,l)) if(l.size() && (maxFiles<=0||nf<maxFiles)){ ch.Add(l.c_str()); nf++; }
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){ printf("  not hypana\n"); return; }
  int fm=-1; if(ch.GetBranch("flash_match_nu_slice_ID"))
    ch.SetBranchAddress("flash_match_nu_slice_ID",&fm);
  std::vector<std::vector<float>>* win[3]={0,0,0};
  std::vector<int>* vw=0;
  for(int p=0;p<3;++p) ch.SetBranchAddress(Form("ct_test_window_plane%d",p),&win[p]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires",&vw);
  double lsx=0,lsy=0,lsz=0,lex=0,ley=0,lez=0;
  auto B=[&](const char* n,void* q){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,q); };
  B("mc_sigmazero_lambda_start_x",&lsx); B("mc_sigmazero_lambda_start_y",&lsy);
  B("mc_sigmazero_lambda_start_z",&lsz); B("mc_sigmazero_lambda_end_x",&lex);
  B("mc_sigmazero_lambda_end_y",&ley);   B("mc_sigmazero_lambda_end_z",&lez);

  const char* BIN[4]={"<2 cm","2-5 cm","5-10 cm",">10 cm"};
  long n[4]={0,0,0,0}, pass[3][4]={{0}};
  long nNoShr=0;

  for(Long64_t ie=0; ie<ch.GetEntries(); ++ie){
    ch.GetEntry(ie); H.convert();
    const bool isSig=(H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
    if(signalOnly != isSig) continue;
    if(ch.GetBranch("flash_match_nu_slice_ID") && fm<0) continue;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int nn=(int)H.trk_score.size(); if(nn<1) continue;
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int k){ return H.trk_score[k]>TRK_SHR_SPLIT; };
    std::vector<double> dv(nn,1e9); std::vector<TVector3> st(nn);
    for(int k=0;k<nn;++k){
      if(isTrk(k)) st[k].SetXYZ(H.trk_sx[k],H.trk_sy[k],H.trk_sz[k]);
      else         st[k].SetXYZ(H.shr_sx[k],H.shr_sy[k],H.shr_sz[k]);
      dv[k]=(st[k]-vtx).Mag(); }
    int ntr=0,nsh=0; for(int k=0;k<nn;++k)(isTrk(k)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;
    int mu=-1; for(int k=0;k<nn;++k)
      if(isTrk(k)&&H.trk_score[k]>MUON_TRK_SCORE&&dv[k]<MUON_TRK_DIST&&
         H.trk_len[k]>MUON_TRK_LEN&&H.trk_llr[k]>MUON_LLR)
        if(mu==-1||H.trk_len[k]>H.trk_len[mu]) mu=k;
    if(mu<0) continue;
    if(!(ntr<NTRACKS_MAX)) continue;
    std::vector<int> cand; for(int k=0;k<nn;++k) if(k!=mu&&isTrk(k)) cand.push_back(k);
    std::vector<double> dvv(nn),nhv(nn);
    for(int k=0;k<nn;++k){ dvv[k]=dv[k];
      nhv[k]=(k<(int)H.pfnplanehits_Y.size())?H.pfnplanehits_Y[k]:0; }
    int pr=-1,pi=-1;
    if(!PidAssign::assign(cand,H.trk_score,H.trk_llr,H.trk_chipr,H.trk_chipi,
                          H.trk_len,dvv,nhv,pr,pi)) continue;
    int s1=-1; double b1=1e9;
    for(int k=0;k<nn;++k){ if(k==mu||k==pr||k==pi||isTrk(k)) continue;
      if(dv[k]<b1){ b1=dv[k]; s1=k; } }
    if(!vw||vw->size()<3) continue;
    if(s1<0) nNoShr++;

    const double fl=std::sqrt((lex-lsx)*(lex-lsx)+(ley-lsy)*(ley-lsy)+(lez-lsz)*(lez-lsz));
    const int b = fl<2?0:(fl<5?1:(fl<10?2:3));
    n[b]++;

    for(int mode=0; mode<3; ++mode){
      bool any=false;
      for(int p=0;p<3;++p){
        if(!win[p]||win[p]->empty()||vw->at(p)<0) continue;
        const int w0=vw->at(p)-(int)win[p]->size()/2;
        const Island::Cell sMu{Island::wireOf(p,st[mu].Y(),st[mu].Z()),Island::tickOf(st[mu].X())};
        const Island::Cell sPr{Island::wireOf(p,st[pr].Y(),st[pr].Z()),Island::tickOf(st[pr].X())};
        const Island::Cell sPi{Island::wireOf(p,st[pi].Y(),st[pi].Z()),Island::tickOf(st[pi].X())};
        Island::Cell sSh{0,0}; const bool haveShr=(s1>=0);
        if(haveShr) sSh=Island::Cell{Island::wireOf(p,st[s1].Y(),st[s1].Z()),Island::tickOf(st[s1].X())};
        std::vector<std::vector<float>> filt(win[p]->size());
        for(size_t w=0;w<win[p]->size();++w) for(float t:(*win[p])[w]){
          bool drop=false;
          if(mode>0 && haveShr){
            const int aw=w0+(int)w;
            auto d2=[&](const Island::Cell& c){
              const double dw=aw-c.w, dt=(t-c.t)*0.055/0.3; return dw*dw+dt*dt; };
            const double ds=d2(sSh);
            if(mode==1) drop = (ds<d2(sMu) && ds<d2(sPr) && ds<d2(sPi));
            else        drop = (ds<25.0);            // within ~5 wires of the shower seed
          }
          if(!drop) filt[w].push_back(t);
        }
        bool ok=false;
        if(Island::planeTest(filt,w0,p,
                             st[mu].X(),st[mu].Y(),st[mu].Z(),
                             st[pr].X(),st[pr].Y(),st[pr].Z(),
                             st[pi].X(),st[pi].Y(),st[pi].Z(), ok)) any=true;
      }
      if(any) pass[mode][b]++;
    }
  }

  printf("\n  %s : island pass rate (>=1 plane) by true Lambda flight\n",
         signalOnly?"SIGNAL":"BACKGROUND");
  printf("  %-9s %6s %11s %14s %15s\n","flight","n","all hits","drop shower","shower+radius");
  long N=0,P[3]={0,0,0};
  for(int b=0;b<4;++b){ N+=n[b]; for(int m=0;m<3;++m) P[m]+=pass[m][b];
    printf("  %-9s %6ld %10.1f%% %13.1f%% %14.1f%%\n",BIN[b],n[b],
      n[b]?100.*pass[0][b]/n[b]:0., n[b]?100.*pass[1][b]/n[b]:0., n[b]?100.*pass[2][b]/n[b]:0.); }
  printf("  %-9s %6ld %10.1f%% %13.1f%% %14.1f%%\n","ALL",N,
    N?100.*P[0]/N:0., N?100.*P[1]/N:0., N?100.*P[2]/N:0.);
  printf("  (events with no shower candidate: %ld)\n",nNoShr);
  fflush(stdout);
}
