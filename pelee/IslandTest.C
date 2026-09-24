// ============================================================================
//  IslandTest.C -- the connectedness test as the Lambda CCQE analysis defines it,
//  rather than the generic island statistics our BDT currently carries.
//
//  Their version (Hyperon_Technote_v2 sec. 7): build a wire x time grid of
//  filtered activity, SEED clusters at the muon / proton / pion track starts,
//  grow them, and require the proton and pion clusters to MERGE while the muon
//  cluster stays SEPARATE -- on at least one plane. Our ct_* variables instead
//  count islands in a window around the PRIMARY VERTEX with no reference to which
//  island holds which candidate, which is why the ablation finds them worthless.
//
//  Coordinate transforms derived and validated here (ctcalib.C):
//     wire = (A*y + B*z + C)/0.3   A,B,C per plane, 98.7% of stored vertex
//                                  wires reproduced within 2 wires
//     tick = 18.1509*x + 818.1     implied drift 0.0551 cm/tick, matching the
//                                  known 0.1098 cm/us at 0.5 us/tick
//
//  Run:  root -l -b -q 'IslandTest.C("lists/sig_genie.list",400)'
// ============================================================================
#include "HypInput.h"
#include "PidAssign.h"
#include "SigmaFV.h"
#include "TChain.h"
#include "TVector3.h"
#include <fstream>
#include <vector>
#include <cmath>

static const double MUON_TRK_SCORE=0.8, MUON_TRK_DIST=20.0, MUON_TRK_LEN=10.0, MUON_LLR=0.7;
static const int    NSHOWERS_MAX=4, NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5, TOPO_MIN=0.20;
// grid connectivity: two occupied cells join if within CDW wires and CDT ticks
static const int    CDW=3;
static const double CDT=40.0;
static const int    MINCLUST=3;     // quality cut: minimum cluster size

// wire = (A*y + B*z + C)/0.3 ; tick = MT*x + CT
static const double WA[3]={-0.86603,0.86603,0.0}, WB[3]={0.5,0.5,1.0}, WC[3]={101.475,99.748,0.0};
static const double MT=18.1509, CT=818.1;

struct Cell { int w; float t; };

// Grow clusters from seeds; returns a label per seed (-1 if the seed found no
// activity). Labels are shared when seeds land in the same connected component.
static void clusterFromSeeds(const std::vector<Cell>& pts,
                             const std::vector<Cell>& seeds,
                             std::vector<int>& seedLabel,
                             std::vector<int>& clustSize)
{
  const int n=pts.size();
  seedLabel.assign(seeds.size(),-1); clustSize.clear();
  if(n==0) return;
  auto conn=[&](const Cell&a,const Cell&b){
    return std::fabs(a.t-b.t)<=CDT && std::abs(a.w-b.w)<=CDW; };
  std::vector<int> lab(n,-1); int nl=0;
  for(int i=0;i<n;++i){ if(lab[i]>=0) continue;
    std::vector<int> st{i}; lab[i]=nl;
    while(!st.empty()){ int q=st.back(); st.pop_back();
      for(int j=0;j<n;++j) if(lab[j]<0 && conn(pts[j],pts[q])){ lab[j]=nl; st.push_back(j);} }
    ++nl; }
  clustSize.assign(nl,0);
  for(int i=0;i<n;++i) clustSize[lab[i]]++;
  // attach each seed to the nearest occupied cell
  for(size_t s=0;s<seeds.size();++s){
    double best=1e18; int bi=-1;
    for(int i=0;i<n;++i){
      double dw=(pts[i].w-seeds[s].w), dt=(pts[i].t-seeds[s].t)/MT*0.3/0.3;
      double d=dw*dw + (dt*0.055/0.3)*(dt*0.055/0.3);   // wires vs ~cm-equivalent ticks
      if(d<best){best=d;bi=i;} }
    if(bi>=0 && best < 400.0) seedLabel[s]=lab[bi];     // within ~20 wires
  }
}

void IslandTest(const char* list="lists/sig_genie.list", int maxFiles=400,
                bool signalOnly=true)
{
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf=0;
  while(std::getline(in,l)) if(l.size() && (maxFiles<=0||nf<maxFiles)){ ch.Add(l.c_str()); nf++; }
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){ printf("  not hypana\n"); return; }
  int fmSlice=-1; if(ch.GetBranch("flash_match_nu_slice_ID"))
    ch.SetBranchAddress("flash_match_nu_slice_ID",&fmSlice);
  std::vector<std::vector<float>> *win[3]={0,0,0};
  std::vector<int>* vw=0;
  ch.SetBranchAddress("ct_test_window_plane0",&win[0]);
  ch.SetBranchAddress("ct_test_window_plane1",&win[1]);
  ch.SetBranchAddress("ct_test_window_plane2",&win[2]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires",&vw);

  long nEv=0, nPassAny=0, nPlaneOK[3]={0,0,0}, nSeeded=0;
  const Long64_t N=ch.GetEntries();
  printf("  %d files, %lld events, signalOnly=%d\n",nf,N,(int)signalOnly); fflush(stdout);

  for(Long64_t ie=0; ie<N; ++ie){
    ch.GetEntry(ie); H.convert();
    const bool isSig=(H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
    if(signalOnly != isSig) continue;
    if(ch.GetBranch("flash_match_nu_slice_ID") && fmSlice<0) continue;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int n=(int)H.trk_score.size(); if(n<1) continue;
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int i){ return H.trk_score[i]>TRK_SHR_SPLIT; };
    std::vector<double> dv(n,1e9); std::vector<TVector3> st(n);
    for(int i=0;i<n;++i){
      if(isTrk(i)) st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else         st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]);
      dv[i]=(st[i]-vtx).Mag(); }
    int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;
    int mu=-1; for(int i=0;i<n;++i)
      if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dv[i]<MUON_TRK_DIST&&
         H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
        if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
    if(mu<0) continue;
    if(!(ntr<NTRACKS_MAX)) continue;
    std::vector<int> cand; for(int i=0;i<n;++i) if(i!=mu&&isTrk(i)) cand.push_back(i);
    std::vector<double> dvv(n),nhv(n);
    for(int i=0;i<n;++i){ dvv[i]=dv[i];
      nhv[i]=(i<(int)H.pfnplanehits_Y.size())?H.pfnplanehits_Y[i]:0; }
    int pr=-1,pi=-1;
    if(!PidAssign::assign(cand,H.trk_score,H.trk_llr,H.trk_chipr,H.trk_chipi,
                          H.trk_len,dvv,nhv,pr,pi)) continue;
    if(!vw||vw->size()<3) continue;
    nEv++;

    bool anyPlane=false, seededSomewhere=false;
    for(int p=0;p<3;++p){
      if(!win[p]||win[p]->empty()) continue;
      const int nw=win[p]->size(), w0=vw->at(p)-nw/2;
      std::vector<Cell> pts;
      for(int i=0;i<nw;++i) for(float t:(*win[p])[i]) pts.push_back({w0+i,t});
      if((int)pts.size()<MINCLUST) continue;
      auto seedOf=[&](int idx){
        Cell c; c.w=(int)std::lround((WA[p]*st[idx].Y()+WB[p]*st[idx].Z()+WC[p])/0.3);
        c.t=(float)(MT*st[idx].X()+CT); return c; };
      std::vector<Cell> seeds={seedOf(mu),seedOf(pr),seedOf(pi)};
      std::vector<int> lab, csz;
      clusterFromSeeds(pts,seeds,lab,csz);
      if(lab[0]<0||lab[1]<0||lab[2]<0) continue;
      seededSomewhere=true;
      // quality: clusters must be big enough
      if(csz[lab[0]]<MINCLUST||csz[lab[1]]<MINCLUST||csz[lab[2]]<MINCLUST) continue;
      const bool pPiMerged = (lab[1]==lab[2]);
      const bool muSeparate= (lab[0]!=lab[1] && lab[0]!=lab[2]);
      if(pPiMerged && muSeparate){ nPlaneOK[p]++; anyPlane=true; }
    }
    if(seededSomewhere) nSeeded++;
    if(anyPlane) nPassAny++;
  }
  printf("\n  events at preselection with a p/pi pair : %ld\n",nEv);
  printf("  all three seeds found activity           : %ld  (%.1f%%)\n",nSeeded,nEv?100.*nSeeded/nEv:0.);
  for(int p=0;p<3;++p)
    printf("  plane %d  p+pi merged AND mu separate     : %ld  (%.1f%%)\n",p,nPlaneOK[p],nEv?100.*nPlaneOK[p]/nEv:0.);
  printf("  PASS on at least one plane               : %ld  (%.1f%%)\n",nPassAny,nEv?100.*nPassAny/nEv:0.);
  fflush(stdout);
}
