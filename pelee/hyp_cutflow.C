// hyp_cutflow.C -- run the Sigma0 candidate logic on HyperonProduction ntuples
// through HypInput, and measure for the FIRST time how often each candidate is
// the particle it claims to be. The PeLEE truth record cannot answer this: it
// stops at the primary vertex and never contains the Lambda decay products.
#include "HypInput.h"
#include "PidAssign.h"
#include "SigmaFV.h"
#include "TVector3.h"
#include <fstream>

static const double MUON_TRK_SCORE=0.8, MUON_TRK_DIST=20.0, MUON_TRK_LEN=10.0, MUON_LLR=0.7;
static const int    NSHOWERS_MAX=4, NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5, TOPO_MIN=0.20;
static const double LLR_SCALE=0.591, DCHI_SCALE=124.5;

void hyp_cutflow(const char* list="lists/sig_genie.list",
                 Long64_t maxev=-1){
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l;
  while(std::getline(in,l)) if(l.size()) ch.Add(l.c_str());
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H;
  if(!H.attach(&ch)){ printf("  not a hypana tree\n"); return; }

  const char* CUT[8]={"AllEvents","FiducialV","CosmicVeto","NShowers",
                      "MuonCandidate","ProtonCandidate","PionCandidate","NTracksCut"};
  long nAll[8]={0}, nSig[8]={0};
  long muOK=0,muN=0, prOK=0,prN=0, piOK=0,piN=0, shOK=0,shN=0;
  std::map<int,long> piTruth, prTruth, shTruth;

  Long64_t N=ch.GetEntries(); if(maxev>0&&maxev<N) N=maxev;
  for(Long64_t ie=0;ie<N;++ie){
    ch.GetEntry(ie); H.convert();
    const bool isSig = (H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
    auto tick=[&](int c){ nAll[c]++; if(isSig) nSig[c]++; };
    tick(0);
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;  tick(1);
    if(topo<TOPO_MIN) continue;                                   tick(2);

    const int n=(int)H.trk_score.size();
    std::vector<TVector3> start(n); std::vector<double> dvtx(n,1e9);
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrack=[&](int i){ return H.trk_score[i]>TRK_SHR_SPLIT; };
    for(int i=0;i<n;++i){
      if(isTrack(i)) start[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else           start[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]);
      dvtx[i]=(start[i]-vtx).Mag();
    }
    int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrack(i)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;                      tick(3);

    int mu=-1;
    for(int i=0;i<n;++i)
      if(isTrack(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dvtx[i]<MUON_TRK_DIST&&
         H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
        if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
    if(mu<0) continue;                                            tick(4);

    // joint proton/pion assignment (see PidAssign.h); the old rule took the two
    // extremes of a single monotonic discriminant and mixed photons/muons into
    // the pion slot.
    std::vector<int> cand;
    for(int i=0;i<n;++i) if(i!=mu&&isTrack(i)) cand.push_back(i);
    if(cand.empty()) continue;                                    tick(5);
    std::vector<double> dvv(n), nhv(n);
    for(int i=0;i<n;++i){ dvv[i]=dvtx[i];
      nhv[i]=(i<(int)H.pfnplanehits_Y.size())?H.pfnplanehits_Y[i]:0; }
    int pr=-1, pi=-1;
    if(!PidAssign::assign(cand,H.trk_score,H.trk_llr,H.trk_chipr,H.trk_chipi,
                          H.trk_len,dvv,nhv,pr,pi)) continue;      tick(6);
    if(!(ntr<NTRACKS_MAX)) continue;                              tick(7);

    int s1=-1; double b1=1e9;
    for(int i=0;i<n;++i){ if(i==mu||i==pr||i==pi||isTrack(i))continue;
      if(dvtx[i]<b1){b1=dvtx[i];s1=i;} }

    if(isSig){   // candidate purity, signal events only
      muN++; if(H.bkt_pdg[mu]==-13) muOK++;
      prN++; if(H.bkt_pdg[pr]==2212) prOK++;  prTruth[H.bkt_pdg[pr]]++;
      piN++; if(H.bkt_pdg[pi]==-211) piOK++;  piTruth[H.bkt_pdg[pi]]++;
      if(s1>-1){ shN++; if(H.bkt_pdg[s1]==22) shOK++; shTruth[H.bkt_pdg[s1]]++; }
    }
  }
  printf("\n  === cut flow (%lld events, hypana list: %s) ===\n",N,list);
  printf("  %-16s %10s %10s %9s\n","cut","all","signal","eff[%]");
  for(int c=0;c<8;++c)
    printf("  %-16s %10ld %10ld %9.2f\n",CUT[c],nAll[c],nSig[c],
           nSig[0]?100.*nSig[c]/nSig[0]:0.);
  printf("\n  === candidate purity on SIGNAL events (impossible in PeLEE) ===\n");
  auto pc=[&](const char*n,long ok,long t){ printf("  %-22s %6ld / %-6ld = %5.1f%%\n",n,ok,t,t?100.*ok/t:0.); };
  pc("muon is a true mu+",muOK,muN);
  pc("proton is a true p",prOK,prN);
  pc("pion is a true pi-",piOK,piN);
  pc("shower1 is a true gamma",shOK,shN);
  auto dump=[&](const char*n,std::map<int,long>&m,long t){
    printf("\n  %s truth composition:\n",n);
    std::vector<std::pair<long,int>> v; for(auto&p:m) v.push_back({p.second,p.first});
    std::sort(v.rbegin(),v.rend());
    for(size_t k=0;k<v.size()&&k<7;++k) printf("      pdg %6d : %5.1f%%\n",v[k].second,t?100.*v[k].first/t:0.); };
  dump("pion candidate",piTruth,piN);
  dump("proton candidate",prTruth,prN);
  dump("shower1 candidate",shTruth,shN);
}
