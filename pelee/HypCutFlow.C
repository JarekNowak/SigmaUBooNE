// ============================================================================
//  HypCutFlow.C -- POT-weighted cut flow for the hypana (HyperonProduction)
//  chain, in the same shape as the PeLEE table that SigmaSelection.C prints.
//
//  Runs BOTH samples in one pass and scales each by its own generated POT:
//    signal      lists/sig_genie.list   true Sigma0 CC anti-numu only
//    background  lists/bkg_overlay.list everything that is not that
//  matching exactly how DumpBDTHyp.C splits the two, so the numbers line up
//  with the BDT training tree.
//
//  POT is summed per sample from ana/MetaTree, so the scales are computed here
//  rather than hard-coded and stay correct if the file lists change.
//
//  Also reports candidate purity on signal events at the final stage, which is
//  the thing only hypana truth can answer.
//
//  Run:  root -l -b -q HypCutFlow.C
// ============================================================================
#include "HypInput.h"
#include "PidAssign.h"
#include "SigmaFV.h"
#include "TChain.h"
#include "TFile.h"
#include "TVector3.h"
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cmath>

static const double DATAPOT=1.1526e21;
static const double MUON_TRK_SCORE=0.8, MUON_TRK_DIST=20.0, MUON_TRK_LEN=10.0, MUON_LLR=0.7;
static const int    NSHOWERS_MAX=4, NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5, TOPO_MIN=0.20;

//  Two chains. The default drops the proton and pion REQUIREMENTS: they are not
//  PID cuts but multiplicity demands (>=1 and >=2 non-muon tracks), and measured
//  against the POT-weighted flow they remove more signal than background
//  (-32.7%/-19.2% and -52.8%/-40.5%), costing 54% of the significance. Candidates
//  are still formed when two non-muon tracks exist; when they are not, the event is
//  KEPT and the mass variables are simply unavailable -- a flag for the classifier
//  rather than a cut. Set requireProtonPion=true to restore the old chain.
static const int NC_OPT=7;
static const char* CUT_OPT[NC_OPT]={"AllEvents","SliceID","FiducialV","CosmicVeto",
                                    "NShowers","MuonCandidate","NTracksCut"};
static const int NC_REQ=9;
static const char* CUT_REQ[NC_REQ]={"AllEvents","SliceID","FiducialV","CosmicVeto","NShowers",
                                    "MuonCandidate","ProtonCandidate","PionCandidate","NTracksCut"};
static const int NC=9;   // array size; only the first NCuse rows are filled/printed

static std::vector<std::string> loadList(const char* p){
  std::ifstream in(p); std::string l; std::vector<std::string> v;
  while(std::getline(in,l)) if(l.size()) v.push_back(l);
  return v;
}

//  sigPOT/bkgPOT default to the values measured over these exact lists by summing
//  ana/MetaTree (as DumpBDTHyp.C does). Pass <=0 to re-scan the files instead --
//  correct but slow, since it opens every file just to read a one-entry tree.
//  bkgPOT was 2.7795e21 for the 3491-file Run 3 list; on 2026-09-24 the RHC overlay
//  extension (Runs 1/2/3/4a, 10,533 files, 4.1906e21) was appended to
//  lists/bkg_overlay.list, giving 14,024 files and 6.9701e21. The old list is kept
//  as lists/bkg_overlay_run3.list -- pass it with bkgPOT=2.7795e21 to reproduce
//  earlier numbers. Re-measure both after the residue pass adds its files.
void HypCutFlow(const char* sigList="lists/sig_genie.list",
                const char* bkgList="lists/bkg_overlay.list",
                double sigPOT=5.9770e22, double bkgPOT=6.9701e21,
                bool requireProtonPion=false)
{
  const int NCuse = requireProtonPion ? NC_REQ : NC_OPT;
  const char** CUT = requireProtonPion ? CUT_REQ : CUT_OPT;
  long nPair=0;   // events at the final stage that DO have a p/pi pair
  double raw[2][NC]={{0}}, wgt[2][NC]={{0}}, scale[2]={0,0};
  long muOK=0,prOK=0,piOK=0,shOK=0,nCand=0;
  std::map<int,long> piTruth;

  for(int pass=0; pass<2; ++pass){
    const bool wantSignal = (pass==0);
    auto files = loadList(wantSignal? sigList : bkgList);

    // ── generated POT: use the supplied value, else scan the files ──────────
    double pot = wantSignal ? sigPOT : bkgPOT;
    std::vector<std::string> good;
    if(pot>0){ good=files; }
    else {
      pot=0;
      for(auto& fn : files){
        TFile* f=TFile::Open(fn.c_str(),"READ");
        if(!f||f->IsZombie()){ if(f) delete f; continue; }
        if(!f->Get("ana/OutputTree")){ f->Close(); delete f; continue; }
        TTree* mt=(TTree*)f->Get("ana/MetaTree");
        if(mt&&mt->GetBranch("POT")){ double p=0; mt->SetBranchAddress("POT",&p); mt->GetEntry(0); pot+=p; }
        f->Close(); delete f; good.push_back(fn);
      }
    }
    scale[pass] = pot>0 ? DATAPOT/pot : 0.0;
    printf("  %-10s %zu files, POT %.4e, scale to data POT %.4e\n",
           wantSignal?"signal":"background", good.size(), pot, scale[pass]);
    fflush(stdout);

    TChain ch("ana/OutputTree");
    for(auto& fn : good) ch.Add(fn.c_str());
    double topo=0; ch.SetBranchAddress("topological_score",&topo);
    HypInput H; if(!H.attach(&ch)){ printf("  not a hypana tree\n"); return; }
    // Neutrino-slice requirement. PeLEE applies this implicitly: 58.3% of its
    // overlay events have nslice==0 / zero PFPs and no usable vertex, so they
    // die at the fiducial cut. hypana keeps them with a Pandora vertex attached,
    // which is why its FiducialV kept 73% against PeLEE's 30%. The matching
    // quantity is the FLASH-matched neutrino slice: 57.0% of hypana events have
    // none, against PeLEE's 58.3%. Bound after attach() -- see HypInput.h.
    int fmSliceID=-1;
    const bool hasSliceBranch = ch.GetBranch("flash_match_nu_slice_ID")!=nullptr;
    if(hasSliceBranch) ch.SetBranchAddress("flash_match_nu_slice_ID",&fmSliceID);

    const Long64_t N=ch.GetEntries();
    for(Long64_t ie=0; ie<N; ++ie){
      ch.GetEntry(ie); H.convert();
      const bool isSig=(H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
      if(wantSignal != isSig) continue;          // same split as DumpBDTHyp.C
      auto tick=[&](int c){ raw[pass][c]+=1.0; wgt[pass][c]+=scale[pass]; };

      tick(0);
      if(hasSliceBranch && fmSliceID<0) continue;                   tick(1);
      if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;  tick(2);
      if(topo<TOPO_MIN) continue;                                   tick(3);

      const int n=(int)H.trk_score.size();
      if(n<1) continue;
      TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
      auto isTrk=[&](int i){ return H.trk_score[i]>TRK_SHR_SPLIT; };
      std::vector<TVector3> st(n); std::vector<double> dv(n,1e9);
      for(int i=0;i<n;++i){
        if(isTrk(i)) st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
        else         st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]);
        dv[i]=(st[i]-vtx).Mag();
      }
      int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
      if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;                      tick(4);

      int mu=-1;
      for(int i=0;i<n;++i)
        if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dv[i]<MUON_TRK_DIST&&
           H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
          if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
      if(mu<0) continue;                                            tick(5);

      std::vector<int> cand;
      for(int i=0;i<n;++i) if(i!=mu&&isTrk(i)) cand.push_back(i);
      std::vector<double> dvv(n), nhv(n);
      for(int i=0;i<n;++i){ dvv[i]=dv[i];
        nhv[i]=(i<(int)H.pfnplanehits_Y.size())?H.pfnplanehits_Y[i]:0; }
      int pr=-1, pi=-1;
      const bool havePair = PidAssign::assign(cand,H.trk_score,H.trk_llr,H.trk_chipr,
                                              H.trk_chipi,H.trk_len,dvv,nhv,pr,pi);
      if(requireProtonPion){
        if(cand.empty()) continue;                                  tick(6);
        if(!havePair) continue;                                     tick(7);
        if(!(ntr<NTRACKS_MAX)) continue;                            tick(8);
      } else {
        if(!(ntr<NTRACKS_MAX)) continue;                            tick(6);
      }
      if(wantSignal && havePair) nPair++;

      if(wantSignal && havePair){           // candidate purity, where a pair exists
        int s1=-1; double b1=1e9;
        for(int i=0;i<n;++i){ if(i==mu||i==pr||i==pi||isTrk(i)) continue;
          if(dv[i]<b1){ b1=dv[i]; s1=i; } }
        nCand++;
        if(H.bkt_pdg[mu]==-13)  muOK++;
        if(H.bkt_pdg[pr]==2212) prOK++;
        if(H.bkt_pdg[pi]==-211) piOK++;
        piTruth[H.bkt_pdg[pi]]++;
        if(s1>-1 && H.bkt_pdg[s1]==22) shOK++;
      }
    }
  }

  printf("\n  =============== POT-weighted cut flow (scaled to %.4e POT) ===============\n",DATAPOT);
  printf("  %-16s %9s %9s %11s %9s %9s %11s\n",
         "Cut","RawSig","RawBkg","S","B","Eff[%]","Pur[%]");
  for(int c=0;c<NCuse;++c){
    const double S=wgt[0][c], B=wgt[1][c];
    printf("  %-16s %9.0f %9.0f %11.3f %9.1f %9.2f %11.4f\n",
           CUT[c], raw[0][c], raw[1][c], S, B,
           wgt[0][0]>0?100.*S/wgt[0][0]:0.,
           (S+B)>0?100.*S/(S+B):0.);
  }
  printf("\n  %-16s %9s\n","Cut","S/sqrt(S+B)");
  for(int c=0;c<NCuse;++c){
    const double S=wgt[0][c], B=wgt[1][c];
    printf("  %-16s %9.4f\n",CUT[c],(S+B)>0?S/std::sqrt(S+B):0.);
  }
  printf("\n  signal events at the final stage : %.0f\n",raw[0][NCuse-1]);
  printf("  of those, a p/pi pair was formed : %ld  (%.1f%%)\n",
         nPair, raw[0][NCuse-1]>0?100.*nPair/raw[0][NCuse-1]:0.);
  printf("\n  === candidate purity where a pair exists (%ld events) ===\n",nCand);
  auto pc=[&](const char* nm,long ok){ printf("  %-24s %5ld / %-5ld = %5.1f%%\n",
                nm,ok,nCand,nCand?100.*ok/nCand:0.); };
  pc("muon is a true mu+",muOK);
  pc("proton is a true p",prOK);
  pc("pion is a true pi-",piOK);
  pc("shower1 is a true gamma",shOK);
  printf("\n  pion candidate truth composition:\n");
  std::vector<std::pair<long,int>> v; for(auto&p:piTruth) v.push_back({p.second,p.first});
  std::sort(v.rbegin(),v.rend());
  for(size_t k=0;k<v.size()&&k<7;++k)
    printf("      pdg %6d : %5.1f%%\n",v[k].second,nCand?100.*v[k].first/nCand:0.);
  fflush(stdout);
}
