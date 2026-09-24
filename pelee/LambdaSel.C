// ============================================================================
//  LambdaSel.C -- the Lambda CCQE selection (Thorpe & Nowak,
//  Hyperon_Technote_v2) reimplemented as closely as our ntuples allow, and run
//  on OUR samples.
//
//  PURPOSE IS CROSS-CHECK, NOT PHYSICS. This selects Lambdas, not Sigma0s. If it
//  reproduces their published yields on our files, then our sample
//  normalisation, POT bookkeeping, reconstruction chain and truth handling all
//  validate against an independent, published result. If it does not, the
//  disagreement localises the problem -- and every Sigma0 number in this
//  repository inherits whatever is wrong.
//
//  THEIR SIGNAL DEFINITION (technote sec 1.3), implemented exactly:
//    1. direct hyperon production, the hyperon leaving the nucleus being a Lambda,
//       anti-numu CC:  nubar_mu + Ar -> mu+ + Lambda + X
//    2. Lambda -> p + pi-                                    (BF 63.9%)
//    3. true p momentum > 300 MeV/c and true pi- momentum > 100 MeV/c
//    4. true vertex in the fiducial volume
//
//  THEIR PUBLISHED RESULT, RHC GENIE, normalised to 1.3e21 POT:
//    4.9 +- 0.6 signal, 0.5 +- 0.2 background, efficiency 6.8 +- 0.8%,
//    purity 91 (+4/-5) %.
//  Scaled to our data POT of 1.1526e21 that is 4.34 signal.
//
//  DELIBERATE DIFFERENCES, each forced by what we have:
//    * no photon anywhere -- this is the Lambda channel;
//    * the pair selector is the provisional score in LambdaStyle.h, not their
//      trained pair BDT, so pair-choice efficiency is a lower bound on theirs;
//    * their decay-analysis BDT is replaced by explicit cuts on the same
//      quantities (mass, pointing angle, decay length);
//    * we add the flash-matched slice requirement, which their PeLEE-independent
//      chain gets implicitly.
//  All of these can only make our efficiency WORSE than theirs, never better.
//  So efficiency below theirs is expected; efficiency ABOVE theirs, or a wildly
//  different background, would indicate a real problem.
//
//  Run:  root -l -b -q 'LambdaSel.C("lists/sig_genie.list",1500,true)'
//        root -l -b -q 'LambdaSel.C("lists/bkg_overlay.list",400,false)'
// ============================================================================
#include "HypInput.h"
#include "SigmaFV.h"
#include "Island.h"
#include "LambdaStyle.h"
#include "TChain.h"
#include "TVector3.h"
#include <fstream>
#include <vector>
#include <cstdio>

using namespace LambdaStyle;

static const double MUON_TRK_SCORE = 0.8;
static const double MUON_TRK_DIST  = 20.0;
static const double MUON_TRK_LEN   = 10.0;
static const double MUON_LLR       = 0.7;
static const int    NTRACKS_MAX    = 5;
static const double TRK_SHR_SPLIT  = 0.5;
static const double TOPO_MIN       = 0.20;
// Their signal thresholds
static const double P_MOM_MIN      = 0.300;   // GeV/c
static const double PI_MOM_MIN     = 0.100;   // GeV/c
// Decay analysis, standing in for their BDT
static const double LAM_MASS_WIN   = 0.040;
static const double COS_ALPHA_MIN  = 0.90;
static const double DECAY_LEN_MIN  = 0.5;

static const int NC = 9;
static const char* CUTNAME[NC] = {
  "AllEvents", "SliceID", "FiducialV", "CosmicVeto", "NTracks",
  "MuonCandidate", "PairCandidate", "LambdaMass", "PointingAlpha" };

void LambdaSel(const char* list = "lists/sig_genie.list",
               int maxFiles = 1500, bool signalOnly = true,
               bool applyIsland = true)
{
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf = 0;
  while(std::getline(in,l)) if(l.size() && (maxFiles<=0 || nf<maxFiles)){ ch.Add(l.c_str()); nf++; }

  double topo = 0; ch.SetBranchAddress("topological_score", &topo);
  HypInput H; if(!H.attach(&ch)){ printf("  not a hypana tree\n"); return; }
  // bind extras AFTER attach()
  int fmSlice = -1;
  const bool haveSlice = ch.GetBranch("flash_match_nu_slice_ID") != nullptr;
  if(haveSlice) ch.SetBranchAddress("flash_match_nu_slice_ID", &fmSlice);
  std::vector<double>* d3 = nullptr;
  const bool haveDedx = ch.GetBranch("trk_three_plane_dedx") != nullptr;
  if(haveDedx) ch.SetBranchAddress("trk_three_plane_dedx", &d3);
  std::vector<int>*    dpdg = nullptr;
  std::vector<double>* dmom = nullptr;
  if(ch.GetBranch("mc_decay_pdg")) ch.SetBranchAddress("mc_decay_pdg", &dpdg);
  if(ch.GetBranch("mc_decay_mom")) ch.SetBranchAddress("mc_decay_mom", &dmom);
  std::vector<std::vector<float>>* win[3] = {0,0,0};
  std::vector<int>* vw = nullptr;
  for(int p=0;p<3;++p) ch.SetBranchAddress(Form("ct_test_window_plane%d",p), &win[p]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires", &vw);

  long nCut[NC] = {0};
  long nIsland = 0;
  const PairCuts pc;

  const Long64_t N = ch.GetEntries();
  printf("\n  %s : %d files, %lld events\n", signalOnly?"LAMBDA SIGNAL":"BACKGROUND", nf, N);
  fflush(stdout);

  for(Long64_t ie = 0; ie < N; ++ie){
    ch.GetEntry(ie); H.convert();

    // ── their signal definition ──────────────────────────────────────────────
    bool isSig = (H.h_hyperon_pdg==3122 && H.nu_pdg==-14 && H.ccnc==0);
    if(isSig){
      bool hasP=false, hasPi=false; double pP=0, pPi=0;
      if(dpdg && dmom) for(size_t j=0;j<dpdg->size() && j<dmom->size(); ++j){
        if(dpdg->at(j)== 2212){ hasP =true; pP =dmom->at(j); }
        if(dpdg->at(j)==-211 ){ hasPi=true; pPi=dmom->at(j); } }
      if(!(hasP && hasPi))                       isSig = false;   // Lambda -> p pi-
      else if(!(pP > P_MOM_MIN && pPi > PI_MOM_MIN)) isSig = false;   // thresholds
      // criterion 4 (true vertex in FV) is applied below, with the reco FV cut,
      // since the denominator must count only in-FV signal
      else if(!SigmaFV::isFiducial(H.true_vtx_x,H.true_vtx_y,H.true_vtx_z)) isSig = false;
    }
    if(signalOnly != isSig) continue;
    nCut[0]++;

    if(haveSlice && fmSlice < 0) continue;                              nCut[1]++;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;         nCut[2]++;
    if(topo < TOPO_MIN) continue;                                       nCut[3]++;

    const int n = (int)H.trk_score.size(); if(n < 1) continue;
    const TVector3 vtx(H.vtx_x, H.vtx_y, H.vtx_z);
    auto isTrk = [&](int i){ return H.trk_score[i] > TRK_SHR_SPLIT; };
    std::vector<TVector3> st(n), dir(n);
    std::vector<double> dv(n, 1e9);
    for(int i=0;i<n;++i){
      if(isTrk(i)) st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else         st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]);
      dir[i].SetXYZ(H.trk_dx[i],H.trk_dy[i],H.trk_dz[i]);
      dv[i] = (st[i]-vtx).Mag();
    }
    int ntr=0; for(int i=0;i<n;++i) if(isTrk(i)) ntr++;
    // their topology is muon + V: at least three tracks, and no shower requirement
    if(!(ntr >= 3 && ntr < NTRACKS_MAX)) continue;                      nCut[4]++;

    int mu = -1;
    for(int i=0;i<n;++i)
      if(isTrk(i) && H.trk_score[i]>MUON_TRK_SCORE && dv[i]<MUON_TRK_DIST &&
         H.trk_len[i]>MUON_TRK_LEN && H.trk_llr[i]>MUON_LLR)
        if(mu==-1 || H.trk_len[i]>H.trk_len[mu]) mu = i;
    if(mu < 0) continue;                                                nCut[5]++;

    int pr = -1, pi = -1; double bestScore = -1e30;
    for(int a=0;a<n;++a){
      if(a==mu || !isTrk(a)) continue;
      for(int b=0;b<n;++b){
        if(b==mu || b==a || !isTrk(b)) continue;
        PairFeat f;
        f.sep    = (st[a]-st[b]).Mag();
        f.scoreP = H.trk_score[a];  f.scorePi = H.trk_score[b];
        f.llrP   = H.trk_llr[a];    f.llrPi   = H.trk_llr[b];
        if(haveDedx && d3 && a<(int)d3->size() && std::isfinite(d3->at(a))) f.dedxP  = d3->at(a);
        if(haveDedx && d3 && b<(int)d3->size() && std::isfinite(d3->at(b))) f.dedxPi = d3->at(b);
        if(!passPrecuts(f, pc)) continue;
        const double s = provisionalPairScore(f);
        if(s > bestScore){ bestScore = s; pr = a; pi = b; }
      }
    }
    if(pr < 0 || pi < 0) continue;                                      nCut[6]++;

    const double pMom = protonMomFromRange(H.trk_len[pr]);
    const double iMom = H.trk_rmumom[pi];
    const DecayObs obs = decayObservables(vtx,
                            st[pr], dir[pr], H.trk_len[pr],
                            st[pi], dir[pi], H.trk_len[pi],
                            pMom, iMom);
    if(!obs.ok) continue;
    if(!(obs.mass > 0) || std::fabs(obs.mass - M_LAMBDA) > LAM_MASS_WIN) continue;  nCut[7]++;
    if(!(obs.decayLength > DECAY_LEN_MIN)) continue;
    if(!(obs.cosAlpha > COS_ALPHA_MIN)) continue;                       nCut[8]++;

    if(applyIsland){
      bool ok2 = false;
      if(vw && vw->size() > 2){
        for(int p=0;p<3;++p){
          if(!win[p] || win[p]->empty() || vw->at(p) < 0) continue;
          bool ok = false;
          if(Island::planeTest(*win[p], vw->at(p)-(int)win[p]->size()/2, p,
                               st[mu].X(),st[mu].Y(),st[mu].Z(),
                               st[pr].X(),st[pr].Y(),st[pr].Z(),
                               st[pi].X(),st[pi].Y(),st[pi].Z(), ok)){ ok2 = true; break; }
        }
      }
      if(ok2) nIsland++;
    }
  }

  printf("\n  %-16s %10s %9s\n", "cut", "events", "rel [%]");
  for(int c=0;c<NC;++c)
    printf("  %-16s %10ld %8.1f\n", CUTNAME[c], nCut[c],
           (c>0 && nCut[c-1]) ? 100.0*nCut[c]/nCut[c-1] : 100.0);
  if(applyIsland) printf("  %-16s %10ld %8.1f\n", "Island", nIsland,
           nCut[NC-1] ? 100.0*nIsland/nCut[NC-1] : 0.0);
  printf("\n  absolute efficiency from the signal definition : %.2f%%   (theirs 6.8%%)\n",
         nCut[0] ? 100.0*(applyIsland?nIsland:nCut[NC-1])/nCut[0] : 0.0);
  fflush(stdout);
}
