// ============================================================================
//  SigmaLambdaSel.C -- an ALTERNATIVE Sigma0 selection that follows the
//  Lambda CCQE chain (Thorpe & Nowak, Hyperon_Technote_v2) as closely as the
//  channel allows, instead of the topology+PID chain in SigmaSelection.C.
//
//  MOTIVATION (all measured, see the memory notes):
//    At equal efficiency (6.7% vs 6.8%) their purity is 91% and ours 2.8%.
//    Stacking their cuts onto our existing selection tops out near 1.9%, so the
//    gap is structural rather than a working point. The three causes:
//      1. our preselection destroys the displaced vertex, their main purifier;
//      2. we have no pair precuts, so no event is ever rejected for lacking a
//         decay-like pair -- we made the pair optional entirely;
//      3. our signal needs an extra object, the photon.
//
//  WHAT IS DIFFERENT HERE, relative to SigmaSelection.C:
//    * The proton/pion pair is REQUIRED, and chosen by scoring every ordered
//      pair rather than ranking single tracks. Their selector's topological
//      input is the separation of the two track STARTS. Ours used the distance
//      to the PRIMARY vertex, with a negative coefficient on the pion, which
//      penalises exactly the displaced decays the method depends on.
//    * Precuts can reject an event outright when no pair is decay-like.
//    * The decay vertex is fitted (a V fit), giving a decay length and the
//      pointing angle alpha. Our previous displacement variables were built from
//      a closest-approach midpoint and measured uncorrelated with truth.
//    * The connectedness island test enters as a cut, not a BDT input, and is
//      applied where it is strongest.
//    * The photon is an ADDITION at the end (Sigma0 -> Lambda gamma), not the
//      organising principle.
//
//  OPERATING POINT (2026-09-12, full-statistics scan in IslandAna.C on the cache
//  written by IslandProcess.C): no pair precuts, |m - m_Lambda| < 20 MeV,
//  cos alpha > 0.95, decay length > 0.5 cm, island on >= 2 planes with the tuned
//  connectivity in Island.h (2 wires x 20 ticks). Cross-validated by event parity
//  (IslandXVal.C, IslandTune.C::IslandPoint). Expected at data POT: S = 1.39,
//  B = 4.1 (10 raw MC events; parity halves 4 and 6), purity 25%. The one-plane
//  point (islandPlanes=1) is S = 1.99, B = 19.5 and is the statistically safer
//  background estimate.
//
//  This macro measures the chain; it does not yet write histograms or a BDT
//  tree. Run it on signal and background separately and compare the tables.
//
//  Run:  root -l -b -q 'SigmaLambdaSel.C("lists/sig_genie.list",700,true)'
//        root -l -b -q 'SigmaLambdaSel.C("lists/bkg_overlay.list",250,false)'
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

// ── Muon identification: their requirements, unchanged from our chain ────────
static const double MUON_TRK_SCORE = 0.8;
static const double MUON_TRK_DIST  = 20.0;   // cm from the primary vertex
static const double MUON_TRK_LEN   = 10.0;   // cm
static const double MUON_LLR       = 0.7;
// ── Preselection ────────────────────────────────────────────────────────────
static const int    NSHOWERS_MAX   = 4;
static const int    NTRACKS_MAX    = 5;
static const double TRK_SHR_SPLIT  = 0.5;
static const double TOPO_MIN       = 0.20;
// ── Decay analysis, standing in for their decay BDT ─────────────────────────
static const double LAM_MASS_WIN   = 0.020;  // GeV about the Lambda mass (scan 2026-09-12, was 0.040)
static const double COS_ALPHA_MIN  = 0.95;   // Lambda must point at its decay vertex (scan 2026-09-12, was 0.90)
static const double DECAY_LEN_MIN  = 0.5;    // cm, a decay that actually flew

static const int NC = 10;
static const char* CUTNAME[NC] = {
  "AllEvents", "SliceID", "FiducialV", "CosmicVeto", "NShowers",
  "MuonCandidate", "PairCandidate", "LambdaMass", "PointingAlpha", "Island" };

void SigmaLambdaSel(const char* list = "lists/sig_genie.list",
                    int maxFiles = 700, bool signalOnly = true,
                    bool requirePhoton = true, bool applyIsland = true,
                    bool usePrecuts = false, int islandPlanes = 2,
                    bool applyMichel = false, bool photonGap = false)
{
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf = 0;
  while(std::getline(in,l)) if(l.size() && (maxFiles<=0 || nf<maxFiles)){ ch.Add(l.c_str()); nf++; }

  double topo = 0; ch.SetBranchAddress("topological_score", &topo);
  HypInput H; if(!H.attach(&ch)){ printf("  not a hypana tree\n"); return; }
  // NOTE: extra branches must be bound AFTER attach() -- binding first leaves
  // them empty on a TChain.
  int fmSlice = -1;
  const bool haveSlice = ch.GetBranch("flash_match_nu_slice_ID") != nullptr;
  if(haveSlice) ch.SetBranchAddress("flash_match_nu_slice_ID", &fmSlice);
  std::vector<double>* d3 = nullptr;
  const bool haveDedx = ch.GetBranch("trk_three_plane_dedx") != nullptr;
  if(haveDedx) ch.SetBranchAddress("trk_three_plane_dedx", &d3);
  std::vector<std::vector<float>>* win[3] = {0,0,0};
  std::vector<int>* vw = nullptr;
  for(int p=0;p<3;++p) ch.SetBranchAddress(Form("ct_test_window_plane%d",p), &win[p]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires", &vw);
  // truth Lambda flight, for diagnostics only (filled for signal)
  double lsx=0,lsy=0,lsz=0,lex=0,ley=0,lez=0;
  auto B=[&](const char* n, void* q){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,q); };
  std::vector<double> *tex=nullptr,*tey=nullptr,*tez=nullptr;   // muon end, Michel tag
  B("trk_end_x",&tex); B("trk_end_y",&tey); B("trk_end_z",&tez);
  B("mc_sigmazero_lambda_start_x",&lsx); B("mc_sigmazero_lambda_start_y",&lsy);
  B("mc_sigmazero_lambda_start_z",&lsz); B("mc_sigmazero_lambda_end_x",&lex);
  B("mc_sigmazero_lambda_end_y",&ley);   B("mc_sigmazero_lambda_end_z",&lez);

  long nCut[NC] = {0};
  long nPairTried = 0, nPairPassed = 0;
  // diagnostics on the surviving sample
  long nLongFlight = 0, nSurv = 0;
  const PairCuts pc;

  const Long64_t N = ch.GetEntries();
  printf("\n  %s : %d files, %lld events\n", signalOnly?"SIGNAL":"BACKGROUND", nf, N);
  printf("  requirePhoton=%d  applyIsland=%d  usePrecuts=%d  islandPlanes=%d  massWin=%.0f MeV  cosAlpha>%.2f\n",
         (int)requirePhoton, (int)applyIsland, (int)usePrecuts, islandPlanes, 1000*LAM_MASS_WIN, COS_ALPHA_MIN);
  fflush(stdout);

  for(Long64_t ie = 0; ie < N; ++ie){
    ch.GetEntry(ie); H.convert();
    const bool isSig = (H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
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
    int ntr=0, nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
    if(!(nsh>0 && nsh<NSHOWERS_MAX)) continue;
    if(!(ntr < NTRACKS_MAX)) continue;                                  nCut[4]++;

    int mu = -1;
    for(int i=0;i<n;++i)
      if(isTrk(i) && H.trk_score[i]>MUON_TRK_SCORE && dv[i]<MUON_TRK_DIST &&
         H.trk_len[i]>MUON_TRK_LEN && H.trk_llr[i]>MUON_LLR)
        if(mu==-1 || H.trk_len[i]>H.trk_len[mu]) mu = i;
    if(mu < 0) continue;                                                nCut[5]++;

    // ── Pair selection: score EVERY ordered pair, keep the best that passes ──
    // Precuts measured as dead weight as a VETO on the full samples (IslandAna.C,
    // 2026-09-12): same FoM without them, ~20% more signal. They are kept as a
    // PREFERENCE: the best precut-passing pair is taken when one exists, else the
    // best pair of all. This is the logic IslandProcess.C caches, so the numbers
    // here reproduce IslandAna.C. usePrecuts=true restores the veto.
    int pr = -1, pi = -1, prAny = -1, piAny = -1; double bestScore = -1e30, bestAny = -1e30;
    for(int a=0;a<n;++a){
      if(a==mu || !isTrk(a)) continue;
      for(int b=0;b<n;++b){
        if(b==mu || b==a || !isTrk(b)) continue;
        PairFeat f;
        f.sep     = (st[a]-st[b]).Mag();
        f.scoreP  = H.trk_score[a];  f.scorePi = H.trk_score[b];
        f.llrP    = H.trk_llr[a];    f.llrPi   = H.trk_llr[b];
        if(haveDedx && d3 && a<(int)d3->size() && std::isfinite(d3->at(a))) f.dedxP  = d3->at(a);
        if(haveDedx && d3 && b<(int)d3->size() && std::isfinite(d3->at(b))) f.dedxPi = d3->at(b);
        nPairTried++;
        const double s = provisionalPairScore(f);
        if(s > bestAny){ bestAny = s; prAny = a; piAny = b; }
        if(!passPrecuts(f, pc)) continue;
        nPairPassed++;
        if(s > bestScore){ bestScore = s; pr = a; pi = b; }
      }
    }
    if(pr < 0 && !usePrecuts){ pr = prAny; pi = piAny; }
    if(pr < 0 || pi < 0) continue;                                      nCut[6]++;

    // ── Decay analysis on the fitted V ──────────────────────────────────────
    const double pMom = protonMomFromRange(H.trk_len[pr]);
    const double iMom = H.trk_rmumom[pi];        // muon-hypothesis, deliberately
    const DecayObs obs = decayObservables(vtx,
                            st[pr], dir[pr], H.trk_len[pr],
                            st[pi], dir[pi], H.trk_len[pi],
                            pMom, iMom);
    if(!obs.ok) continue;
    if(!(obs.mass > 0) || std::fabs(obs.mass - M_LAMBDA) > LAM_MASS_WIN) continue;  nCut[7]++;
    if(!(obs.decayLength > DECAY_LEN_MIN)) continue;
    if(!(obs.cosAlpha > COS_ALPHA_MIN)) continue;                       nCut[8]++;

    // ── Photon: the Sigma0 addition ─────────────────────────────────────────
    int s1 = -1; double bs = 1e9;
    for(int i=0;i<n;++i){ if(i==mu||i==pr||i==pi||isTrk(i)) continue;
      if(dv[i] < bs){ bs = dv[i]; s1 = i; } }
    if(requirePhoton && s1 < 0) continue;
    // Optional handles measured in PhotonMichel.C (2026-09-12), off by default:
    //  * photonGap: the photon must convert away from the vertex (>2 cm). Keeps
    //    95% of signal, removes 21% of background at the one-plane point.
    //  * applyMichel: a contained muon must have another PFP starting within 5 cm
    //    of its end. What it removes is 87% true mu- (wrong-sign numu CC); the
    //    signal muon is a mu+ and always decays. +15% FoM on the one-plane point,
    //    nothing measurable on the two-plane point (8 MC background events).
    if(photonGap && s1 >= 0 && dv[s1] < 2.0) continue;
    if(applyMichel && tex && mu < (int)tex->size()){
      const TVector3 mEnd(tex->at(mu), tey->at(mu), tez->at(mu));
      if(SigmaFV::isFiducial(mEnd.X(), mEnd.Y(), mEnd.Z())){
        bool tagged=false;
        for(int i=0;i<n;++i){ if(i==mu||i==pr||i==pi) continue; if((st[i]-mEnd).Mag() < 5.0){ tagged=true; break; } }
        if(!tagged) continue;
      }
    }

    // ── Connectedness ───────────────────────────────────────────────────────
    // The test is evaluated on every plane and `islandPlanes` of them must pass.
    // Two planes is 6.8x rejection against 3.7x for one (IslandAna.C).
    bool islandOK = true;
    if(applyIsland){
      int nPass = 0;
      if(vw && vw->size() > 2){
        for(int p=0;p<3;++p){
          if(!win[p] || win[p]->empty() || vw->at(p) < 0) continue;
          bool ok = false;
          if(Island::planeTest(*win[p], vw->at(p)-(int)win[p]->size()/2, p,
                               st[mu].X(),st[mu].Y(),st[mu].Z(),
                               st[pr].X(),st[pr].Y(),st[pr].Z(),
                               st[pi].X(),st[pi].Y(),st[pi].Z(), ok)) nPass++;
        }
      }
      islandOK = (nPass >= islandPlanes);
    }
    if(!islandOK) continue;                                             nCut[9]++;

    nSurv++;
    const double fl = std::sqrt((lex-lsx)*(lex-lsx)+(ley-lsy)*(ley-lsy)+(lez-lsz)*(lez-lsz));
    if(fl > 2.0) nLongFlight++;
  }

  printf("\n  %-16s %10s %10s\n", "cut", "events", "rel [%]");
  for(int c=0;c<NC;++c)
    printf("  %-16s %10ld %9.1f\n", CUTNAME[c], nCut[c],
           (c>0 && nCut[c-1]) ? 100.0*nCut[c]/nCut[c-1] : 100.0);
  printf("\n  absolute efficiency from AllEvents : %.2f%%\n",
         nCut[0] ? 100.0*nCut[NC-1]/nCut[0] : 0.0);
  printf("  ordered pairs tried %ld, passing precuts %ld (%.1f%%)\n",
         nPairTried, nPairPassed, nPairTried ? 100.0*nPairPassed/nPairTried : 0.0);
  if(signalOnly && nSurv)
    printf("  surviving signal with true Lambda flight > 2 cm : %.1f%%  (was 43%% in the old chain)\n",
           100.0*nLongFlight/nSurv);
  fflush(stdout);
}
