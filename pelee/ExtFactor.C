// ============================================================================
//  ExtFactor.C -- factorised EXT (beam-off cosmic) background at the Lambda-style
//  operating point.
//
//  The hypana EXT carries 150 data events per event (ext_hypana README), so it
//  cannot be followed to the operating point itself: none survive, and zero means
//  < 346 events at 90% CL against an MC background of ~4. Instead
//    EXT(op) = EXT(MuonCandidate) x P(op | MuonCandidate, cosmic)
//  with the first factor from the hypana EXT cache and the second from the
//  overlay MC's own cosmic-dominated candidates (CosmicTag.C: cosFrac > 0.5, the
//  selected slice's hits mostly unmatched to MC = overlaid cosmic data). The two
//  samples are compared step by step wherever the EXT still has events, which is
//  the test of the factorisation.
//
//  Operating point as IslandOpPoint.C: best pair, |m - m_Lambda| < 20 MeV, decay
//  length > 0.5 cm, cos(alpha) > 0.95, island on >= 2 planes, photon gap > 2 cm.
//
//  Run:  root -l -b -q 'ExtFactor.C("isl_ext.root","isl_bkg.root","isl_bkg_cos.root")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

namespace {
static const int NSTEP = 6;
static const char* STEP[NSTEP] = { "MuonCandidate", "pair", "Lambda mass +-20 MeV",
                                   "decay>0.5, cos(alpha)>0.95", "island >= 2 planes", "photon gap > 2 cm" };
struct Row { double w[NSTEP] = {0}, w2[NSTEP] = {0}; long n[NSTEP] = {0}; };

// fill the cumulative step counts of one cache; cosCut<0: no cosmic selection
void fill(const char* fn, const char* friendFn, double cosCut, Row& r){
  TFile f(fn); TTree* t = (TTree*)f.Get("isl");
  if (!t) { printf("  no isl tree in %s\n", fn); return; }
  float cosFrac = -1;
  if (friendFn && friendFn[0]) { t->AddFriend("cos", friendFn); t->SetBranchAddress("cosFrac", &cosFrac); }
  int hasPair=0, islN=0; double w=0, mass=-1, len=-1, cosA=-2, phDist=-1;
  t->SetBranchAddress("hasPair",&hasPair); t->SetBranchAddress("islNPlanes",&islN);
  t->SetBranchAddress("weight",&w);        t->SetBranchAddress("mass",&mass);
  t->SetBranchAddress("decayLen",&len);    t->SetBranchAddress("cosAlpha",&cosA);
  t->SetBranchAddress("phDist",&phDist);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (cosCut >= 0 && !(cosFrac > cosCut)) continue;
    bool pass[NSTEP] = { true, (bool)hasPair, std::fabs(mass - 1.115683) < 0.020,
                         len > 0.5 && cosA > 0.95, islN >= 2, phDist > 2.0 };
    for (int k = 0; k < NSTEP; ++k) {
      bool ok = true; for (int j = 0; j <= k; ++j) ok = ok && pass[j];
      if (!ok) break;
      r.w[k] += w; r.w2[k] += w * w; r.n[k]++;
    }
  }
}
}

void ExtFactor(const char* extFile = "isl_ext.root", const char* bkgFile = "isl_bkg.root",
               const char* cosFile = "isl_bkg_cos.root", double cosCut = 0.5)
{
  Row ext, mcCos, mcAll;
  fill(extFile, "", -1, ext);
  fill(bkgFile, cosFile, cosCut, mcCos);
  fill(bkgFile, "", -1, mcAll);

  printf("\n  %-28s | %-24s | %-24s | %-14s\n", "step (cumulative)", "EXT (data-scaled)",
         "MC cosmic-dominated", "MC all bkg");
  printf("  %-28s | %10s %13s | %10s %13s | %14s\n", "", "raw", "P(step|prev)", "raw", "P(step|prev)", "weighted");
  for (int k = 0; k < NSTEP; ++k) {
    auto P = [&](const Row& r){ return (k && r.n[k-1]) ? (double)r.n[k] / r.n[k-1] : 1.0; };
    auto E = [&](const Row& r){ const double p = P(r); return (k && r.n[k-1]) ? std::sqrt(p*(1-p)/r.n[k-1]) : 0.0; };
    printf("  %-28s | %10ld %6.3f+-%5.3f | %10ld %6.3f+-%5.3f | %14.2f\n", STEP[k],
           ext.n[k], P(ext), E(ext), mcCos.n[k], P(mcCos), E(mcCos), mcAll.w[k]);
  }

  // factorised estimate: EXT at MuonCandidate x MC-cosmic pass fraction to the end
  const double extMu = ext.w[0], extMuErr = std::sqrt(ext.w2[0]);
  const long nMu = mcCos.n[0], nOp = mcCos.n[NSTEP-1];
  const double p = nMu ? (double)nOp / nMu : 0;
  printf("\n  EXT at MuonCandidate        %.0f +- %.0f  (%ld raw)\n", extMu, extMuErr, ext.n[0]);
  printf("  P(op | MuonCand), MC cosmic  %ld / %ld = %.2e\n", nOp, nMu, p);
  if (nOp > 0) {
    const double rel = std::sqrt(1.0/nOp + (ext.n[0] ? 1.0/ext.n[0] : 0));
    printf("  => EXT at the operating point  %.2f +- %.2f\n", extMu * p, extMu * p * rel);
  } else {
    printf("  => EXT at the operating point  < %.2f (90%% CL on the MC cosmic pass fraction, 2.30/%ld)\n",
           extMu * 2.30 / std::max(nMu, 1L), nMu);
  }
  printf("  (MC background at the operating point %.2f; its cosmic-dominated part %.2f, %ld raw)\n",
         mcAll.w[NSTEP-1], mcCos.w[NSTEP-1], nOp);

  // Finer factorisation: cosmic pass fractions up to the last step that still holds
  // >= minCos MC cosmic events, then the ALL-background conditional pass fraction for
  // the remaining steps. ASSUMPTION: past that step a cosmic-induced candidate passes
  // the island and photon-gap tests at the same rate as the neutrino background.
  const long minCos = 5;
  int kc = 0; for (int k = 1; k < NSTEP; ++k) if (mcCos.n[k] >= minCos) kc = k;
  if (kc < NSTEP - 1 && mcAll.w[kc] > 0) {
    const double pc = (double)mcCos.n[kc] / nMu;
    const double pa = mcAll.w[NSTEP-1] / mcAll.w[kc];
    // binomial errors on the raw counts (bkg weights are uniform)
    const double nAk = mcAll.w2[kc] > 0 ? mcAll.w[kc]*mcAll.w[kc]/mcAll.w2[kc] : 0;   // effective raw
    const double rel = std::sqrt((ext.n[0] ? 1.0/ext.n[0] : 0) + (1-pc)/(pc*nMu)
                                 + (nAk > 0 ? (1-pa)/(pa*nAk) : 0));
    printf("\n  finer factorisation (cosmic MC to \"%s\", %ld events; then all-bkg %.4f):\n",
           STEP[kc], mcCos.n[kc], pa);
    printf("  => EXT at the operating point  %.2f +- %.2f   [assumes the remaining steps pass cosmics\n"
           "     at the neutrino-background rate; MC cosmic alone gives 0/%ld there]\n",
           extMu * pc * pa, extMu * pc * pa * rel, mcCos.n[kc]);
  }
}
