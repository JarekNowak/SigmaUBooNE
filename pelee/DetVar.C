// ============================================================================
//  DetVar.C -- detector systematics of the PeLEE chain from the run 3 RHC
//  detector-variation samples, compared with their CV on the events they share.
//
//  Input: the `dv` per-event trees that SigmaSelection.C writes in detector-
//  variation mode (its detvarFile argument), one file per variation:
//    SigmaHists_dv_<name>.root
//  The variations re-simulate the SAME generated events as the CV, but some lost
//  files, so each one is compared with the CV on the (run, sub, evt, E_nu) it
//  shares with it. Over the shared events the yields are strongly correlated; the
//  statistical error on var/CV comes from the events whose selection changed (both
//  sides take the CV event's weight, so the POT scales of the files drop out):
//    var(D) = sum_events (w*pass_var - w*pass_cv)^2,  D = Y_var - Y_cv.
//  CAUTION: these samples are on the OLD NuMI flux and their Sigma0 rate per POT
//  is 1.5x the current-flux overlay's. Only the ratios var/CV are used.
//  Each variation is taken as a one-sided 1-sigma shift and the |var/CV - 1| are
//  added in quadrature.
//
//  Run:  root -l -b -q DetVar.C
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include <unordered_map>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

namespace {
struct Ev { int sig, cmax; float w; };
struct Key { int r, s, e; float E;
  bool operator==(const Key& o) const { return r==o.r && s==o.s && e==o.e && E==o.E; } };
struct KeyHash { size_t operator()(const Key& k) const {
  return std::hash<long long>()(((long long)k.r << 40) ^ ((long long)k.s << 20) ^ k.e)
         ^ std::hash<float>()(k.E); } };
using EvMap = std::unordered_map<Key, Ev, KeyHash>;

bool load(const char* fn, EvMap& m, long& dup){
  TFile f(fn); TTree* t = nullptr; f.GetObject("dv", t);
  if (!t) { printf("  %s: no dv tree\n", fn); return false; }
  int r, s, e, sig, cmax; float E, w;
  t->SetBranchAddress("run",&r);  t->SetBranchAddress("sub",&s);  t->SetBranchAddress("evt",&e);
  t->SetBranchAddress("nuE",&E);  t->SetBranchAddress("isSig",&sig);
  t->SetBranchAddress("cmax",&cmax); t->SetBranchAddress("w",&w);
  m.reserve(t->GetEntries());
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    if (!m.emplace(Key{r,s,e,E}, Ev{sig,cmax,w}).second) dup++;
  }
  return true;
}
}

void DetVar(const char* prefix = "SigmaHists_dv_")
{
  const std::vector<std::string> VARS = { "LYAttenuation","LYDown","LYRayleigh","Recombination",
    "SCE","WireModThetaXZ","WireModThetaYZ","WireModX","WireModYZ" };
  // stages reported: index into the SigmaSelection cut flow
  const int NS = 3; const int ST[NS] = { 4, 7, 8 };
  const char* STN[NS] = { "MuonCandidate", "NTracksCut", "BDTG" };

  EvMap cv; long dupCV = 0;
  if (!load(Form("%sCV.root", prefix), cv, dupCV)) return;
  printf("\n  CV: %zu events (%ld repeated keys dropped)\n", cv.size(), dupCV);

  double q2[NS][2] = {{0}};   // quadrature sums of |ratio-1|, [stage][sig,bkg]
  double s2[NS][2] = {{0}};   // the same for the statistical errors: what MC noise alone gives
  printf("\n  %-15s %8s |", "variation", "shared");
  for (int k = 0; k < NS; k++) printf("  %-24s|", STN[k]);
  printf("\n  %-15s %8s |", "", "");
  for (int k = 0; k < NS; k++) printf("  %11s %11s |", "S var/CV", "B var/CV");
  printf("\n");

  for (const auto& v : VARS) {
    EvMap var; long dupV = 0;
    if (!load(Form("%s%s.root", prefix, v.c_str()), var, dupV)) continue;
    // yields over shared events: [stage][sig/bkg]
    double ycv[NS][2] = {{0}}, yv[NS][2] = {{0}}, vd[NS][2] = {{0}};
    long shared = 0, sigFlip = 0;
    for (const auto& p : var) {
      auto it = cv.find(p.first);
      if (it == cv.end()) continue;
      shared++;
      const Ev& a = it->second; const Ev& b = p.second;
      if (a.sig != b.sig) sigFlip++;
      const int c = a.sig ? 0 : 1;           // truth class from the CV (identical by construction)
      for (int k = 0; k < NS; k++) {
        // the CV weight on both sides: same generated event, same truth weights, and
        // the variation's own w carries ITS POT scale, which differs where files were
        // lost (WireModThetaYZ keeps 81% of the events -> a spurious x1.23)
        const double pa = (a.cmax >= ST[k]) ? a.w : 0.0, pb = (b.cmax >= ST[k]) ? a.w : 0.0;
        ycv[k][c] += pa; yv[k][c] += pb; vd[k][c] += (pb - pa) * (pb - pa);
      }
    }
    printf("  %-15s %8ld |", v.c_str(), shared);
    for (int k = 0; k < NS; k++) {
      for (int c = 0; c < 2; c++) {
        const double r = ycv[k][c] > 0 ? yv[k][c] / ycv[k][c] : 0;
        const double er = ycv[k][c] > 0 ? std::sqrt(vd[k][c]) / ycv[k][c] : 0;
        printf(" %5.3f+-%4.3f", r, er);
        q2[k][c] += (r - 1) * (r - 1); s2[k][c] += er * er;
      }
      printf(" |");
    }
    if (sigFlip) printf("  (%ld truth-class mismatches)", sigFlip);
    printf("\n");
  }
  //  A variation consistent with no effect still adds its statistical scatter, so the
  //  quadrature sum is compared with the sum of the statistical errors: where the two
  //  are equal the detector effect is not resolved and the noise level is a bound.
  printf("\n  quadrature sum |var/CV-1|   (sum of statistical errors, i.e. MC noise alone):\n");
  for (int k = 0; k < NS; k++)
    printf("    %-14s  signal %5.1f%% (%4.1f%%)   background %5.1f%% (%4.1f%%)\n", STN[k],
           100 * std::sqrt(q2[k][0]), 100 * std::sqrt(s2[k][0]),
           100 * std::sqrt(q2[k][1]), 100 * std::sqrt(s2[k][1]));
  // CV raw content at each stage, for scale
  long n[NS][2] = {{0}};
  for (const auto& p : cv) for (int k = 0; k < NS; k++) if (p.second.cmax >= ST[k]) n[k][p.second.sig ? 0 : 1]++;
  printf("\n  CV raw events: ");
  for (int k = 0; k < NS; k++) printf("%s S %ld B %ld   ", STN[k], n[k][0], n[k][1]);
  printf("\n");
}
