// CompareBDT.C -- out-of-sample comparison of two BDTInput files (cross-validated BdtScore):
// weighted AUC, B/B0 at fixed signal efficiency (all bkg and per component), FoM peaks.
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
struct E { float s, w; int smp; bool sig; };
static std::vector<E> load(const char* fn){
  std::vector<E> v; TFile f(fn);
  for (int isSig = 0; isSig < 2; ++isSig) {
    TTree* t = (TTree*)f.Get(isSig ? "TreeS" : "TreeB");
    float w, sc, smp; t->SetBranchAddress("Weight",&w); t->SetBranchAddress("BdtScore",&sc);
    t->SetBranchAddress("Sample",&smp);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); v.push_back({sc, w, (int)smp, (bool)isSig}); }
  }
  return v;
}
static void report(const char* fn, const char* tag){
  auto v = load(fn);
  std::sort(v.begin(), v.end(), [](const E& a, const E& b){ return a.s > b.s; });
  double S0 = 0, B0 = 0, C0[6] = {0};
  for (auto& e : v) { if (e.sig) S0 += e.w; else { B0 += e.w; C0[e.smp] += e.w; } }
  // AUC: integral of signal eff d(bkg eff) going down in score
  double s = 0, b = 0, auc = 0, C[6] = {0};
  const double effs[3] = {0.2, 0.4, 0.5};   int ie = 0;   // ascending: signal eff grows as the cut drops
  double bAt[3] = {0}, cAt[3][6] = {{0}}, cutAt[3] = {0};
  double best[3] = {-1,-1,-1}, bc[3] = {0}, bs[3] = {0}, bb[3] = {0};
  const double fsys[3] = {0.0, 0.02, 0.05};
  for (size_t i = 0; i < v.size(); ++i) {
    const E& e = v[i];
    if (e.sig) { s += e.w; } else { auc += (s / S0) * (e.w / B0); b += e.w; C[e.smp] += e.w; }
    while (ie < 3 && s / S0 >= effs[ie]) { bAt[ie] = b / B0; for (int k = 0; k < 6; k++) cAt[ie][k] = C0[k] > 0 ? C[k] / C0[k] : 0; cutAt[ie] = e.s; ie++; }
    if (i + 1 < v.size() && v[i + 1].s == e.s) continue;
    if (s > 0 && b > 0 && s / S0 > 0.02)
      for (int k = 0; k < 3; k++) { const double fm = s / std::sqrt(s + b + fsys[k]*fsys[k]*b*b);
        if (fm > best[k]) { best[k] = fm; bc[k] = e.s; bs[k] = s; bb[k] = b; } }
  }
  printf("\n  === %s (%s) ===  S0 %.3f  B0 %.1f (hyp %.1f ovl %.1f dirt %.1f EXT %.1f)\n", tag, fn, S0, B0, C0[1], C0[2], C0[3], C0[5]);
  printf("  AUC %.4f\n", auc);
  for (int k = 0; k < 3; k++)
    printf("  at %2.0f%% signal eff (BDTG > %+.3f): B/B0 %.4f   overlay %.4f  dirt %.4f  EXT %.4f  hyperons %.4f\n",
           100*effs[k], cutAt[k], bAt[k], cAt[k][2], cAt[k][3], cAt[k][5], cAt[k][1]);
  for (int k = 0; k < 3; k++)
    printf("  FoM peak, f=%.0f%%: BDTG > %+.3f  S %.3f  B %.1f  purity %.3f%%  FoM %.4f\n",
           100*fsys[k], bc[k], bs[k], bb[k], 100*bs[k]/(bs[k]+bb[k]), best[k]);
}
void CompareBDT(const char* a = "BDTInput.root.oldbdt", const char* b = "BDTInput.root"){
  report(a, "A");
  report(b, "B");
}
