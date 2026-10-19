// ============================================================================
//  Sensitivity.C -- expected 90% CL upper limits on the Sigma0 rate for the
//  PeLEE BDT chain (SigmaSelection.C) and the Lambda-style chain
//  (SigmaLambdaSel.C), in units of each generator's prediction: mu = 1 is the
//  model. Replaces the BNB-era ../CrossSection.C for planning purposes.
//
//  One-bin counting experiment, n ~ Pois(mu*S*e + b):
//    b  lognormal, mean B, fractional width sqrt(rMC^2 + f^2): rMC the MC (and
//       EXT) statistics of B, f the background systematic, which is scanned
//    e  lognormal, mean 1, fractional width sqrt(g^2 + rS^2): g the signal
//       normalisation systematic (flux, efficiency, POT), rS the signal MC stats
//  The limit is Bayesian with a flat prior in mu >= 0 and b, e marginalised,
//  cross-checked with hybrid CLs (CLs = P(N<=n|mu) / P(N<=n|0), same
//  marginalisation). For a known b the two are identical (Helene); the
//  self-test at start-up checks both against the PDG table.
//  Expected limit = the limit at the background-only median n, with the 16-84%
//  band from the quantiles of Pois(B) (the limit is monotonic in n). The BDT-cut
//  scan uses the Asimov n = B instead, which is smooth in the cut.
//  Z_A is the median discovery significance with an uncertain background
//  (Cowan, Cranmer, Gross, Vitells, eq. 97 of arXiv:1007.1727).
//
//  mu_UL = n_UL(B, sigma_b) / S, so the generator only rescales the limit: the
//  cut that minimises it does not depend on the signal normalisation.
//
//  Inputs: BDTInput.root (cross-validated BdtScore, S and B per event) and
//  SigmaHists.root (GenieSignal, the overlay's own signal through the same
//  selection) for the PeLEE chain; the Lambda-style numbers are the README
//  operating point and are typed in below.
//
//  Run:  root -l -b -q Sensitivity.C+
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "Math/SpecFuncMathCore.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

namespace {
const double CL = 0.90;

// ── Lambda-style operating point (README "Operating point", ExtFactor.C) ──────
//  two planes, photon gap on. GENIE from the hyperon sample, which holds 1.32x
//  the plain overlay's in-FV Sigma0 per POT (README "The three signal samples
//  disagree"); the GENIE column below uses the overlay normalisation, the one
//  the PeLEE chain is measured with, so the two chains face the same model.
const double LS_SG_HYP = 1.33,  LS_SG_ERR = 0.16;   // 69 raw
const double LS_SN     = 0.142, LS_SN_ERR = 0.009;  // 245 raw, flux-weighted
const double LS_BMC    = 3.76,  LS_BMC_ERR = 0.78;  // 23 raw
const double LS_EXT    = 0.62,  LS_EXT_ERR = 0.27;  // factorised
//  NuWro Sigma0 lost to FSI (mostly to Lambda): IslandProcess.C on sig_nuwro.list
//  with signalOnly=false, flux-weighted, 37 raw. The PeLEE B carries the same
//  component as "hyperons".
const double LS_HYP    = 0.0207, LS_HYP_ERR = 0.0034;
const double GENIE_HYP_OVER_OVERLAY = 110.1 / 83.3;
//  Dirt is not in hypana format. PeLEE dirt/overlay is 8.7% at MuonCandidate
//  and 0.6% after the BDT; the upper end bounds the Lambda-style dirt.
const double LS_DIRT_MAX = 0.087 * LS_BMC;

// Lognormal nuisance, mean 1, fractional width r: 41 nodes in z = -5..5.
struct Nuis { std::vector<double> x, w; };
Nuis lognormal(double r){
  Nuis n;
  if (r <= 0) { n.x = {1.0}; n.w = {1.0}; return n; }
  const double s = std::sqrt(std::log(1 + r * r));
  double sum = 0;
  for (int i = 0; i <= 40; ++i) {
    const double z = -5.0 + 0.25 * i, w = std::exp(-0.5 * z * z);
    n.x.push_back(std::exp(s * z - 0.5 * s * s)); n.w.push_back(w); sum += w;
  }
  for (auto& w : n.w) w /= sum;
  return n;
}

struct Model { double B; Nuis b, e; };
Model model(double B, double rB, double rS){ return { B, lognormal(rB), lognormal(rS) }; }

// marginal likelihood of n at s signal events, relative to Pois(n|n)
double like(double n, double s, const Model& m){
  const double ref = n > 0 ? n * std::log(n) - n : 0.0;
  double L = 0;
  for (size_t i = 0; i < m.b.x.size(); ++i)
    for (size_t j = 0; j < m.e.x.size(); ++j) {
      const double lam = s * m.e.x[j] + m.B * m.b.x[i];
      L += m.b.w[i] * m.e.w[j] * std::exp((n > 0 ? n * std::log(lam) : 0.0) - lam - ref);
    }
  return L;
}

// marginal P(N <= n), continuous in n through the regularised gamma function
double cdf(double n, double s, const Model& m){
  double P = 0;
  for (size_t i = 0; i < m.b.x.size(); ++i)
    for (size_t j = 0; j < m.e.x.size(); ++j)
      P += m.b.w[i] * m.e.w[j] * ROOT::Math::inc_gamma_c(n + 1, s * m.e.x[j] + m.B * m.b.x[i]);
  return P;
}

// Bayesian upper limit in signal events, flat prior
double ulBayes(double n, const Model& m){
  // start at the -3 sigma signal-efficiency node (node 8), double until the
  // posterior at the edge is negligible
  double smax = (n + 10 * std::sqrt(n + 1) + 10) / m.e.x[std::min<size_t>(8, m.e.x.size() - 1)];
  const int NG = 4000;
  for (int iter = 0; iter < 10; ++iter, smax *= 2) {
    std::vector<double> L(NG + 1);
    double Lmax = 0;
    for (int k = 0; k <= NG; ++k) { L[k] = like(n, smax * k / NG, m); Lmax = std::max(Lmax, L[k]); }
    if (L[NG] > 1e-7 * Lmax) continue;
    std::vector<double> C(NG + 1, 0.0);
    for (int k = 1; k <= NG; ++k) C[k] = C[k - 1] + 0.5 * (L[k] + L[k - 1]);
    const double target = CL * C[NG];
    for (int k = 1; k <= NG; ++k)
      if (C[k] >= target) return smax * (k - 1 + (target - C[k - 1]) / (C[k] - C[k - 1])) / NG;
  }
  return -1;
}

// hybrid CLs upper limit in signal events
double ulCLs(double n, const Model& m){
  const double cb = cdf(n, 0, m);
  double lo = 0, hi = n + 10 * std::sqrt(n + 1) + 10;
  while (cdf(n, hi, m) / cb > 1 - CL) hi *= 2;
  for (int it = 0; it < 60; ++it) { const double mid = 0.5 * (lo + hi); (cdf(n, mid, m) / cb > 1 - CL ? lo : hi) = mid; }
  return 0.5 * (lo + hi);
}

int quantile(double B, double q){ int n = 0; while (ROOT::Math::inc_gamma_c(n + 1, B) < q) ++n; return n; }

double zA(double s, double b, double sb){
  if (s <= 0 || b <= 0) return 0;
  if (sb <= 0) return std::sqrt(2 * ((s + b) * std::log(1 + s / b) - s));
  const double v = sb * sb;
  return std::sqrt(2 * ((s + b) * std::log((s + b) * (b + v) / (b * b + (s + b) * v))
                        - b * b / v * std::log(1 + v * s / (b * (b + v)))));
}

// expected limits in signal events: 16/50/84% band (Bayes) and the median CLs
struct Exp { double lo, med, hi, cls; };
Exp expected(double B, double rB, double rS){
  const Model m = model(B, rB, rS);
  const int n16 = quantile(B, 0.1587), n50 = quantile(B, 0.5), n84 = quantile(B, 0.8413);
  return { ulBayes(n16, m), ulBayes(n50, m), ulBayes(n84, m), ulCLs(n50, m) };
}

struct Ev { float s, w; int smp; };
struct Yield { double S = 0, S2 = 0, B = 0, B2 = 0, C[6] = {0}; long rawS = 0, raw[6] = {0}; };
//  Dirt and EXT carry 0.66 and 0.75 events per raw event, so at tight cuts a
//  component can run out of raw events while its true yield is not small: a
//  component with none left gets the error of one event of its own weight.
Yield yieldAt(const std::vector<Ev>& sig, const std::vector<Ev>& bkg, double cut){
  Yield y; double wMax[6] = {0};
  for (auto& e : sig) if (e.s > cut) { y.S += e.w; y.S2 += e.w * e.w; y.rawS++; }
  for (auto& e : bkg) {
    wMax[e.smp] = std::max(wMax[e.smp], (double)e.w);
    if (e.s > cut) { y.B += e.w; y.B2 += e.w * e.w; y.C[e.smp] += e.w; y.raw[e.smp]++; }
  }
  for (int k : {3, 5}) if (y.raw[k] == 0) y.B2 += wMax[k] * wMax[k];
  return y;
}

struct Row { const char* name; double sN, rSN, sG, rSG, B, rB; };

// nPredN/nPredG > 0: also print mu_UL x the model's in-FV signal count, the
// limit on the number of signal events, which is model-independent as far as
// the efficiency is
void printChain(const Row& r, double g, const double* fs, int nf, double nPredN = 0, double nPredG = 0){
  printf("\n  %s   S(NuWro) %.3f  S(GENIE) %.3f  B %.2f (MC stat %.1f%%)  g %.0f%%\n",
         r.name, r.sN, r.sG, r.B, 100 * r.rB, 100 * g);
  printf("    %4s  %-26s %9s %8s   %7s %7s", "f", "mu_UL NuWro: median [16,84]", "GENIE", "CLs/Bay", "Z NuWro", "Z GENIE");
  if (nPredN > 0) printf("   N_UL NuWro / GENIE");
  printf("\n");
  for (int k = 0; k < nf; ++k) {
    const double rB = std::sqrt(r.rB * r.rB + fs[k] * fs[k]);
    const Exp eN = expected(r.B, rB, std::sqrt(g * g + r.rSN * r.rSN));
    const Exp eG = expected(r.B, rB, std::sqrt(g * g + r.rSG * r.rSG));
    printf("    %3.0f%%  %7.1f  [%6.1f, %6.1f]     %9.2f %8.3f   %7.3f %7.3f", 100 * fs[k],
           eN.med / r.sN, eN.lo / r.sN, eN.hi / r.sN, eG.med / r.sG, eN.cls / eN.med,
           zA(r.sN, r.B, rB * r.B), zA(r.sG, r.B, rB * r.B));
    if (nPredN > 0) printf("   %6.0f / %6.0f", eN.med / r.sN * nPredN, eG.med / r.sG * nPredG);
    printf("\n");
  }
}
}

void Sensitivity(const char* bdtInput = "BDTInput.root", const char* hists = "SigmaHists.root",
                 double bdtCut = 0.759, double g = 0.20)
{
  // ── 1. self-test against the PDG table (90% CL, b = 0: 2.30 3.89 5.32 6.68 7.99)
  printf("\n== self-test, b = 0 (PDG: 2.30 3.89 5.32 6.68 7.99) and b = 3, n = 3 ==\n  Bayes:");
  const Model m0 = model(0, 0, 0);
  for (int n = 0; n <= 4; ++n) printf(" %.3f", ulBayes(n, m0));
  printf("\n  CLs:  ");
  for (int n = 0; n <= 4; ++n) printf(" %.3f", ulCLs(n, m0));
  const Model m3 = model(3, 0, 0);
  printf("\n  b=3 n=3: Bayes %.3f  CLs %.3f (identical for known b)\n", ulBayes(3, m3), ulCLs(3, m3));
  const Model m3s = model(3, 0.3, 0.2);
  printf("  b=3 +-30%%, e +-20%%, n=3: Bayes %.3f  CLs %.3f\n", ulBayes(3, m3s), ulCLs(3, m3s));

  // ── 2. PeLEE chain: events from BDTInput.root, GENIE signal from SigmaHists.root
  std::vector<Ev> sig, bkg;
  {
    TFile f(bdtInput);
    for (int isS = 0; isS < 2; ++isS) {
      TTree* t = (TTree*)f.Get(isS ? "TreeS" : "TreeB");
      if (!t) { printf("no %s in %s\n", isS ? "TreeS" : "TreeB", bdtInput); return; }
      float w, sc, smp;
      t->SetBranchAddress("Weight", &w); t->SetBranchAddress("BdtScore", &sc); t->SetBranchAddress("Sample", &smp);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); (isS ? sig : bkg).push_back({sc, w, (int)smp}); }
    }
  }
  double genie = -1, genieErr = 0, nPredG = 0, nPredN = 0;   // nPred: in-FV signal, AllEvents
  {
    TFile f(hists);
    TH1D* h = (TH1D*)f.Get("GenieSignal");
    if (h) { const int b = h->GetXaxis()->FindBin("BDTG"); genie = h->GetBinContent(b); genieErr = h->GetBinError(b);
             nPredG = h->GetBinContent(h->GetXaxis()->FindBin("AllEvents")); }
    TH1D* hn = (TH1D*)f.Get("AllEvents/NTracksSigmaAllEvents");
    if (hn) nPredN = hn->Integral(0, hn->GetNbinsX() + 1);
  }
  if (genie <= 0) { printf("no GenieSignal in %s: rerun SigmaSelection.C\n", hists); return; }
  const Yield y = yieldAt(sig, bkg, bdtCut);
  const double ratioGN = genie / y.S;   // valid at the cut SigmaSelection.C was run with
  printf("\n== PeLEE BDT chain, BDTG > %.3f ==\n", bdtCut);
  printf("  S NuWro %.3f (%ld raw)  S GENIE overlay %.3f +- %.3f  ->  GENIE/NuWro %.2f\n",
         y.S, y.rawS, genie, genieErr, ratioGN);
  printf("  in-FV signal at AllEvents: NuWro %.2f  GENIE %.2f\n", nPredN, nPredG);
  printf("  B %.2f: hyperons %.2f (%ld)  overlay %.2f (%ld)  dirt %.2f (%ld)  EXT %.2f (%ld) raw; MC stat %.2f%%\n",
         y.B, y.C[1], y.raw[1], y.C[2], y.raw[2], y.C[3], y.raw[3], y.C[5], y.raw[5], 100 * std::sqrt(y.B2) / y.B);

  const double sGls = LS_SG_HYP / GENIE_HYP_OVER_OVERLAY;
  const double bLs = LS_BMC + LS_EXT + LS_HYP;
  const double rLs = std::sqrt(LS_BMC_ERR * LS_BMC_ERR + LS_EXT_ERR * LS_EXT_ERR + LS_HYP_ERR * LS_HYP_ERR) / bLs;
  printf("\n== Lambda-style chain (two planes, gap) ==\n");
  printf("  S NuWro %.3f  S GENIE %.3f (hyperon sample %.2f / %.3f)  B %.2f (%.2f MC + %.2f EXT + %.3f hyperons), MC stat %.1f%%\n",
         LS_SN, sGls, LS_SG_HYP, GENIE_HYP_OVER_OVERLAY, bLs, LS_BMC, LS_EXT, LS_HYP, 100 * rLs);

  // ── 3. expected limits per chain and background systematic
  const double fs[4] = {0.0, 0.07, 0.15, 0.30};
  printf("\n== Expected 90%% CL limits on mu (rate / model), background-only median ==\n");
  const Row pe = { "PeLEE BDT", y.S, std::sqrt(y.S2) / y.S, genie, genieErr / genie, y.B, std::sqrt(y.B2) / y.B };
  const Row ls = { "Lambda-style", LS_SN, LS_SN_ERR / LS_SN, sGls, LS_SG_ERR / LS_SG_HYP, bLs, rLs };
  printChain(pe, g, fs, 4, nPredN, nPredG);
  printChain(ls, g, fs, 4);

  // ── 4. what moves the Lambda-style limit, f = 15%
  //  At B ~ 4 the median n is an integer, so a small change in B can move the
  //  median limit either way (a lower B at the same n gives a HIGHER limit). The
  //  Asimov column (n = B) shows the smooth dependence.
  printf("\n== Lambda-style variants, f = 15%%: mu_UL (NuWro / GENIE) ==\n");
  printf("  %-44s %8s %4s   %-15s   %-15s\n", "", "B", "n50", "median", "Asimov n = B");
  auto lsVar = [&](const char* tag, double B, double rB, double gg, double sG){
    const double rBt = std::sqrt(rB * rB + 0.15 * 0.15);
    const Exp eN = expected(B, rBt, std::hypot(gg, LS_SN_ERR / LS_SN));
    const Exp eG = expected(B, rBt, std::hypot(gg, LS_SG_ERR / LS_SG_HYP));
    const double aN = ulBayes(B, model(B, rBt, std::hypot(gg, LS_SN_ERR / LS_SN)));
    const double aG = ulBayes(B, model(B, rBt, std::hypot(gg, LS_SG_ERR / LS_SG_HYP)));
    printf("  %-44s %8.2f %4d   %6.1f / %5.2f   %6.1f / %5.2f\n", tag, B, quantile(B, 0.5),
           eN.med / LS_SN, eG.med / sG, aN / LS_SN, aG / sG);
  };
  lsVar("nominal", bLs, rLs, g, sGls);
  lsVar("without the EXT estimate", LS_BMC + LS_HYP, std::hypot(LS_BMC_ERR, LS_HYP_ERR) / (LS_BMC + LS_HYP), g, sGls);
  lsVar("with dirt at its bound (+0.33 +- 0.33)", bLs + LS_DIRT_MAX,
        std::hypot(rLs * bLs, LS_DIRT_MAX) / (bLs + LS_DIRT_MAX), g, sGls);
  lsVar("no signal normalisation systematic (g = 0)", bLs, rLs, 0.0, sGls);
  lsVar("GENIE at the hyperon-sample normalisation", bLs, rLs, g, LS_SG_HYP);
  lsVar("background known exactly (rMC = 0)", bLs, 0.0, g, sGls);

  // ── 5. PeLEE BDT cut chosen by the limit, against the counting FoM
  printf("\n== PeLEE BDT cut that minimises the Asimov limit (NuWro mu; GENIE = NuWro / %.2f) ==\n", ratioGN);
  printf("  %4s   %-44s   %-30s   %-20s\n", "f", "limit-optimal: cut S B (MC stat, EXT raw) mu", "FoM-optimal: cut S B mu", "at the op point: mu");
  const double fsc[5] = {0.0, 0.05, 0.10, 0.15, 0.30};
  std::vector<double> cuts;
  for (double c = 0.50; c < 0.951; c += 0.005) cuts.push_back(c);
  std::vector<Yield> ys; for (double c : cuts) ys.push_back(yieldAt(sig, bkg, c));
  auto muAt = [&](const Yield& q, double f){
    const Model m = model(q.B, std::sqrt(q.B2 / (q.B * q.B) + f * f), std::sqrt(g * g + q.S2 / (q.S * q.S)));
    return ulBayes(q.B, m) / q.S; };
  for (double f : fsc) {
    int iBest = -1, iFom = -1; double muBest = 1e99, fomBest = -1;
    for (size_t i = 0; i < cuts.size(); ++i) {
      const Yield& q = ys[i];
      if (q.rawS < 50 || q.B <= 0) continue;
      const double mu = muAt(q, f);
      if (mu < muBest) { muBest = mu; iBest = i; }
      const double fom = q.S / std::sqrt(q.S + q.B + f * f * q.B * q.B);
      if (fom > fomBest) { fomBest = fom; iFom = i; }
    }
    const Yield &qb = ys[iBest], &qf = ys[iFom];
    printf("  %3.0f%%   %.3f %5.2f %6.1f (%4.1f%%, %3ld) %6.1f          %.3f %5.2f %6.1f %6.1f        %6.1f\n", 100 * f,
           cuts[iBest], qb.S, qb.B, 100 * std::sqrt(qb.B2) / qb.B, qb.raw[5], muBest,
           cuts[iFom], qf.S, qf.B, muAt(qf, f), muAt(y, f));
  }
}
