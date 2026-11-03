// ============================================================================
//  GenieHypXsecWeight.h -- per-event weight correcting GENIE's antineutrino Delta S = 1
//  quasi-elastic hyperon cross section (genie::PaisQELLambdaPXSec, QEL-CC-LAMBDA:
//  nubar p -> l+ Lambda, nubar p -> l+ Sigma0, nubar n -> l+ Sigma-) to the full
//  contraction L_mu_nu H^mu_nu with the same SU(3) form factors. See README.md.
//
//  w(E_rest, u) = [dsigma/dQ2]_correct / [dsigma/dQ2]_GENIE for a free nucleon, tabulated in
//  log10 E_rest (antineutrino energy in the struck-nucleon rest frame, 0.2-20 GeV) and
//  u = (Q2 - Q2min)/(Q2max - Q2min) at E_rest. Two tables: charged muon (nubar_mu) and
//  electron (nubar_e). Bilinear interpolation between bin centres, clamped to the table,
//  capped at WMAX (near Q2max GENIE's expression is clipped to 0 and the ratio diverges where
//  GENIE has no events: <= 0.16% of the correct cross section).
//
//  Header-only, ROOT only.
//    GenieHypXsecWeight::Weights W;  W.load("hyp_xsec_weights_mu.root", "hyp_xsec_weights_e.root");
//    double w = W.weight(nuPdg, hyperonPdg, nu, hitNucleon, lepton);    // four-vectors in GeV
//  Apply ONLY to GENIE QEL-CC-LAMBDA events (see README.md, "Which events"); returns 1 otherwise.
// ============================================================================
#pragma once
#include "TFile.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

namespace GenieHypXsecWeight {

static const double WMAX = 5.0;
inline double hyperonMass(int pdg) { return pdg == 3122 ? 1.115683 : pdg == 3212 ? 1.192642 : 1.197449; }
inline double nucleonMass(int pdg) { return pdg == 3112 ? 0.93957 : 0.93827; }
inline bool isHyperonQE(int pdg) { return pdg == 3122 || pdg == 3212 || pdg == 3112; }

// One table set (one charged-lepton mass).
struct Table {
  std::map<int, TH2D*> h;
  double ml = 0.10566;
  bool load(const char* fn, double leptonMass) {
    ml = leptonMass;
    TFile* f = TFile::Open(fn, "READ");
    if (!f || f->IsZombie()) { printf("GenieHypXsecWeight: cannot open %s\n", fn); return false; }
    for (int p : {3122, 3212, 3112}) {
      TH2D* x = (TH2D*)f->Get(Form("w_%d", p));
      if (x) { x = (TH2D*)x->Clone(Form("ghxw_%d_%g", p, leptonMass)); x->SetDirectory(nullptr); h[p] = x; }
    }
    delete f;
    return h.size() == 3;
  }
  // E_rest [GeV], Q2 = -(nu - l)^2 [GeV^2]
  double weight(int pdg, double Erest, double Q2) const {
    auto it = h.find(pdg); if (it == h.end() || Erest <= 0) return 1.;
    const double Mn = nucleonMass(pdg), MY = hyperonMass(pdg);
    const double s = Mn * Mn + 2 * Mn * Erest, rs = std::sqrt(s);
    const double El = (s + ml * ml - MY * MY) / (2 * rs), pl2 = El * El - ml * ml;
    if (pl2 <= 0) return 1.;                                   // below threshold at rest
    const double pl = std::sqrt(pl2), Ein = (s - Mn * Mn) / (2 * rs);
    const double q2lo = -ml * ml + 2 * Ein * (El - pl), q2hi = -ml * ml + 2 * Ein * (El + pl);
    double u = (Q2 - q2lo) / (q2hi - q2lo); u = std::min(std::max(u, 0.0), 1.0);
    TH2D* t = it->second;
    const TAxis* ax = t->GetXaxis(); const TAxis* ay = t->GetYaxis();
    double x = std::log10(Erest);
    x = std::min(std::max(x, ax->GetBinCenter(1)), ax->GetBinCenter(ax->GetNbins()));
    u = std::min(std::max(u, ay->GetBinCenter(1)), ay->GetBinCenter(ay->GetNbins()));
    return std::min(t->Interpolate(x, u), WMAX);
  }
};

// Both flavours; picks the table from the neutrino PDG code.
struct Weights {
  Table mu, e;
  bool haveE = false;
  bool load(const char* fmu = "hyp_xsec_weights_mu.root", const char* fe = "hyp_xsec_weights_e.root") {
    bool ok = mu.load(fmu, 0.10566);
    haveE = fe && *fe && e.load(fe, 0.000511);
    if (!haveE) printf("GenieHypXsecWeight: no electron table; nubar_e events use the muon table\n");
    return ok;
  }
  // E_rest and Q2 already computed (hit nucleon at rest: E_rest = E_nu).
  double weight(int nuPdg, int hyperonPdg, double Erest, double Q2) const {
    if (!isHyperonQE(hyperonPdg) || nuPdg >= 0) return 1.;   // GENIE has this channel for nubar only
    const Table& t = (std::abs(nuPdg) == 12 && haveE) ? e : mu;
    return t.weight(hyperonPdg, Erest, Q2);
  }
  // From four-vectors: incoming nubar, struck nucleon as generated (off-shell is fine), outgoing lepton.
  double weight(int nuPdg, int hyperonPdg, const TLorentzVector& nu, const TLorentzVector& nucleon,
                const TLorentzVector& lepton) const {
    TLorentzVector n = nu; n.Boost(-nucleon.BoostVector());
    return weight(nuPdg, hyperonPdg, n.E(), -(nu - lepton).M2());
  }
};

}  // namespace GenieHypXsecWeight
