// ============================================================================
//  HypXsecWeight.h -- per-event weight correcting GENIE's numubar hyperon QE
//  cross section (PaisQELLambdaPXSec) to the full trace matrix element.
//  Table from HypXsecWeights.py (hyp_xsec_weights.root): w(log10 E_rest, u) per hyperon,
//  E_rest = antineutrino energy in the struck-nucleon rest frame,
//  u = (Q2 - Q2min)/(Q2max - Q2min) for a free nucleon at E_rest (clamped to [0,1]).
//  Without the struck nucleon (analysis ntuples) pass a nucleon at rest.
// ============================================================================
#pragma once
#include "TFile.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include <map>
#include <cmath>
#include <cstdio>

namespace HypXsecWeight {

static const double MMU = 0.10566;
static const double WMAX = 5.0;

inline double hyperonMass(int pdg) { return pdg == 3122 ? 1.115683 : pdg == 3212 ? 1.192642 : 1.197449; }
inline double nucleonMass(int pdg) { return pdg == 3112 ? 0.93957 : 0.93827; }

struct Table {
  std::map<int, TH2D*> h;
  bool load(const char* fn = "hyp_xsec_weights.root") {
    TFile* f = TFile::Open(fn);
    if (!f || f->IsZombie()) { printf("  HypXsecWeight: cannot open %s\n", fn); return false; }
    for (int p : {3122, 3212, 3112}) {
      TH2D* x = (TH2D*)f->Get(Form("w_%d", p));
      if (x) { x->SetDirectory(nullptr); h[p] = x; }
    }
    delete f;
    return h.size() == 3;
  }
  // E_rest and u from the free-nucleon kinematics, then bilinear interpolation.
  double weight(int pdg, double Erest, double Q2) const {
    auto it = h.find(pdg); if (it == h.end() || Erest <= 0) return 1.;
    const double Mn = nucleonMass(pdg), MY = hyperonMass(pdg);
    const double s = Mn * Mn + 2 * Mn * Erest, rs = std::sqrt(s);
    const double El = (s + MMU * MMU - MY * MY) / (2 * rs), pl2 = El * El - MMU * MMU;
    if (pl2 <= 0) return 1.;
    const double pl = std::sqrt(pl2), Ein = (s - Mn * Mn) / (2 * rs);
    const double q2lo = -MMU * MMU + 2 * Ein * (El - pl), q2hi = -MMU * MMU + 2 * Ein * (El + pl);
    double u = (Q2 - q2lo) / (q2hi - q2lo); u = std::min(std::max(u, 0.0), 1.0);
    const TH2D* t = it->second;
    double x = std::log10(Erest);
    const TAxis* ax = t->GetXaxis(); const TAxis* ay = t->GetYaxis();
    x = std::min(std::max(x, ax->GetBinCenter(1)), ax->GetBinCenter(ax->GetNbins()));
    u = std::min(std::max(u, ay->GetBinCenter(1)), ay->GetBinCenter(ay->GetNbins()));
    // Near u = 1 GENIE's Pais expression goes negative and is clipped to 0, so the ratio
    // diverges where GENIE has (almost) no events: <= 0.3% of the corrected sigma. Cap it.
    return std::min(const_cast<TH2D*>(t)->Interpolate(x, u), WMAX);
  }
  // From four-vectors: nubar, struck nucleon (as generated), mu+.
  double weight(int pdg, const TLorentzVector& nu, const TLorentzVector& nuc, const TLorentzVector& mu) const {
    TLorentzVector n = nu; n.Boost(-nuc.BoostVector());
    return weight(pdg, n.E(), -(nu - mu).M2());
  }
};

}  // namespace HypXsecWeight
