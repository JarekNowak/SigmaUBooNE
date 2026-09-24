#ifndef LambdaStyle_h
#define LambdaStyle_h
// ============================================================================
//  LambdaStyle -- the physics pieces of the Lambda CCQE selection
//  (Thorpe & Nowak, Hyperon_Technote_v2), adapted for Sigma0 -> Lambda gamma.
//
//  WHY a second selection exists. Compared at equal efficiency (6.7% vs 6.8%)
//  their analysis reaches 91% purity where ours reaches 2.8%. The gap is not the
//  operating point -- stacking their whole cut sequence onto our events tops out
//  near 1.9%. Three causes were measured:
//    1. our preselection destroys the displaced vertex, their main purifier
//       (Lambda flight > 2 cm: 55.9% of signal before preselection, 43% after);
//    2. we have no equivalent of their pair precuts;
//    3. our signal needs an extra object, the photon.
//  This header addresses (1) and (2). (3) is intrinsic and handled by treating
//  the photon as an ADDITION to their chain rather than rebuilding around it.
//
//  THE ARCHITECTURAL FIX. PidAssign.h scores each track for a role using, among
//  other things, log10(distance to the PRIMARY vertex), with a coefficient of
//  -0.2758 on the pion. That actively penalises displaced tracks -- precisely the
//  signature of a Lambda that flew. Their pair selector never uses the distance
//  to the primary vertex. Its only topological input is the separation between
//  the two track STARTS: a displaced Lambda still has its proton and pion
//  starting close to EACH OTHER, just far from the neutrino. Everything here is
//  built on that quantity instead.
// ============================================================================
#include "TVector3.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace LambdaStyle {

static const double M_PROTON = 0.938272;
static const double M_PION   = 0.139570;
static const double M_LAMBDA = 1.115683;
static const double M_SIGMA0 = 1.192642;

// ── Range-based momenta ─────────────────────────────────────────────────────
//  Proton: their parameterisation, E_k = 29.9317 * R^0.586304 MeV, validated
//  against our truth (383 MeV/c predicted vs 410 measured at R = 4.8 cm).
//
//  Pion: deliberately NOT their formula. Fitting p(R) on our own truth removes
//  the -20% muon-hypothesis bias but roughly DOUBLES the Lambda mass width
//  (7.8 -> 12.5 MeV half-IQR), because the pion carries most of the 38 MeV decay
//  energy and rescaling it amplifies its spread. The biased-but-tight
//  muon-hypothesis range momentum gives the sharper peak, which is what a mass
//  window or a BDT actually uses. Pass trk_range_muon_mom straight through.
inline double protonMomFromRange(double R_cm){
  if(!(R_cm > 0)) return -1.0;
  const double Ek = 29.9317 * std::pow(R_cm, 0.586304) / 1000.0;   // GeV
  return std::sqrt(Ek*Ek + 2.0*M_PROTON*Ek);
}

// ── Secondary vertex: fit a V to the two decay tracks ────────────────────────
//  Reconstructed tracks do not meet at a point, so the decay vertex needs
//  extrapolating. Their fitter varies a vertex v and two arm directions to match
//  the tracks. We fit v only, holding the arms at the measured track directions:
//  that captures the vertex (which is what the decay length and pointing angle
//  need) and stays robust on the short tracks this channel produces.
//
//  Metric: for points sampled along each track, the distance to the
//  corresponding arm (a ray from v). Points lying UPSTREAM of v are penalised by
//  `upstreamPenalty`, since the decay products cannot precede the decay.
struct VFitResult {
  TVector3 vertex;
  double   chi2       = -1.0;   // mean squared distance, cm^2
  bool     ok         = false;
};

inline double vFitMetric(const TVector3& v,
                         const TVector3& s1, const TVector3& d1, double L1,
                         const TVector3& s2, const TVector3& d2, double L2,
                         int nSample, double upstreamPenalty)
{
  double sum = 0.0; int n = 0;
  for(int leg = 0; leg < 2; ++leg){
    const TVector3& s = leg ? s2 : s1;
    const TVector3& d = leg ? d2 : d1;
    const double    L = leg ? L2 : L1;
    for(int i = 0; i < nSample; ++i){
      const TVector3 x = s + d * (L * i / double(std::max(1, nSample - 1)));
      const TVector3 w = x - v;
      const double   t = w.Dot(d);              // projection onto the arm
      double dist;
      if(t >= 0.0) dist = (w - d*t).Mag();      // downstream: perpendicular distance
      else         dist = w.Mag() * upstreamPenalty;   // upstream: penalised
      sum += dist*dist; ++n;
    }
  }
  return n ? sum / n : 1e9;
}

inline VFitResult fitV(const TVector3& s1, const TVector3& d1, double L1,
                       const TVector3& s2, const TVector3& d2, double L2,
                       int nSample = 12, double upstreamPenalty = 3.0)
{
  VFitResult r;
  if(d1.Mag() <= 0 || d2.Mag() <= 0) return r;
  const TVector3 u1 = d1.Unit(), u2 = d2.Unit();
  // seed at the midpoint of the two track starts (their initial value)
  TVector3 v = 0.5 * (s1 + s2);
  double best = vFitMetric(v, s1,u1,L1, s2,u2,L2, nSample, upstreamPenalty);
  // coordinate descent with a shrinking step: robust and cheap, no minimiser
  double step = 4.0;                                   // cm
  for(int iter = 0; iter < 40 && step > 0.02; ++iter){
    bool improved = false;
    for(int ax = 0; ax < 3; ++ax){
      for(int sgn = -1; sgn <= 1; sgn += 2){
        TVector3 trial = v;
        trial[ax] += sgn * step;
        const double m = vFitMetric(trial, s1,u1,L1, s2,u2,L2, nSample, upstreamPenalty);
        if(m < best){ best = m; v = trial; improved = true; }
      }
    }
    if(!improved) step *= 0.5;
  }
  r.vertex = v; r.chi2 = best; r.ok = true;
  return r;
}

// ── Decay observables built on the fitted vertex ─────────────────────────────
//  decayLength : |v - primary vertex|
//  alpha       : angle between the reconstructed Lambda momentum and the line
//                from the primary vertex to the decay vertex. If the vertex is
//                right and the decay real, the Lambda must point along that line.
//                This is their strongest decay variable and it only works once
//                the vertex is fitted rather than guessed.
struct DecayObs {
  double mass        = -1.0;
  double decayLength = -1.0;
  double alpha       = -2.0;   // radians; -2 flags "undefined"
  double cosAlpha    = -2.0;
  double vfitChi2    = -1.0;
  TVector3 lambdaMom;
  bool   ok          = false;
};

inline DecayObs decayObservables(const TVector3& primVtx,
                                 const TVector3& pStart, const TVector3& pDir, double pLen,
                                 const TVector3& iStart, const TVector3& iDir, double iLen,
                                 double pMom, double iMom)
{
  DecayObs o;
  if(pDir.Mag() <= 0 || iDir.Mag() <= 0) return o;
  if(!(pMom > 0) || !(iMom > 0)) return o;
  const TVector3 up = pDir.Unit(), ui = iDir.Unit();
  const TVector3 P  = up*pMom + ui*iMom;
  const double   E  = std::sqrt(pMom*pMom + M_PROTON*M_PROTON)
                    + std::sqrt(iMom*iMom + M_PION*M_PION);
  const double m2 = E*E - P.Mag2();
  if(m2 > 0) o.mass = std::sqrt(m2);
  o.lambdaMom = P;

  const VFitResult vf = fitV(pStart, up, pLen, iStart, ui, iLen);
  if(!vf.ok) return o;
  o.vfitChi2 = vf.chi2;
  const TVector3 flight = vf.vertex - primVtx;
  o.decayLength = flight.Mag();
  if(flight.Mag() > 0 && P.Mag() > 0){
    o.cosAlpha = P.Unit().Dot(flight.Unit());
    o.cosAlpha = std::max(-1.0, std::min(1.0, o.cosAlpha));
    o.alpha    = std::acos(o.cosAlpha);
  }
  o.ok = true;
  return o;
}

// ── Pair selector ────────────────────────────────────────────────────────────
//  Their 7 inputs, in their order. Note what is ABSENT: any distance to the
//  primary vertex. `sep` is the separation of the two track starts.
//
//  `threePlaneDedx` is their short-track particle identifier (a wire-angle
//  weighted mean dE/dx). hypana provides trk_three_plane_dedx directly; where it
//  is unavailable pass -1 and the corresponding terms drop out.
struct PairFeat {
  double sep      = -1;   // |start_p - start_pi|, cm
  double scoreP   = -1;   // Pandora track/shower score, proton leg
  double scorePi  = -1;
  double dedxP    = -1;   // three-plane mean dE/dx, proton leg
  double dedxPi   = -1;
  double llrP     = -9;   // LLR PID, proton leg
  double llrPi    = -9;
};

//  Precuts. Their purpose is to let an event be REJECTED outright when no pair
//  of tracks is decay-like, which the original ranking could not do: with three
//  tracks and one muon, the remaining two were always nominated however
//  implausible. Thresholds are ours, chosen to be loose.
struct PairCuts {
  double maxLlrProton = 0.30;   // proton leg must not look minimum-ionising
  double minLlrPion   = -0.20;  // pion leg must not look proton-like
  double maxSep       = 8.0;    // cm, the two starts share a decay vertex
};

inline bool passPrecuts(const PairFeat& f, const PairCuts& c){
  if(f.sep < 0 || f.sep > c.maxSep)        return false;
  if(f.llrP  > c.maxLlrProton)             return false;
  if(f.llrPi < c.minLlrPion)               return false;
  return true;
}

//  Provisional pair score, used until the pair-level BDT is trained by
//  DumpPairs.C + a TMVA job. Built from the same 7 quantities and the physics
//  they encode: the two starts should be close; the proton leg should be
//  heavily ionising and proton-like; the pion leg minimum-ionising and not
//  proton-like. Deliberately simple and monotonic so its behaviour is legible;
//  replace wholesale once the BDT exists.
inline double provisionalPairScore(const PairFeat& f){
  double s = 0.0;
  s -= 0.60 * f.sep;                                    // shared decay vertex
  s -= 2.00 * f.llrP;                                   // proton leg: low LLR
  s += 2.00 * f.llrPi;                                  // pion leg: high LLR
  if(f.dedxP  > 0) s += 0.35 * std::min(f.dedxP, 12.0); // proton leg ionises
  if(f.dedxPi > 0) s -= 0.35 * std::min(f.dedxPi, 12.0);// pion leg does not
  if(f.scoreP  > 0) s += 0.50 * f.scoreP;               // both should be track-like
  if(f.scorePi > 0) s += 0.50 * f.scorePi;
  return s;
}

} // namespace LambdaStyle
#endif
