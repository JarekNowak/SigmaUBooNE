#ifndef PidAssign_h
#define PidAssign_h
// ============================================================================
//  PidAssign -- joint proton/pion candidate assignment for the Sigma0 selection.
//
//  WHY: the old rule ranked non-muon tracks by a single monotonic discriminant
//    D = LLR/0.591 + (chi2_p - chi2_pi)/124.5,  proton = min D, pion = max D
//  and took the extremes. Measured against HyperonProduction truth on the GENIE
//  run3b sample (1133 selected signal events), that gives a pion candidate which
//  is a true pi- only 33.5% of the time. Its contaminants are photons (20.4%),
//  muons (15.4%) and protons (12.5%): shower-like and MIP-like PFPs both sit at
//  the high-D end, so a monotonic rank cannot separate them from a real pion,
//  whose ionisation sits in a BAND between minimum-ionising and proton-like.
//
//  Instead score each track separately for "is the pion" and "is the proton"
//  with a quadratic logistic model, then choose the ORDERED PAIR maximising
//  Sproton(a) + Spion(b). Selecting a better pair among the PFPs already present
//  costs no efficiency -- the event is kept either way, only the indices change.
//
//  Measured out-of-sample (2-fold by event parity), same 1133 events:
//               pion       proton     both
//    old D      33.5%      61.2%      29.7%
//    this       42.5%      62.3%      35.6%
//    ceiling    51.4%      66.2%            (events where the track is present)
//  i.e. 65% -> 83% of the achievable pion purity, at no cost in signal.
//
//  Inputs are only variables carried by BOTH ntuple formats (PeLEE trk_score_v,
//  trk_llr_pid_score_v, trk_pid_chipr_v, trk_pid_chipi_v, trk_len_v, start
//  position, hit count; and their hypana equivalents), so the same assignment
//  runs in the PeLEE selection and in the hypana studies. Adding the hypana-only
//  wiggliness and three-plane dE/dx was measured to gain nothing (both 35.8%).
//
//  Coefficients: fit on the GENIE run3b signal sample; the deployed all-data fit
//  reproduces the cross-validated numbers (42.7% / 62.1% in-sample), so it is not
//  overfit. Refit with fit_final.py if the sample changes.
// ============================================================================
#include <cmath>
#include <vector>
#include <cstddef>

namespace PidAssign {

// feature order: sc, llr, lpr, lpi, lln, ldv, snh, lpi^2, lpr^2, llr^2, sc^2, lpi*lpr, llr*lpi
static const double PI_B = -0.065954;
static const double PI_C[13] = {-0.036953,  2.437298,  2.971687, -3.526364,  1.956626,
                                -0.275767, -0.112188,  0.610087, -0.917478, -3.641604,
                                -0.813744,  0.306769, -1.322442};
static const double PR_B = -0.851471;
static const double PR_C[13] = { 0.733462, -1.199681, -0.538603,  0.588946, -0.116288,
                                 0.033157, -0.003959,  0.358927, -0.363737,  0.780561,
                                 1.370665, -0.215181, -1.219626};

// Sentinel-safe log10 with the same clamps the fit used.
static inline double L10(double v, double lo){ return std::log10(v > lo ? v : lo); }

// Build the 13 features for one track. chi2 sentinels (<= -900) map to 300, the
// value the fit used for "no PID information".
static inline void features(double score, double llr, double chipr, double chipi,
                            double len, double dvtx, double nhitsY, double f[13])
{
  const double cpr = (chipr > -900.0 && std::isfinite(chipr)) ? chipr : 300.0;
  const double cpi = (chipi > -900.0 && std::isfinite(chipi)) ? chipi : 300.0;
  const double sc  = std::isfinite(score) ? score : 0.0;
  const double lr  = std::isfinite(llr)   ? llr   : 0.0;
  const double lpr = L10(cpr, 1e-2), lpi = L10(cpi, 1e-2);
  const double lln = L10(std::isfinite(len)?len:0.0, 1e-2);
  const double ldv = L10(std::isfinite(dvtx)?dvtx:0.0, 1e-2);
  const double snh = std::sqrt(nhitsY > 0 ? nhitsY : 0.0);
  f[0]=sc; f[1]=lr; f[2]=lpr; f[3]=lpi; f[4]=lln; f[5]=ldv; f[6]=snh;
  f[7]=lpi*lpi; f[8]=lpr*lpr; f[9]=lr*lr; f[10]=sc*sc; f[11]=lpi*lpr; f[12]=lr*lpi;
}
static inline double dot(const double* c, double b, const double f[13]){
  double s=b; for(int i=0;i<13;++i) s += c[i]*f[i]; return s;
}
static inline double pionScore(double score,double llr,double chipr,double chipi,
                               double len,double dvtx,double nhitsY){
  double f[13]; features(score,llr,chipr,chipi,len,dvtx,nhitsY,f); return dot(PI_C,PI_B,f);
}
static inline double protonScore(double score,double llr,double chipr,double chipi,
                                 double len,double dvtx,double nhitsY){
  double f[13]; features(score,llr,chipr,chipi,len,dvtx,nhitsY,f); return dot(PR_C,PR_B,f);
}

// Choose the ordered pair (proton, pion) from `cand` maximising Sp(a)+Spi(b).
// `cand` holds the indices of the non-muon track-like PFPs. Returns false (and
// leaves pr/pi untouched) when fewer than two candidates exist -- the same
// condition under which the old rule failed to find a pion.
template <class Vec>
inline bool assign(const std::vector<int>& cand,
                   const Vec& score, const Vec& llr, const Vec& chipr, const Vec& chipi,
                   const Vec& len, const std::vector<double>& dvtx,
                   const std::vector<double>& nhitsY,
                   int& pr, int& pi)
{
  const size_t n = cand.size();
  if(n < 2) return false;
  std::vector<double> sp(n), si(n);
  for(size_t k=0;k<n;++k){
    const int i = cand[k];
    sp[k] = protonScore(score[i],llr[i],chipr[i],chipi[i],len[i],dvtx[i],nhitsY[i]);
    si[k] = pionScore  (score[i],llr[i],chipr[i],chipi[i],len[i],dvtx[i],nhitsY[i]);
  }
  double best = -1e300; int ba=-1, bb=-1;
  for(size_t a=0;a<n;++a) for(size_t b=0;b<n;++b){
    if(a==b) continue;
    const double s = sp[a] + si[b];
    if(s > best){ best=s; ba=(int)a; bb=(int)b; }
  }
  if(ba<0) return false;
  pr = cand[ba]; pi = cand[bb];
  return true;
}

} // namespace PidAssign
#endif
