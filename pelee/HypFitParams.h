#ifndef HypFitParams_h
#define HypFitParams_h
// ============================================================================
//  Measured parameters for the Sigma0 global hypothesis fit (Phase 0, 2026-08-22).
//
//  All numbers measured on truth-matched PFPs in the NuWro HyperonProduction
//  sample (700 files, 1323 signal events: Sigma0, nubar_mu CC, Lambda -> p pi-,
//  muon > 0.15 GeV/c, reco vertex AND true vertex in FV).
//
//  Re-measure before use on a different sample or after a reco change.
// ============================================================================
namespace HypFit {

  // ---- mass constraints ---------------------------------------------------
  //  NOTE: centre the chi2 on the MEASURED means, not the PDG masses. The
  //  reconstruction is biased (+7.7 MeV on Lambda, -11.4 MeV on Sigma0) and the
  //  biases are a substantial fraction of the widths, so using the PDG values
  //  would give every hypothesis a built-in offset.
  //    Lambda bias is consistent with taking the pion momentum from range under
  //    the MUON hypothesis; Sigma0 bias with the photon energy under-measurement
  //    (reco 0.0792 vs true 0.0861 GeV).
  constexpr double M_LAMBDA_PDG  = 1.115683;
  constexpr double M_SIGMA0_PDG  = 1.192642;
  constexpr double M_LAMBDA_RECO = 1.1234;   // gaussian core, n=397
  constexpr double SIG_LAMBDA    = 0.0272;   // 27.2 MeV
  constexpr double M_SIGMA0_RECO = 1.1812;   // gaussian core, n=243
  constexpr double SIG_SIGMA0    = 0.0454;   // 45.4 MeV  (~1.7x wider, as expected)

  // ---- topology -----------------------------------------------------------
  //  p-pi distance of closest approach between the two track LINES.
  //  Median 0.35 cm, 16-84% = 0.08 - 1.36 cm. Strongly peaked, use a one-sided
  //  scale or a template; do NOT use the Lambda decay DISTANCE, which is
  //  uncorrelated with truth (r = -0.014) -- see sigma0-lambda-displacement-dead.
  constexpr double DCA_PPI_MEDIAN = 0.35;    // cm
  constexpr double DCA_PPI_SCALE  = 1.00;    // cm, ~84th percentile

  //  Muon is the longest of the three tracks 92.1% of the time -> soft prior
  //  worth ln(0.921/0.079) = 2.46 in log-likelihood. NOT a cut.
  constexpr double P_MU_LONGEST   = 0.921;

  // ---- photon -------------------------------------------------------------
  constexpr double EGAMMA_RECO_MEAN = 0.0792; // GeV
  constexpr double EGAMMA_TRUE_MEAN = 0.0861; // GeV; reco low by ~8%

  // ---- per-role PID templates: median [16%, 84%] --------------------------
  //  role order: 0 = mu+, 1 = proton, 2 = pi-, 3 = gamma
  //  n = 1213 / 808 / 650 / 814
  //
  //  var                 mu+                proton              pi-                 gamma
  //  trk_shr_score       1.00[0.92,1.00]    0.90[0.59,0.96]     0.80[0.38,0.93]     0.11[0.02,0.58]
  //  trk_llrpid          0.96[0.87,0.98]   -0.29[-0.76,0.12]    0.14[-0.13,0.62]    0.56[0.26,0.74]
  //  trk_pid_chipr     180.1[131,222]      16.4[3.96,118]      96.2[55.1,207]     191.7[140,291]
  //  three_plane_dedx    2.47[2.10,3.71]   10.28[6.01,24.5]     7.52[3.68,40.4]     3.30[2.05,6.53]
  //  rms_deflection      0.01[0.00,0.02]    0.01[0.00,0.06]     0.01[0.00,0.08]     0.08[0.02,0.19]
  //  shr_dedx_plane2     2.42[2.08,3.47]    8.83[4.09,17.1]     5.57[2.46,17.9]     3.06[1.18,6.03]
  //  moliere_avg         1.78[0.93,3.75]    6.18[2.14,19.8]     8.95[3.13,19.8]    13.26[7.62,24.3]
  //  trk_length        138.7[46.8,294]     12.82[3.38,39.0]    10.82[4.27,22.6]    14.15[7.71,25.4]
  //
  //  READING: the muon (length 139 vs ~12, LLR 0.96) and the photon (score 0.11,
  //  wiggliness 8x any track, Moliere 13.3) are strongly identified; the proton is
  //  strong via chipr 16 vs 96-192. The PION HAS NO UNIQUE SIGNATURE -- it lies
  //  between the others on every variable. It must be identified by EXCLUSION plus
  //  the Lambda mass constraint once the other three are pinned, which is exactly
  //  what the sequential selection cannot do (pion candidate is 34% pure there).

  // ---- enabling fact ------------------------------------------------------
  //  98.9% of PFPs carry BOTH a track fit and a shower fit, so a pi- classified
  //  as a shower still has every trk_* variable. The relaxed candidate pool needs
  //  no re-reconstruction, only a change of eligibility rule.
  constexpr double FRAC_PFP_BOTH_FITS = 0.989;

} // namespace HypFit
#endif
