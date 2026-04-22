// ================================================================
//  CrossSection.C
//
//  Total cross-section upper limit for Sigma0 quasi-elastic
//  production from muon antineutrinos in MicroBooNE, NuMI RHC beam
//
//       nu_mu_bar + Ar → mu+ + Sigma0 + X
//       Sigma0 → Lambda + gamma
//
//  Two methods are computed and compared:
//    1. Feldman-Cousins (purely statistical)
//    2. Bayesian with Gaussian systematic convolution (recommended)
//
//  Usage:
//    root -l -b -q CrossSection.C
// ================================================================

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

// ================================================================
//  Simple progress printer  (ANSI colours, elapsed time)
// ================================================================
static auto t_start_cs = std::chrono::steady_clock::now();

static void PrintStep(int step, const char* msg)
{
  auto   now     = std::chrono::steady_clock::now();
  double elapsed = std::chrono::duration<double>(now - t_start_cs).count();
  std::cout << "\033[1;36m[Step " << std::setw(2) << step << "]\033[0m "
            << std::left << std::setw(52) << msg
            << "  (+" << std::fixed << std::setprecision(2) << elapsed << " s)"
            << std::endl;
}

static void PrintDone(const char* extra = "")
{
  if (extra && extra[0])
    std::cout << "          \033[0;32m✔  " << extra << "\033[0m" << std::endl;
  else
    std::cout << "          \033[0;32m✔  done\033[0m" << std::endl;
}

#include "TMath.h"
#include "TFeldmanCousins.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TGraph.h"
#include "TAxis.h"

// ================================================================
//  Global nuisance parameters — accessed inside TF1 lambdas
// ================================================================
static double g_mu_nom;   // nominal signal expectation = eff * Phi * n_Ar
static double g_sigma_n;  // absolute systematic uncertainty on mu_nom

// ================================================================
//  Likelihood marginalised over Gaussian efficiency nuisance
//
//  P(N_obs | sigma_xsec) =
//    integral over mu of:
//      Poisson(N_obs | mu) * Gauss(mu | mu_nom(sigma_xsec), sigma_n)
//
//  where mu_0 = sigma_xsec * exposure, exposure = eff * Phi * n_Ar
//  and sigma_n = mu_0 * syst_frac  (absolute uncertainty on mu).
//
//  This integral has a closed form only for N_obs=0; the general case
//  is evaluated numerically for correctness with any observed count.
//
//  Integration limits: extend to mu_0 + 10*sigma_n on the high side
//  so that the Poisson peak at mu ~ N_obs is always captured even
//  when N_obs >> mu_0.
// ================================================================
double MarginalLikelihood(double sigma_xsec, double exposure, double syst_frac,
                          int n_obs = 0)
{
  double mu_0    = sigma_xsec * exposure;   // central signal expectation
  double sigma_n = mu_0 * syst_frac;        // absolute uncertainty on mu

  if (sigma_n < 1e-12) {
    // No systematic: pure Poisson
    return TMath::PoissonI(n_obs, mu_0);
  }

  // Integration range: cover both the Gaussian peak (at mu_0) and
  // the Poisson peak (near mu ~ n_obs), plus generous tails.
  const int    N_steps = 5000;
  double mu_lo = TMath::Max(0.0,
                   TMath::Min(mu_0, (double)n_obs) - 8.0 * sigma_n);
  double mu_hi = TMath::Max(mu_0, (double)n_obs)   + 8.0 * sigma_n;
  double dmu   = (mu_hi - mu_lo) / N_steps;
  double result = 0.0;

  for (int i = 0; i < N_steps; i++) {
    double mu    = mu_lo + (i + 0.5) * dmu;
    double pois  = TMath::PoissonI(n_obs, mu);             // Poisson(n_obs|mu)
    double gauss = TMath::Gaus(mu, mu_0, sigma_n, true);   // normalised Gauss
    result += pois * gauss * dmu;
  }
  return result;
}

// ================================================================
//  Compute 90% CL Bayesian upper limit on sigma_xsec
//  with flat prior on sigma_xsec >= 0.
//
//  n_obs is passed through to MarginalLikelihood so the posterior
//  is correct for any observed count, not just zero.
//
//  sigma_max is set adaptively: start from a broad range and warn
//  if the CL is reached only near the boundary (truncation).
// ================================================================
double BayesianUpperLimit(double exposure, double syst_frac,
                          int n_obs = 0, double CL = 0.90)
{
  // Upper scan boundary: generous fraction of the Poisson best-fit
  // cross section so the posterior is always fully covered.
  // For n_obs=0 this gives ~50/exposure; for n_obs>0 we extend further.
  const int    N = 100000;
  double       sigma_max = (50.0 + 20.0 * n_obs) / exposure;
  const double dsigma    = sigma_max / N;

  double norm     = 0.0;
  double cumul    = 0.0;
  double sigma_ul = sigma_max;

  std::vector<double> posterior(N);
  for (int i = 0; i < N; i++) {
    double sigma  = (i + 0.5) * dsigma;
    posterior[i]  = MarginalLikelihood(sigma, exposure, syst_frac, n_obs);
    norm         += posterior[i] * dsigma;
  }

  bool found = false;
  for (int i = 0; i < N; i++) {
    double sigma = (i + 0.5) * dsigma;
    cumul += posterior[i] * dsigma / norm;
    if (cumul >= CL) {
      sigma_ul = sigma;
      found = true;
      break;
    }
  }
  if (!found)
    std::cerr << "\033[1;33m  WARNING: Bayesian posterior truncated — "
                 "increase sigma_max (currently "
              << sigma_max << " cm^2)\033[0m" << std::endl;
  return sigma_ul;
}

// ================================================================
//  Compute the equal-tail Bayesian lower limit on sigma_xsec.
//
//  For a central interval at confidence level CL the lower edge is
//  the (1-CL)/2 quantile of the posterior.  When n_obs = 0 the
//  posterior has all its weight near zero and the lower limit is 0;
//  it only becomes meaningfully non-zero for n_obs >= 1.
// ================================================================
double BayesianLowerLimit(double exposure, double syst_frac,
                          int n_obs = 0, double CL = 0.90)
{
  if (n_obs == 0) return 0.0;

  const int    N        = 100000;
  double       sigma_max = (50.0 + 20.0 * n_obs) / exposure;
  const double dsigma   = sigma_max / N;
  const double alpha_lo = (1.0 - CL) / 2.0;  // lower-tail probability

  // Build normalisation
  double norm = 0.0;
  std::vector<double> posterior(N);
  for (int i = 0; i < N; i++) {
    double sigma  = (i + 0.5) * dsigma;
    posterior[i]  = MarginalLikelihood(sigma, exposure, syst_frac, n_obs);
    norm         += posterior[i] * dsigma;
  }

  // Find lower quantile
  double cumul   = 0.0;
  double sigma_ll = 0.0;
  for (int i = 0; i < N; i++) {
    cumul += posterior[i] * dsigma / norm;
    if (cumul >= alpha_lo) {
      sigma_ll = (i + 0.5) * dsigma;
      break;
    }
  }
  return sigma_ll;
}

// ================================================================
void CrossSection()
{
  t_start_cs = std::chrono::steady_clock::now();   // reset timer

  std::cout << "\n\033[1;35m";
  std::cout << "=============================================================\n";
  std::cout << "  Sigma0 QE cross-section upper limit  —  MicroBooNE NuMI RHC\n";
  std::cout << "=============================================================\033[0m\n\n";

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // --------------------------------------------------------------
  // 1. Physical constants
  // --------------------------------------------------------------
  PrintStep(1, "Loading physical constants ...");
  const double N_A        = 6.02214076e23;  // Avogadro [mol^-1]
  const double rho_LAr    = 1.3836;         // LAr density at 87 K, 1 atm [g/cm^3]
  const double A_Ar       = 39.948;         // Ar molar mass [g/mol]
  const double cm2_to_pb  = 1.0e36;         // 1 cm^2 = 10^36 pb

  PrintDone();

  // --------------------------------------------------------------
  // 2. Experimental inputs
  // --------------------------------------------------------------
  PrintStep(2, "Setting experimental inputs ...");
  const double PoT         = 4.90e20;   // Protons-on-Target (NuMI run 3b)
  const int    N_obs       = 1.0;         // Observed candidate events (data)
  const double N_bkg       = 0.0;       // Predicted background after all cuts
  const double N_pred_init = 17.86;      // MC predicted signal in FV before selection  ← REPLACE with actual MC
  const double N_pred_sel  = 0.3;//*1.5677;       // MC predicted signal after all cuts           ← REPLACE with actual MC
  const double syst_frac   = 0.30;      // Total fractional systematic uncertainty (flux⊕eff⊕xsec⊕det⊕POT)
  const double CL          = 0.90;      // Confidence level (primary)
  const double CL68        = 0.68;      // Confidence level (secondary)

  // Signal selection efficiency (ratio of selected to initial MC events)
  const double efficiency  = N_pred_sel / N_pred_init;
  PrintDone(Form("eff = %.3f%%,  syst = %.0f%%", efficiency*100, syst_frac*100));

  // --------------------------------------------------------------
  // 3. MicroBooNE Fiducial Volume
  // --------------------------------------------------------------
  PrintStep(3, "Computing fiducial volume ...");
  //  x: [3.00, 253.35] cm
  //  y: [-112.53, 114.47] cm
  //  z: [3.1, 1033.0] cm  minus dead region [675.1, 775.1] cm
  const double x_lo = 3.00,    x_hi = 253.35;
  const double y_lo = -112.53, y_hi = 114.47;
  const double z_lo = 3.10,    z_hi = 1033.00;
  const double z_dead_lo = 675.1, z_dead_hi = 775.1;

  const double dx   = x_hi - x_lo;                        // 250.35 cm
  const double dy   = y_hi - y_lo;                        // 227.00 cm
  const double dz   = (z_hi - z_lo) - (z_dead_hi - z_dead_lo);  // 929.90 cm
  const double V_FV = dx * dy * dz;                       // [cm^3]

  PrintDone(Form("V_FV = %.2f cm^3", V_FV));

  // --------------------------------------------------------------
  // 4. Number of target Ar nuclei in FV
  // --------------------------------------------------------------
  PrintStep(4, "Computing number of target Ar nuclei in FV ...");
  const double mass_FV = rho_LAr * V_FV;             // [g]
  const double n_Ar    = (mass_FV / A_Ar) * N_A;    // number of Ar nuclei in FV

  PrintDone(Form("n_Ar = %.4e Ar nuclei", n_Ar));

  // --------------------------------------------------------------
  // 5. Integrated muon antineutrino flux at MicroBooNE (NuMI RHC)
  //
  //    *** IMPORTANT: Replace flux_per_POT with the integral of
  //        your numubar flux histogram from the PPFX flux files,
  //        in units of [numubar / cm^2 / POT], evaluated over
  //        the full energy range of the simulation.
  //
  //        Typical source:
  //          MicroBooNE_NuMI_RHC_flux_unisim.root  (numubar histogram)
  //          integral = h_numubar->Integral("width") * bin_width
  //
  //        Representative value for off-axis NuMI RHC at ~100m
  //        off-axis (MicroBooNE position): ~1--5 x 10^-8 /cm^2/POT
  //        The placeholder below gives N_pred consistent with input.
  //
  //    flux_per_POT = N_pred_init / (sigma_MC_ref * n_Ar * PoT)
  //    where sigma_MC_ref is the GENIE reference cross-section.
  //    If you know sigma_MC_ref, you can derive flux_per_POT here.
  // --------------------------------------------------------------
  PrintStep(5, "Computing integrated flux and exposure ...");
  const double flux_per_POT = 3.52298e-10;// 3.0e-8;              // [numubar / cm^2 / POT]  <-- REPLACE
  const double Phi          = flux_per_POT * PoT;  // integrated flux [numubar / cm^2]

  // Effective exposure:  ε × Φ × n_Ar  [events per (cm^2 per Ar nucleus)]
  const double exposure = efficiency * Phi * n_Ar;  // [cm^-2]

  PrintDone(Form("Phi = %.4e /cm^2  |  exposure = %.4e", Phi, exposure));

  // --------------------------------------------------------------
  // 6. Purely statistical upper limit — Feldman-Cousins
  // --------------------------------------------------------------
  PrintStep(6, "Running Feldman-Cousins (stat only) ...");
  TFeldmanCousins fc(CL);
  fc.SetMuStep(0.001);
  double N_UL_FC      = fc.CalculateUpperLimit(N_obs, N_bkg);
  double N_LL_FC      = fc.CalculateLowerLimit(N_obs, N_bkg);
  double sigma_UL_FC  = N_UL_FC / exposure;
  double sigma_LL_FC  = N_LL_FC / exposure;

  TFeldmanCousins fc68(CL68);
  fc68.SetMuStep(0.001);
  double N_UL_FC_68     = fc68.CalculateUpperLimit(N_obs, N_bkg);
  double N_LL_FC_68     = fc68.CalculateLowerLimit(N_obs, N_bkg);
  double sigma_UL_FC_68 = N_UL_FC_68 / exposure;
  double sigma_LL_FC_68 = N_LL_FC_68 / exposure;
  PrintDone(Form("FC 90%%: [%.3f, %.3f]   FC 68%%: [%.3f, %.3f]",
                 N_LL_FC, N_UL_FC, N_LL_FC_68, N_UL_FC_68));

  // --------------------------------------------------------------
  // 7. Bayesian upper limit with marginalised Gaussian systematic
  // --------------------------------------------------------------
  PrintStep(7, "Running Bayesian integration (5000 steps × 100000 scan) ...");
  double sigma_UL_Bayes    = BayesianUpperLimit(exposure, syst_frac, N_obs, CL);
  double sigma_LL_Bayes    = BayesianLowerLimit(exposure, syst_frac, N_obs, CL);
  double N_UL_Bayes        = sigma_UL_Bayes * exposure;
  double N_LL_Bayes        = sigma_LL_Bayes * exposure;
  double sigma_UL_Bayes_68 = BayesianUpperLimit(exposure, syst_frac, N_obs, CL68);
  double sigma_LL_Bayes_68 = BayesianLowerLimit(exposure, syst_frac, N_obs, CL68);
  double N_UL_Bayes_68     = sigma_UL_Bayes_68 * exposure;
  double N_LL_Bayes_68     = sigma_LL_Bayes_68 * exposure;
  PrintDone(Form("Bayes 90%%: [%.3f, %.3f]   Bayes 68%%: [%.3f, %.3f]",
                 N_LL_Bayes, N_UL_Bayes, N_LL_Bayes_68, N_UL_Bayes_68));

  // --------------------------------------------------------------
  // 9. Reference cross-section implied by MC prediction
  // --------------------------------------------------------------
  PrintStep(8, "Computing MC reference cross-section ...");
  double sigma_MC_ref = N_pred_init / (Phi * n_Ar);  // [cm^2/Ar]

  PrintDone(Form("sigma_MC_ref = %.4e cm^2/Ar  (%.4f pb/Ar)", sigma_MC_ref, sigma_MC_ref*cm2_to_pb));

  // --------------------------------------------------------------
  // 10. Sensitivity: expected upper limit (Asimov, N_obs = 0)
  // --------------------------------------------------------------
  PrintStep(9, "Computing expected (Asimov) upper limit ...");
  // Asimov sensitivity: expected upper limit assuming N_obs = 0
  double N_exp_UL_FC  = fc.CalculateUpperLimit(0, N_bkg);
  double sigma_exp_UL = BayesianUpperLimit(exposure, syst_frac, 0, CL);  // N_obs=0 Asimov

  // Ratio of limits to MC reference
  double ratio_FC_UL       = sigma_UL_FC       / sigma_MC_ref;
  double ratio_FC_LL       = (N_obs > 0) ? sigma_LL_FC       / sigma_MC_ref : 0.0;
  double ratio_Bayes_UL    = sigma_UL_Bayes    / sigma_MC_ref;
  double ratio_Bayes_LL    = (N_obs > 0) ? sigma_LL_Bayes    / sigma_MC_ref : 0.0;
  double ratio_FC_UL_68    = sigma_UL_FC_68    / sigma_MC_ref;
  double ratio_FC_LL_68    = (N_obs > 0) ? sigma_LL_FC_68    / sigma_MC_ref : 0.0;
  double ratio_Bayes_UL_68 = sigma_UL_Bayes_68 / sigma_MC_ref;
  double ratio_Bayes_LL_68 = (N_obs > 0) ? sigma_LL_Bayes_68 / sigma_MC_ref : 0.0;
  // Keep backward-compatible aliases used in PrintDone below
  double ratio_FC    = ratio_FC_UL;
  double ratio_Bayes = ratio_Bayes_UL;
  double ratio_FC_68    = ratio_FC_UL_68;
  double ratio_Bayes_68 = ratio_Bayes_UL_68;
  PrintDone(Form("sigma_exp_UL(Bayes,90%%) = %.4e cm^2/Ar  |  UL/sigma_MC = %.2f (Bayes 90%%)  %.2f (Bayes 68%%)",
                 sigma_exp_UL, ratio_Bayes, ratio_Bayes_68));

  // --------------------------------------------------------------
  // 11. Uncertainty breakdown
  // --------------------------------------------------------------
  PrintStep(10, "Computing systematic uncertainty breakdown ...");
  //    (example breakdown — adjust to your analysis)
  double bk_flux    = 0.15;  // flux normalisation uncertainty
  double bk_eff     = 0.20;  // selection efficiency uncertainty  (== syst_eff above)
  double bk_xsec    = 0.10;  // cross-section model (nuclear effects)
  double bk_det     = 0.12;  // detector response (recombination, SCE)
  double bk_pot     = 0.02;  // POT counting uncertainty
  // Quadrature sum
  double syst_total   = TMath::Sqrt(bk_flux*bk_flux + bk_eff*bk_eff +
                                    bk_xsec*bk_xsec + bk_det*bk_det +
                                    bk_pot*bk_pot);
  PrintDone(Form("total syst (quadrature) = %.1f%%", syst_total*100));

  // --------------------------------------------------------------
  // 12. Print results
  // --------------------------------------------------------------
  PrintStep(11, "Printing results summary ...");
  std::cout << std::endl;

  std::cout << "\n--- Input parameters ---" << std::endl;
  std::cout << "  PoT                          : " << std::scientific << PoT           << std::endl;
  std::cout << "  Integrated numubar flux      : " << Phi                               << " /cm^2" << std::endl;
  std::cout << "  (flux per POT placeholder)   : " << flux_per_POT                      << " /cm^2/POT" << std::endl;

  std::cout << "\n--- Fiducial volume ---" << std::endl;
  std::cout << std::fixed;
  std::cout << "  dx (x)                       : " << dx   << " cm" << std::endl;
  std::cout << "  dy (y)                       : " << dy   << " cm" << std::endl;
  std::cout << "  dz (z, excl. dead region)    : " << dz   << " cm" << std::endl;
  std::cout << "  V_FV                         : " << V_FV << " cm^3" << std::endl;
  std::cout << "  Mass in FV                   : " << mass_FV / 1e3 << " kg" << std::endl;
  std::cout << "  N_Ar (FV)                    : " << std::scientific << n_Ar  << std::endl;

  std::cout << "\n--- Event counts ---" << std::endl;
  std::cout << std::fixed;
  std::cout << "  N_pred (initial, FV)         : " << N_pred_init  << std::endl;
  std::cout << "  N_pred (after selection)     : " << N_pred_sel   << std::endl;
  std::cout << "  Signal selection efficiency  : " << efficiency*100 << " %" << std::endl;
  std::cout << "  N_background (predicted)     : " << N_bkg        << std::endl;
  std::cout << "  N_observed                   : " << N_obs        << std::endl;
  std::cout << "  Systematic uncertainty       : " << syst_frac*100 << " %" << std::endl;

  std::cout << "\n--- Systematic uncertainty breakdown ---" << std::endl;
  std::cout << "  Flux normalisation           : " << bk_flux*100 << " %" << std::endl;
  std::cout << "  Selection efficiency         : " << bk_eff*100  << " %" << std::endl;
  std::cout << "  Cross-section model          : " << bk_xsec*100 << " %" << std::endl;
  std::cout << "  Detector response            : " << bk_det*100  << " %" << std::endl;
  std::cout << "  POT counting                 : " << bk_pot*100  << " %" << std::endl;
  std::cout << "  Total (in quadrature)        : " << syst_total*100 << " %" << std::endl;

  std::cout << "\n--- Exposure ---" << std::endl;
  std::cout << std::scientific;
  std::cout << "  eff * Phi * n_Ar             : " << exposure << " events / (cm^2/Ar)" << std::endl;
  std::cout << "  sigma_MC_ref                 : " << sigma_MC_ref         << " cm^2/Ar" << std::endl;
  std::cout << "  sigma_MC_ref                 : " << sigma_MC_ref*cm2_to_pb << " pb/Ar"  << std::endl;

  const bool two_sided = (N_obs > 0);
  std::cout << "\n--- Interval on number of signal events ---" << std::endl;
  std::cout << std::fixed;
  if (two_sided) {
    std::cout << "  N [FC,    90% CL, stat]      : [" << N_LL_FC       << ", " << N_UL_FC       << "]" << std::endl;
    std::cout << "  N [FC,    68% CL, stat]      : [" << N_LL_FC_68    << ", " << N_UL_FC_68    << "]" << std::endl;
    std::cout << "  N [Bayes, 90% CL, stat+syst] : [" << N_LL_Bayes    << ", " << N_UL_Bayes    << "]" << std::endl;
    std::cout << "  N [Bayes, 68% CL, stat+syst] : [" << N_LL_Bayes_68 << ", " << N_UL_Bayes_68 << "]" << std::endl;
  } else {
    std::cout << "  N_UL [Feldman-Cousins, 90% CL, stat]     : " << N_UL_FC       << std::endl;
    std::cout << "  N_UL [Feldman-Cousins, 68% CL, stat]     : " << N_UL_FC_68    << std::endl;
    std::cout << "  N_UL [Bayesian,        90% CL, stat+syst]: " << N_UL_Bayes    << std::endl;
    std::cout << "  N_UL [Bayesian,        68% CL, stat+syst]: " << N_UL_Bayes_68 << std::endl;
  }

  std::cout << "\n--- Cross-section interval [cm^2/Ar] ---" << std::endl;
  std::cout << std::scientific;
  if (two_sided) {
    std::cout << "  sigma [FC,    90% CL]  : [" << sigma_LL_FC       << ", " << sigma_UL_FC       << "]" << std::endl;
    std::cout << "  sigma [FC,    68% CL]  : [" << sigma_LL_FC_68    << ", " << sigma_UL_FC_68    << "]" << std::endl;
    std::cout << "  sigma [Bayes, 90% CL]  : [" << sigma_LL_Bayes    << ", " << sigma_UL_Bayes    << "]" << std::endl;
    std::cout << "  sigma [Bayes, 68% CL]  : [" << sigma_LL_Bayes_68 << ", " << sigma_UL_Bayes_68 << "]" << std::endl;
    std::cout << std::endl;
    std::cout << "\n--- Cross-section interval [pb/Ar] ---" << std::endl;
    std::cout << "  sigma [FC,    90% CL]  : [" << sigma_LL_FC*cm2_to_pb       << ", " << sigma_UL_FC*cm2_to_pb       << "]" << std::endl;
    std::cout << "  sigma [FC,    68% CL]  : [" << sigma_LL_FC_68*cm2_to_pb    << ", " << sigma_UL_FC_68*cm2_to_pb    << "]" << std::endl;
    std::cout << "  sigma [Bayes, 90% CL]  : [" << sigma_LL_Bayes*cm2_to_pb    << ", " << sigma_UL_Bayes*cm2_to_pb    << "]" << std::endl;
    std::cout << "  sigma [Bayes, 68% CL]  : [" << sigma_LL_Bayes_68*cm2_to_pb << ", " << sigma_UL_Bayes_68*cm2_to_pb << "]" << std::endl;
  } else {
    std::cout << "  sigma_UL [FC,    90% CL]  : " << sigma_UL_FC       << " cm^2/Ar" << std::endl;
    std::cout << "  sigma_UL [FC,    68% CL]  : " << sigma_UL_FC_68    << " cm^2/Ar" << std::endl;
    std::cout << "  sigma_UL [Bayes, 90% CL]  : " << sigma_UL_Bayes    << " cm^2/Ar" << std::endl;
    std::cout << "  sigma_UL [Bayes, 68% CL]  : " << sigma_UL_Bayes_68 << " cm^2/Ar" << std::endl;
    std::cout << std::endl;
    std::cout << "  sigma_UL [FC,    90% CL]  : " << sigma_UL_FC*cm2_to_pb       << " pb/Ar" << std::endl;
    std::cout << "  sigma_UL [FC,    68% CL]  : " << sigma_UL_FC_68*cm2_to_pb    << " pb/Ar" << std::endl;
    std::cout << "  sigma_UL [Bayes, 90% CL]  : " << sigma_UL_Bayes*cm2_to_pb    << " pb/Ar" << std::endl;
    std::cout << "  sigma_UL [Bayes, 68% CL]  : " << sigma_UL_Bayes_68*cm2_to_pb << " pb/Ar" << std::endl;
  }

  std::cout << "\n--- Limit relative to MC reference cross-section ---" << std::endl;
  std::cout << std::fixed;
  if (two_sided) {
    std::cout << "  sigma/sigma_MC [FC,    90% CL]  : [" << ratio_FC_LL       << ", " << ratio_FC_UL       << "]" << std::endl;
    std::cout << "  sigma/sigma_MC [FC,    68% CL]  : [" << ratio_FC_LL_68    << ", " << ratio_FC_UL_68    << "]" << std::endl;
    std::cout << "  sigma/sigma_MC [Bayes, 90% CL]  : [" << ratio_Bayes_LL    << ", " << ratio_Bayes_UL    << "]" << std::endl;
    std::cout << "  sigma/sigma_MC [Bayes, 68% CL]  : [" << ratio_Bayes_LL_68 << ", " << ratio_Bayes_UL_68 << "]" << std::endl;
  } else {
    std::cout << "  sigma_UL/sigma_MC [FC,    90% CL]   : " << ratio_FC       << std::endl;
    std::cout << "  sigma_UL/sigma_MC [FC,    68% CL]   : " << ratio_FC_68    << std::endl;
    std::cout << "  sigma_UL/sigma_MC [Bayes, 90% CL]   : " << ratio_Bayes    << std::endl;
    std::cout << "  sigma_UL/sigma_MC [Bayes, 68% CL]   : " << ratio_Bayes_68 << std::endl;
  }

  std::cout << "\n=============================================================\n" << std::endl;

  // --------------------------------------------------------------
  // 13. Plot the Bayesian posterior PDF and the upper limit
  // --------------------------------------------------------------
  PrintStep(12, "Building posterior PDF plot ...");
  const int Npts = 2000;
  double sigma_max_plot = 4.0 * sigma_UL_Bayes;
  double dsigma         = sigma_max_plot / Npts;

  // Compute posterior (unnormalised)
  std::vector<double> vSigma(Npts), vPost(Npts);
  double norm_post = 0.0;
  for (int i = 0; i < Npts; i++) {
    vSigma[i] = (i + 0.5) * dsigma;
    vPost[i]  = MarginalLikelihood(vSigma[i], exposure, syst_frac, N_obs);
    norm_post += vPost[i] * dsigma;
  }
  for (int i = 0; i < Npts; i++) vPost[i] /= norm_post;

  TGraph *gPost = new TGraph(Npts, vSigma.data(), vPost.data());
  gPost->SetLineColor(kBlue+1);
  gPost->SetLineWidth(3);
  gPost->SetFillColor(kBlue-9);
  gPost->SetFillStyle(1001);

  // Filled area: [0, UL_90] when upper-limit only; [LL_90, UL_90] when two-sided
  std::vector<double> vSigmaFill, vPostFill;
  double fill90_lo = two_sided ? sigma_LL_Bayes    : 0.0;
  double fill68_lo = two_sided ? sigma_LL_Bayes_68 : 0.0;
  vSigmaFill.push_back(fill90_lo); vPostFill.push_back(0.0);
  for (int i = 0; i < Npts; i++) {
    if (vSigma[i] < fill90_lo) continue;
    if (vSigma[i] > sigma_UL_Bayes) break;
    vSigmaFill.push_back(vSigma[i]);
    vPostFill.push_back(vPost[i]);
  }
  vSigmaFill.push_back(sigma_UL_Bayes); vPostFill.push_back(0.0);
  TGraph *gFill = new TGraph(vSigmaFill.size(), vSigmaFill.data(), vPostFill.data());
  gFill->SetFillColor(kBlue-9);
  gFill->SetFillStyle(1001);
  gFill->SetLineColor(kBlue-9);

  // Display unit: 10^{-40} cm^2/Ar
  // x_display = sigma_Ar * 1e40   (value in units of 1e-40 cm^2/Ar)
  // y_display = PDF / 1e40        (PDF per unit of 1e-40 cm^2/Ar)
  const double unit_scale = 1.0e40;

  TCanvas *c = new TCanvas("c_posterior", "Bayesian Posterior", 900, 650);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.13);

  double x_axis_max = sigma_max_plot * unit_scale;
  TH1D *hFrame = new TH1D("hFrame", "", 100, 0, x_axis_max);
  double y_max = *std::max_element(vPost.begin(), vPost.end()) / unit_scale * 1.25;
  hFrame->SetMaximum(y_max);
  hFrame->SetMinimum(0.0);
  hFrame->GetXaxis()->SetTitle("#sigma (#bar{#nu}_{#mu} Ar #rightarrow #mu^{+} + #Sigma^{0})Br(#Lambda #rightarrow p + #pi^{+})  [10^{-40} cm^{2} / Ar]");
  hFrame->GetYaxis()->SetTitle("Posterior probability density  [10^{40} / cm^{2} Ar^{-1}]");
  hFrame->GetXaxis()->SetTitleSize(0.042);
  hFrame->GetYaxis()->SetTitleSize(0.042);
  hFrame->GetXaxis()->SetLabelSize(0.040);
  hFrame->GetYaxis()->SetLabelSize(0.040);
  hFrame->GetXaxis()->SetTitleOffset(1.2);
  hFrame->Draw("AXIS");

  // Build display graphs: x → 10^{-40} cm^2/Ar,  y → PDF per (10^{-40} cm^2/Ar)
  TGraph *gPost_Ar = new TGraph(Npts);
  for (int i = 0; i < Npts; i++)
    gPost_Ar->SetPoint(i, vSigma[i] * unit_scale,
                          vPost[i]  / unit_scale);
  gPost_Ar->SetLineColor(kBlue+1);
  gPost_Ar->SetLineWidth(3);

  TGraph *gFill_Ar = new TGraph(vSigmaFill.size());
  for (size_t i = 0; i < vSigmaFill.size(); i++)
    gFill_Ar->SetPoint(i, vSigmaFill[i] * unit_scale,
                          vPostFill[i]  / unit_scale);
  gFill_Ar->SetFillColor(kAzure-9);
  gFill_Ar->SetFillStyle(1001);
  gFill_Ar->SetLineColor(kAzure-9);

  gFill_Ar->Draw("F SAME");
  gPost_Ar->Draw("L SAME");

  // 68% CL filled region on posterior
  std::vector<double> vSigmaFill68, vPostFill68;
  vSigmaFill68.push_back(fill68_lo); vPostFill68.push_back(0.0);
  for (int i = 0; i < Npts; i++) {
    if (vSigma[i] < fill68_lo) continue;
    if (vSigma[i] > sigma_UL_Bayes_68) break;
    vSigmaFill68.push_back(vSigma[i]);
    vPostFill68.push_back(vPost[i]);
  }
  vSigmaFill68.push_back(sigma_UL_Bayes_68); vPostFill68.push_back(0.0);
  TGraph *gFill68_Ar = new TGraph(vSigmaFill68.size());
  for (size_t i = 0; i < vSigmaFill68.size(); i++)
    gFill68_Ar->SetPoint(i, vSigmaFill68[i] * unit_scale,
                            vPostFill68[i]  / unit_scale);
  gFill68_Ar->SetFillColor(kAzure-4);
  gFill68_Ar->SetFillStyle(1001);
  gFill68_Ar->SetLineColor(kAzure-4);
  gFill68_Ar->Draw("F SAME");

  // Redraw posterior curve on top so it isn't obscured
  gPost_Ar->Draw("L SAME");

  // Upper limit lines: x in units of 10^{-40} cm^2/Ar
  TLine *lUL_Bayes    = new TLine(sigma_UL_Bayes    * unit_scale, 0, sigma_UL_Bayes    * unit_scale, y_max*0.85);
  TLine *lUL_Bayes_68 = new TLine(sigma_UL_Bayes_68 * unit_scale, 0, sigma_UL_Bayes_68 * unit_scale, y_max*0.85);
  TLine *lUL_FC       = new TLine(sigma_UL_FC        * unit_scale, 0, sigma_UL_FC        * unit_scale, y_max*0.85);
  TLine *lUL_FC_68    = new TLine(sigma_UL_FC_68     * unit_scale, 0, sigma_UL_FC_68     * unit_scale, y_max*0.85);

  lUL_Bayes   ->SetLineColor(kBlue+1); lUL_Bayes   ->SetLineWidth(3); lUL_Bayes   ->SetLineStyle(1);
  lUL_Bayes_68->SetLineColor(kBlue+1); lUL_Bayes_68->SetLineWidth(2); lUL_Bayes_68->SetLineStyle(2);
  lUL_FC      ->SetLineColor(kRed+1);  lUL_FC      ->SetLineWidth(3); lUL_FC      ->SetLineStyle(1);
  lUL_FC_68   ->SetLineColor(kRed+1);  lUL_FC_68   ->SetLineWidth(2); lUL_FC_68   ->SetLineStyle(2);

  lUL_Bayes->Draw(); lUL_Bayes_68->Draw(); lUL_FC->Draw(); lUL_FC_68->Draw();

  // Lower limit lines (only drawn when two-sided)
  TLine *lLL_Bayes    = nullptr;
  TLine *lLL_Bayes_68 = nullptr;
  TLine *lLL_FC       = nullptr;
  TLine *lLL_FC_68    = nullptr;
  if (two_sided) {
    lLL_Bayes    = new TLine(sigma_LL_Bayes    * unit_scale, 0, sigma_LL_Bayes    * unit_scale, y_max*0.85);
    lLL_Bayes_68 = new TLine(sigma_LL_Bayes_68 * unit_scale, 0, sigma_LL_Bayes_68 * unit_scale, y_max*0.85);
    lLL_FC       = new TLine(sigma_LL_FC        * unit_scale, 0, sigma_LL_FC        * unit_scale, y_max*0.85);
    lLL_FC_68    = new TLine(sigma_LL_FC_68     * unit_scale, 0, sigma_LL_FC_68     * unit_scale, y_max*0.85);
    lLL_Bayes   ->SetLineColor(kBlue+1); lLL_Bayes   ->SetLineWidth(3); lLL_Bayes   ->SetLineStyle(1);
    lLL_Bayes_68->SetLineColor(kBlue+1); lLL_Bayes_68->SetLineWidth(2); lLL_Bayes_68->SetLineStyle(2);
    lLL_FC      ->SetLineColor(kRed+1);  lLL_FC      ->SetLineWidth(3); lLL_FC      ->SetLineStyle(1);
    lLL_FC_68   ->SetLineColor(kRed+1);  lLL_FC_68   ->SetLineWidth(2); lLL_FC_68   ->SetLineStyle(2);
    lLL_Bayes->Draw(); lLL_Bayes_68->Draw(); lLL_FC->Draw(); lLL_FC_68->Draw();
  }

  // MC reference line  (x in units of 10^{-40} cm^2/Ar)
  TLine *lMC = new TLine(sigma_MC_ref * unit_scale, 0, sigma_MC_ref * unit_scale, y_max*0.85);
  lMC->SetLineColor(kOrange+1); lMC->SetLineWidth(3); lMC->SetLineStyle(5);
  lMC->Draw();

  // Legend — values quoted in units of 10^{-40} cm^2/Ar
  TLegend *leg = new TLegend(0.38, 0.20, 0.88, 0.62);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->AddEntry(gPost_Ar,   "Bayesian posterior (flat prior)", "L");
  leg->AddEntry(gFill_Ar,   "90% CL region (Bayesian)", "F");
  leg->AddEntry(gFill68_Ar, "68% CL region (Bayesian)", "F");
  if (two_sided) {
    leg->AddEntry(lUL_Bayes,
      Form("#sigma^{Bayes} = [%.2f, %.2f] #times 10^{-40} cm^{2}/Ar (90%% CL)",
           sigma_LL_Bayes * unit_scale, sigma_UL_Bayes * unit_scale), "L");
    leg->AddEntry(lUL_Bayes_68,
      Form("#sigma^{Bayes} = [%.2f, %.2f] #times 10^{-40} cm^{2}/Ar (68%% CL)",
           sigma_LL_Bayes_68 * unit_scale, sigma_UL_Bayes_68 * unit_scale), "L");
    leg->AddEntry(lUL_FC,
      Form("#sigma^{FC}    = [%.2f, %.2f] #times 10^{-40} cm^{2}/Ar (90%% CL, stat only)",
           sigma_LL_FC * unit_scale, sigma_UL_FC * unit_scale), "L");
    leg->AddEntry(lUL_FC_68,
      Form("#sigma^{FC}    = [%.2f, %.2f] #times 10^{-40} cm^{2}/Ar (68%% CL, stat only)",
           sigma_LL_FC_68 * unit_scale, sigma_UL_FC_68 * unit_scale), "L");
  } else {
    leg->AddEntry(lUL_Bayes,
      Form("#sigma_{UL}^{Bayes} = %.2f #times 10^{-40} cm^{2}/Ar (90%% CL)",
           sigma_UL_Bayes * unit_scale), "L");
    leg->AddEntry(lUL_Bayes_68,
      Form("#sigma_{UL}^{Bayes} = %.2f #times 10^{-40} cm^{2}/Ar (68%% CL)",
           sigma_UL_Bayes_68 * unit_scale), "L");
    leg->AddEntry(lUL_FC,
      Form("#sigma_{UL}^{FC}    = %.2f #times 10^{-40} cm^{2}/Ar (90%% CL, stat only)",
           sigma_UL_FC * unit_scale), "L");
    leg->AddEntry(lUL_FC_68,
      Form("#sigma_{UL}^{FC}    = %.2f #times 10^{-40} cm^{2}/Ar (68%% CL, stat only)",
           sigma_UL_FC_68 * unit_scale), "L");
  }
  leg->AddEntry(lMC,
    Form("#sigma_{MC ref}     = %.2f #times 10^{-40} cm^{2}/Ar",
         sigma_MC_ref * unit_scale), "L");
  leg->Draw();

  // Annotation
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.036);
  latex.SetTextColor(kGray+2);
  latex.DrawLatex(0.45, 0.86, "MicroBooNE NuMI RHC,  4.9 #times 10^{20} POT");
  latex.DrawLatex(0.45, 0.81, Form("#bar{#nu}_{#mu} Ar #rightarrow #mu^{+}+#Sigma^{0}, "
                                    "(#Sigma^{0} #rightarrow #gamma + #Lambda, #Lambda #rightarrow p + #pi^{+}) "));
  latex.DrawLatex(0.45, 0.76, Form("N_{obs}=%d,  N_{bkg}=%.0f,  #varepsilon=%.1f%%",
                                      N_obs, N_bkg, efficiency*100));
  latex.DrawLatex(0.45, 0.71, Form("Syst. uncertainty: %.0f%%", syst_frac*100));

  c->Print("Sigma0_CrossSection_UpperLimit.png");
  c->Print("Sigma0_CrossSection_UpperLimit.pdf");

  auto t_end = std::chrono::steady_clock::now();
  double total = std::chrono::duration<double>(t_end - t_start_cs).count();
  PrintDone(Form("Sigma0_CrossSection_UpperLimit.{png,pdf}  saved"));
  std::cout << "\n\033[1;35m  All done in " << std::fixed << std::setprecision(2)
            << total << " s\033[0m\n" << std::endl;
}
