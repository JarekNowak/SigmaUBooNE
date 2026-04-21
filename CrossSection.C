// ================================================================
//  CrossSection.C
//
//  Total cross-section upper limit for Sigma0 quasi-elastic
//  production from muon antineutrinos in MicroBooNE, NuMI RHC beam
//
//       nu_mu_bar + Ar → mu+ + Sigma0 + X
//       Sigma0 → Lambda + gamma
//
//  Three methods are computed and compared:
//    1. Feldman-Cousins (purely statistical)
//    2. Bayesian with Gaussian systematic convolution (recommended)
//    3. TRolke frequentist profile likelihood with systematics
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
#include "TRolke.h"
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
static double g_mu_nom;   // nominal signal expectation = eff * Phi * N_T
static double g_sigma_n;  // absolute systematic uncertainty on mu_nom

// ================================================================
//  Likelihood marginalised over Gaussian efficiency nuisance
//
//  P(N_obs | sigma_xsec) =
//    integral over mu of:
//      Poisson(N_obs | mu) * Gauss(mu | mu_nom(sigma_xsec), sigma_n)
//
//  where mu_0 = sigma_xsec * exposure, exposure = eff * Phi * N_T
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
  const double A_nuc      = 40;             // nucleons per ^40Ar nucleus
  const double cm2_to_pb  = 1.0e36;         // 1 cm^2 = 10^36 pb

  PrintDone();

  // --------------------------------------------------------------
  // 2. Experimental inputs
  // --------------------------------------------------------------
  PrintStep(2, "Setting experimental inputs ...");
  const double PoT         = 4.90e20;   // Protons-on-Target (NuMI run 3b)
  const int    N_obs       = 0.0;         // Observed candidate events (data)
  const double N_bkg       = 0.0;       // Predicted background after all cuts
  const double N_pred_init = 17.86;      // MC predicted signal in FV before selection  ← REPLACE with actual MC
  const double N_pred_sel  = 0.3;       // MC predicted signal after all cuts           ← REPLACE with actual MC
  const double syst_frac   = 0.30;      // Total fractional systematic uncertainty (flux⊕eff⊕xsec⊕det⊕POT)
  const double CL          = 0.90;      // Confidence level

  // Signal selection efficiency (ratio of selected to initial MC events)
  const double efficiency  = N_pred_sel / N_pred_init;
  // TRolke's Gaussian efficiency model accounts only for the
  // efficiency-specific uncertainty; the other systematics are folded
  // into the Bayesian integral separately via syst_frac.
  const double syst_eff    = 0.20;     // efficiency systematic only
  const double eff_err     = efficiency * syst_eff;   // absolute uncertainty

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
  // 4. Number of target nucleons in FV
  // --------------------------------------------------------------
  PrintStep(4, "Computing number of target nucleons ...");
  const double mass_FV = rho_LAr * V_FV;             // [g]
  const double n_Ar    = (mass_FV / A_Ar) * N_A;    // number of Ar nuclei
  const double N_T     = n_Ar * A_nuc;               // total nucleons

  PrintDone(Form("N_T = %.4e nucleons", N_T));

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
  //    flux_per_POT = N_pred_init / (sigma_MC_ref * N_T * PoT)
  //    where sigma_MC_ref is the GENIE reference cross-section.
  //    If you know sigma_MC_ref, you can derive flux_per_POT here.
  // --------------------------------------------------------------
  PrintStep(5, "Computing integrated flux and exposure ...");
    const double flux_per_POT = 3.52298e-10;// 3.0e-8;              // [numubar / cm^2 / POT]  <-- REPLACE
  const double Phi          = flux_per_POT * PoT;  // integrated flux [numubar / cm^2]

  // Effective exposure:  ε × Φ × N_T  [dimensionless / cm^2 × cm^2 = events/cm^2]
  // (units: [numubar/cm^2] × [nucleons] × [dimensionless] → events per (cm^2 per nucleon))
  const double exposure = efficiency * Phi * N_T;  // [cm^-2]

  PrintDone(Form("Phi = %.4e /cm^2  |  exposure = %.4e", Phi, exposure));

  // --------------------------------------------------------------
  // 6. Purely statistical upper limit — Feldman-Cousins
  // --------------------------------------------------------------
  PrintStep(6, "Running Feldman-Cousins (stat only) ...");
  TFeldmanCousins fc(CL);
  fc.SetMuStep(0.001);
  double N_UL_FC    = fc.CalculateUpperLimit(N_obs, N_bkg);
  double sigma_UL_FC = N_UL_FC / exposure;
  PrintDone(Form("N_UL(FC) = %.3f  →  sigma_UL = %.4e cm^2/nucleon", N_UL_FC, sigma_UL_FC));

  // --------------------------------------------------------------
  // 7. Statistical + systematic upper limit — TRolke
  //    Model: Poisson signal, known background (= 0),
  //           Gaussian efficiency with 20% relative uncertainty.
  //
  //  SetBounding(true) is REQUIRED when efficiency is small or
  //  N_obs is large.  Without it, the internal Brent root-finder
  //  searches an unbounded interval; with mu_hat ~ N_obs/eff >> 1
  //  it doubles its range indefinitely and never returns.
  // --------------------------------------------------------------
  PrintStep(7, "Running TRolke (stat + syst, Gaussian efficiency) ...");
  TRolke rolke;
  rolke.SetCL(CL);
  rolke.SetBounding(true);   // CRITICAL: prevents infinite root search
  // SetKnownBkgGaussEff(x, bk, em, sde):
  //   x   = observed events
  //   bk  = known background (exactly 0)
  //   em  = efficiency estimate
  //   sde = std dev of efficiency (absolute, eff-only syst)
  rolke.SetKnownBkgGaussEff(N_obs, N_bkg, efficiency, eff_err);
  double N_UL_Rolke     = rolke.GetUpperLimit();
  double sigma_UL_Rolke = N_UL_Rolke / (Phi * N_T);
  PrintDone(Form("N_UL(Rolke) = %.3f  →  sigma_UL = %.4e cm^2/nucleon",
                 N_UL_Rolke, sigma_UL_Rolke));
    
  // --------------------------------------------------------------
  // 8. Bayesian upper limit with marginalised Gaussian systematic
  // --------------------------------------------------------------
  PrintStep(8, "Running Bayesian integration (5000 steps × 100000 scan) ...");
  double sigma_UL_Bayes = BayesianUpperLimit(exposure, syst_frac, N_obs, CL);
  double N_UL_Bayes     = sigma_UL_Bayes * exposure;
  PrintDone(Form("N_UL(Bayes) = %.3f  →  sigma_UL = %.4e cm^2/nucleon", N_UL_Bayes, sigma_UL_Bayes));

  // --------------------------------------------------------------
  // 9. Reference cross-section implied by MC prediction
  // --------------------------------------------------------------
  PrintStep(9, "Computing MC reference cross-section ...");
  double sigma_MC_ref = N_pred_init / (Phi * N_T);

  PrintDone(Form("sigma_MC_ref = %.4e cm^2/nucleon  (%.4f pb)", sigma_MC_ref, sigma_MC_ref*cm2_to_pb));

  // --------------------------------------------------------------
  // 10. Sensitivity: expected upper limit (Asimov, N_obs = 0)
  // --------------------------------------------------------------
  PrintStep(10, "Computing expected (Asimov) upper limit ...");
  // Asimov sensitivity: expected upper limit assuming N_obs = 0
  double N_exp_UL_FC  = fc.CalculateUpperLimit(0, N_bkg);
  double sigma_exp_UL = BayesianUpperLimit(exposure, syst_frac, 0, CL);  // N_obs=0 Asimov

  // Ratio of upper limit to MC reference
  double ratio_FC    = sigma_UL_FC    / sigma_MC_ref;
  double ratio_Rolke = sigma_UL_Rolke / sigma_MC_ref;
  double ratio_Bayes = sigma_UL_Bayes / sigma_MC_ref;
  PrintDone(Form("sigma_exp_UL = %.4e cm^2  |  UL/sigma_MC = %.2f (Bayes)", sigma_exp_UL, ratio_Bayes));

  // --------------------------------------------------------------
  // 11. Uncertainty breakdown
  // --------------------------------------------------------------
  PrintStep(11, "Computing systematic uncertainty breakdown ...");
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
  PrintStep(12, "Printing results summary ...");
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
  std::cout << "  N_nucleons (FV)              : " << N_T  << std::endl;

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
  std::cout << "  eff * Phi * N_T              : " << exposure << " events / (cm^2/nucleon)" << std::endl;
  std::cout << "  sigma_MC_ref                 : " << sigma_MC_ref         << " cm^2/nucleon" << std::endl;
  std::cout << "  sigma_MC_ref                 : " << sigma_MC_ref*cm2_to_pb << " pb/nucleon"  << std::endl;

  std::cout << "\n--- Upper limits on number of signal events (" << CL*100 << "% CL) ---" << std::endl;
  std::cout << std::fixed;
  std::cout << "  N_UL  [Feldman-Cousins, stat]  : " << N_UL_FC    << std::endl;
  std::cout << "  N_UL  [TRolke, stat+syst]      : " << N_UL_Rolke << std::endl;
  std::cout << "  N_UL  [Bayesian, stat+syst]    : " << N_UL_Bayes << std::endl;

  std::cout << "\n--- Cross-section upper limits (" << CL*100 << "% CL) ---" << std::endl;
  std::cout << std::scientific;
  std::cout << "  sigma_UL [Feldman-Cousins]   : " << sigma_UL_FC    << " cm^2/nucleon"
            << "   (" << sigma_UL_FC   *A_nuc << " cm^2/Ar)" << std::endl;
  std::cout << "  sigma_UL [TRolke]            : " << sigma_UL_Rolke << " cm^2/nucleon"
            << "   (" << sigma_UL_Rolke*A_nuc << " cm^2/Ar)" << std::endl;
  std::cout << "  sigma_UL [Bayesian]          : " << sigma_UL_Bayes << " cm^2/nucleon"
            << "   (" << sigma_UL_Bayes*A_nuc << " cm^2/Ar)" << std::endl;
  std::cout << std::endl;
  std::cout << "  sigma_UL [Feldman-Cousins]   : " << sigma_UL_FC*cm2_to_pb    << " pb/nucleon" << std::endl;
  std::cout << "  sigma_UL [TRolke]            : " << sigma_UL_Rolke*cm2_to_pb << " pb/nucleon" << std::endl;
  std::cout << "  sigma_UL [Bayesian]          : " << sigma_UL_Bayes*cm2_to_pb << " pb/nucleon" << std::endl;

  std::cout << "\n--- Upper limit relative to MC reference cross-section ---" << std::endl;
  std::cout << std::fixed;
  std::cout << "  sigma_UL/sigma_MC [FC]       : " << ratio_FC    << std::endl;
  std::cout << "  sigma_UL/sigma_MC [TRolke]   : " << ratio_Rolke << std::endl;
  std::cout << "  sigma_UL/sigma_MC [Bayesian] : " << ratio_Bayes << std::endl;

  std::cout << "\n=============================================================\n" << std::endl;

  // --------------------------------------------------------------
  // 13. Plot the Bayesian posterior PDF and the upper limit
  // --------------------------------------------------------------
  PrintStep(13, "Building posterior PDF plot ...");
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

  // Filled area up to upper limit (90% region)
  std::vector<double> vSigmaFill, vPostFill;
  vSigmaFill.push_back(0.0); vPostFill.push_back(0.0);
  for (int i = 0; i < Npts; i++) {
    if (vSigma[i] > sigma_UL_Bayes) break;
    vSigmaFill.push_back(vSigma[i]);
    vPostFill.push_back(vPost[i]);
  }
  vSigmaFill.push_back(sigma_UL_Bayes); vPostFill.push_back(0.0);
  TGraph *gFill = new TGraph(vSigmaFill.size(), vSigmaFill.data(), vPostFill.data());
  gFill->SetFillColor(kBlue-9);
  gFill->SetFillStyle(1001);
  gFill->SetLineColor(kBlue-9);

  // Display unit: cm^2 per Ar nucleus  (= cm^2/nucleon * A_nuc)
  // The posterior PDF transforms as P_Ar(s) = P_nuc(s/A_nuc) / A_nuc,
  // so x-values are multiplied by A_nuc and y-values divided by A_nuc,
  // keeping the integral normalised to 1.
  const double nuc_to_Ar = static_cast<double>(A_nuc);  // = 40

  TCanvas *c = new TCanvas("c_posterior", "Bayesian Posterior", 900, 650);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.13);

  // Display unit: 10^{-40} cm^2/Ar
  // x_display = x_cm2_per_Ar * 1e40   (value in units of 1e-40 cm^2/Ar)
  // y_display = y_per_cm2_per_Ar / 1e40  (PDF per unit of 1e-40 cm^2/Ar)
  const double unit_scale = 1.0e40;  // multiply x by this, divide y by this

  // Axis frame — x in units of 10^{-40} cm^2/Ar
  double x_axis_max = sigma_max_plot * nuc_to_Ar * unit_scale;
  TH1D *hFrame = new TH1D("hFrame", "", 100, 0, x_axis_max);
  double y_max = *std::max_element(vPost.begin(), vPost.end()) / nuc_to_Ar / unit_scale * 1.25;
  hFrame->SetMaximum(y_max);
  hFrame->SetMinimum(0.0);
  hFrame->GetXaxis()->SetTitle("#sigma (#bar{#nu}_{#mu} Ar #rightarrow #mu^{+} #Sigma^{0} X)  [10^{-40} cm^{2} / Ar]");
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
    gPost_Ar->SetPoint(i, vSigma[i] * nuc_to_Ar * unit_scale,
                          vPost[i]  / nuc_to_Ar  / unit_scale);
  gPost_Ar->SetLineColor(kBlue+1);
  gPost_Ar->SetLineWidth(3);

  TGraph *gFill_Ar = new TGraph(vSigmaFill.size());
  for (size_t i = 0; i < vSigmaFill.size(); i++)
    gFill_Ar->SetPoint(i, vSigmaFill[i] * nuc_to_Ar * unit_scale,
                          vPostFill[i]  / nuc_to_Ar  / unit_scale);
  gFill_Ar->SetFillColor(kAzure-9);
  gFill_Ar->SetFillStyle(1001);
  gFill_Ar->SetLineColor(kAzure-9);

  gFill_Ar->Draw("F SAME");
  gPost_Ar->Draw("L SAME");

  // Upper limit lines: x in units of 10^{-40} cm^2/Ar
  TLine *lUL_Bayes  = new TLine(sigma_UL_Bayes *nuc_to_Ar*unit_scale, 0, sigma_UL_Bayes *nuc_to_Ar*unit_scale, y_max*0.85);
  TLine *lUL_FC     = new TLine(sigma_UL_FC    *nuc_to_Ar*unit_scale, 0, sigma_UL_FC    *nuc_to_Ar*unit_scale, y_max*0.70);
  TLine *lUL_Rolke  = new TLine(sigma_UL_Rolke *nuc_to_Ar*unit_scale, 0, sigma_UL_Rolke *nuc_to_Ar*unit_scale, y_max*0.70);

  lUL_Bayes->SetLineColor(kBlue+1);  lUL_Bayes->SetLineWidth(2); lUL_Bayes->SetLineStyle(1);
  lUL_FC   ->SetLineColor(kRed+1);   lUL_FC   ->SetLineWidth(2); lUL_FC   ->SetLineStyle(2);
  lUL_Rolke->SetLineColor(kGreen+2); lUL_Rolke->SetLineWidth(2); lUL_Rolke->SetLineStyle(3);

  lUL_Bayes->Draw(); lUL_FC->Draw(); lUL_Rolke->Draw();

  // MC reference line  (x in units of 10^{-40} cm^2/Ar)
  TLine *lMC = new TLine(sigma_MC_ref*nuc_to_Ar*unit_scale, 0, sigma_MC_ref*nuc_to_Ar*unit_scale, y_max*0.55);
  lMC->SetLineColor(kOrange+1); lMC->SetLineWidth(2); lMC->SetLineStyle(5);
  lMC->Draw();

  // Legend — values quoted in units of 10^{-40} cm^2/Ar
  TLegend *leg = new TLegend(0.50, 0.55, 0.92, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.034);
  leg->AddEntry(gPost_Ar, "Bayesian posterior (flat prior)", "L");
  leg->AddEntry(gFill_Ar, "90% CL region (Bayesian)", "F");
  leg->AddEntry(lUL_Bayes,
    Form("#sigma_{UL}^{Bayes} = %.2f #times 10^{-40} cm^{2}/Ar (90%% CL)",
         sigma_UL_Bayes*nuc_to_Ar*unit_scale), "L");
  leg->AddEntry(lUL_FC,
    Form("#sigma_{UL}^{FC}    = %.2f #times 10^{-40} cm^{2}/Ar (stat only)",
         sigma_UL_FC*nuc_to_Ar*unit_scale), "L");
 // leg->AddEntry(lUL_Rolke,
 //   Form("#sigma_{UL}^{Rolke} = %.2f #times 10^{-40} cm^{2}/Ar (stat+syst)",
 //        sigma_UL_Rolke*nuc_to_Ar*unit_scale), "L");
  leg->AddEntry(lMC,
    Form("#sigma_{MC ref}     = %.2f #times 10^{-40} cm^{2}/Ar",
         sigma_MC_ref*nuc_to_Ar*unit_scale), "L");
  leg->Draw();

  // Annotation
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.036);
  latex.SetTextColor(kGray+2);
  latex.DrawLatex(0.15, 0.86, "MicroBooNE NuMI RHC,  4.9 #times 10^{20} POT");
  latex.DrawLatex(0.15, 0.81, Form("#bar{#nu}_{#mu} Ar #rightarrow #mu^{+} #Sigma^{0} X,  "
                                    "N_{obs}=%d,  N_{bkg}=%.0f,  #varepsilon=%.1f%%",
                                    N_obs, N_bkg, efficiency*100));
  latex.DrawLatex(0.15, 0.76, Form("Syst. uncertainty: %.0f%%", syst_frac*100));

  c->Print("Sigma0_CrossSection_UpperLimit.png");
  c->Print("Sigma0_CrossSection_UpperLimit.pdf");

  auto t_end = std::chrono::steady_clock::now();
  double total = std::chrono::duration<double>(t_end - t_start_cs).count();
  PrintDone(Form("Sigma0_CrossSection_UpperLimit.{png,pdf}  saved"));
  std::cout << "\n\033[1;35m  All done in " << std::fixed << std::setprecision(2)
            << total << " s\033[0m\n" << std::endl;
}
