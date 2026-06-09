// ================================================================
//  CrossSection.C
//
//  Total cross-section upper limit for Sigma0 quasi-elastic
//  production from muon antineutrinos in MicroBooNE, BNB run 3b.
//
//       nu_mu_bar + Ar → mu+ + Sigma0 + X
//       Sigma0 → Lambda + gamma,   Lambda → p + pi-
//
//  Two methods are computed and compared:
//    1. Feldman-Cousins (purely statistical)
//    2. Bayesian with Gaussian systematic convolution (recommended)
//
//  Usage:
//    root -b -l -q CrossSection.C
// ================================================================

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>

#include "TMath.h"
#include "TFeldmanCousins.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TGraph.h"

// ================================================================
//  Progress printer
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

// ================================================================
//  Likelihood marginalised over Gaussian efficiency nuisance
//
//  P(N_obs | sigma) = ∫ Poisson(N_obs | mu) × Gauss(mu | mu_0, sigma_n) dmu
//  where mu_0 = sigma × exposure  and  sigma_n = mu_0 × syst_frac.
// ================================================================
double MarginalLikelihood(double sigma_xsec, double exposure, double syst_frac,
                          int n_obs = 0)
{
  const double mu_0    = sigma_xsec * exposure;
  const double sigma_n = mu_0 * syst_frac;

  if (sigma_n < 1e-12)
    return TMath::PoissonI(n_obs, mu_0);

  const int    N_steps = 5000;
  const double mu_lo   = TMath::Max(0.0,
                           TMath::Min(mu_0, (double)n_obs) - 8.0 * sigma_n);
  const double mu_hi   = TMath::Max(mu_0, (double)n_obs)   + 8.0 * sigma_n;
  const double dmu     = (mu_hi - mu_lo) / N_steps;
  double result = 0.0;
  for (int i = 0; i < N_steps; i++) {
    const double mu = mu_lo + (i + 0.5) * dmu;
    result += TMath::PoissonI(n_obs, mu) * TMath::Gaus(mu, mu_0, sigma_n, true) * dmu;
  }
  return result;
}

// ================================================================
//  Bayesian credible interval [LL, UL] at confidence level CL.
//  Flat prior on sigma >= 0.  LL = 0 when n_obs = 0.
// ================================================================
void BayesianInterval(double exposure, double syst_frac, int n_obs, double CL,
                      double &sigma_ll, double &sigma_ul)
{
  const int    N        = 100000;
  const double sigma_max = (50.0 + 20.0 * n_obs) / exposure;
  const double dsigma   = sigma_max / N;

  std::vector<double> post(N);
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    post[i]  = MarginalLikelihood((i + 0.5) * dsigma, exposure, syst_frac, n_obs);
    norm    += post[i] * dsigma;
  }

  const double alpha_lo = (n_obs > 0) ? (1.0 - CL) / 2.0 : 0.0;
  const double alpha_hi = (n_obs > 0) ? (1.0 + CL) / 2.0 : CL;

  sigma_ll = 0.0;
  sigma_ul = sigma_max;
  double cumul = 0.0;
  bool found_ll = (alpha_lo == 0.0), found_ul = false;
  for (int i = 0; i < N; i++) {
    cumul += post[i] * dsigma / norm;
    if (!found_ll && cumul >= alpha_lo) { sigma_ll = (i + 0.5) * dsigma; found_ll = true; }
    if ( found_ll && cumul >= alpha_hi) { sigma_ul = (i + 0.5) * dsigma; found_ul = true; break; }
  }
  if (!found_ul)
    std::cerr << "\033[1;33m  WARNING: posterior truncated at sigma_max="
              << sigma_max << " cm^2\033[0m" << std::endl;
}

// ================================================================
void CrossSection()
{
  t_start_cs = std::chrono::steady_clock::now();

  std::cout << "\n\033[1;35m";
  std::cout << "=============================================================\n";
  std::cout << "  Sigma0 cross-section limit  —  MicroBooNE BNB run 3b\n";
  std::cout << "=============================================================\033[0m\n\n";

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // --------------------------------------------------------------
  // 1. Physical constants
  // --------------------------------------------------------------
  PrintStep(1, "Physical constants ...");
  const double N_A       = 6.02214076e23;  // Avogadro [mol^-1]
  const double rho_LAr   = 1.3836;         // LAr density at 87 K [g/cm^3]
  const double A_Ar      = 39.948;         // Ar molar mass [g/mol]
  const double cm2_to_pb = 1.0e36;         // 1 cm^2 = 10^36 pb
  const double unit      = 1.0e40;         // display unit: 10^{-40} cm^2/Ar
  PrintDone();

  // --------------------------------------------------------------
  // 2. Experimental inputs  ← update N_obs when data are unblinded
  // --------------------------------------------------------------
  PrintStep(2, "Experimental inputs ...");
  const double PoT         = 4.90e20;   // Protons-on-Target (BNB run 3b)
  const int    N_obs       = 1;          // Observed events (placeholder — data not yet unblinded)

  // Selection results from optimised Sigma0 analysis (BDTg > 0.96)
  const double N_pred_init = 17.869;     // MC signal: AllEvents (PoT-scaled)
  const double N_pred_sel  = 1.039;      // MC signal: after BDTg cut  (eff=5.8%, purity=95.1%)
  const double N_bkg       = 0.053;      // Total background after BDTg cut (Hyperons; Bkg=Dirt=0)

  // Systematic uncertainties
  const double syst_frac   = 0.30;      // Total (flux⊕eff⊕xsec⊕det⊕POT) for Bayesian integral
  const double bk_flux     = 0.15;
  const double bk_eff      = 0.20;
  const double bk_xsec     = 0.10;
  const double bk_det      = 0.12;
  const double bk_pot      = 0.02;
  const double syst_total  = TMath::Sqrt(bk_flux*bk_flux + bk_eff*bk_eff +
                                         bk_xsec*bk_xsec + bk_det*bk_det +
                                         bk_pot*bk_pot);

  const double CL   = 0.90;
  const double CL68 = 0.68;

  const double efficiency = N_pred_sel / N_pred_init;
  PrintDone(Form("eff=%.2f%%  N_bkg=%.3f  syst=%.0f%%", efficiency*100, N_bkg, syst_frac*100));

  // --------------------------------------------------------------
  // 3. Fiducial volume
  // --------------------------------------------------------------
  PrintStep(3, "Fiducial volume ...");
  const double dx   = 253.35 -   3.00;          // 250.35 cm
  const double dy   = 114.47 - (-112.53);        // 227.00 cm
  const double dz   = (1033.0 - 3.1) - (775.1 - 675.1);  // 929.90 cm
  const double V_FV = dx * dy * dz;
  PrintDone(Form("V_FV = %.2f cm^3", V_FV));

  // --------------------------------------------------------------
  // 4. Number of Ar nuclei in FV
  // --------------------------------------------------------------
  PrintStep(4, "Target nuclei ...");
  const double mass_FV = rho_LAr * V_FV;
  const double n_Ar    = (mass_FV / A_Ar) * N_A;
  PrintDone(Form("n_Ar = %.4e", n_Ar));

  // --------------------------------------------------------------
  // 5. Flux and exposure
  //    flux_per_POT is back-derived from N_pred_init and sigma_MC_ref;
  //    replace with the PPFX integral when the flux files are available.
  // --------------------------------------------------------------
  PrintStep(5, "Flux and exposure ...");
  const double flux_per_POT = 3.52298e-10;  // [numubar / cm^2 / POT]  ← REPLACE with PPFX value
  const double Phi          = flux_per_POT * PoT;
  const double exposure     = efficiency * Phi * n_Ar;
  const double sigma_MC_ref = N_pred_init / (Phi * n_Ar);
  PrintDone(Form("Phi=%.4e /cm^2  exposure=%.4e  sigma_MC=%.4e cm^2/Ar",
                 Phi, exposure, sigma_MC_ref));

  // --------------------------------------------------------------
  // 6. Feldman-Cousins (stat only), 90% + 68%
  // --------------------------------------------------------------
  PrintStep(6, "Feldman-Cousins ...");
  TFeldmanCousins fc(CL),   fc68(CL68);
  fc.SetMuStep(0.001);  fc68.SetMuStep(0.001);

  const double N_UL_FC      = fc.CalculateUpperLimit(N_obs, N_bkg);
  const double N_LL_FC      = fc.CalculateLowerLimit(N_obs, N_bkg);
  const double N_UL_FC68    = fc68.CalculateUpperLimit(N_obs, N_bkg);
  const double N_LL_FC68    = fc68.CalculateLowerLimit(N_obs, N_bkg);
  const double sUL_FC       = N_UL_FC   / exposure;
  const double sLL_FC       = N_LL_FC   / exposure;
  const double sUL_FC68     = N_UL_FC68 / exposure;
  const double sLL_FC68     = N_LL_FC68 / exposure;
  PrintDone(Form("90%%: [%.3f, %.3f]  68%%: [%.3f, %.3f]",
                 N_LL_FC, N_UL_FC, N_LL_FC68, N_UL_FC68));

  // --------------------------------------------------------------
  // 7. Bayesian (stat+syst), 90% + 68%
  // --------------------------------------------------------------
  PrintStep(7, "Bayesian integration ...");
  double sLL_B, sUL_B, sLL_B68, sUL_B68;
  BayesianInterval(exposure, syst_frac, N_obs, CL,   sLL_B,   sUL_B);
  BayesianInterval(exposure, syst_frac, N_obs, CL68, sLL_B68, sUL_B68);
  const double N_UL_B   = sUL_B   * exposure,  N_LL_B   = sLL_B   * exposure;
  const double N_UL_B68 = sUL_B68 * exposure,  N_LL_B68 = sLL_B68 * exposure;
  PrintDone(Form("90%%: [%.3f, %.3f]  68%%: [%.3f, %.3f]",
                 N_LL_B, N_UL_B, N_LL_B68, N_UL_B68));

  // --------------------------------------------------------------
  // 8. Asimov expected limit (N_obs = 0)
  // --------------------------------------------------------------
  PrintStep(8, "Asimov expected limit ...");
  double sLL_exp, sUL_exp;
  BayesianInterval(exposure, syst_frac, 0, CL, sLL_exp, sUL_exp);
  PrintDone(Form("sigma_exp_UL(Bayes,90%%) = %.4e cm^2/Ar  (%.2f × sigma_MC)",
                 sUL_exp, sUL_exp / sigma_MC_ref));

  // --------------------------------------------------------------
  // 9. Results summary
  // --------------------------------------------------------------
  PrintStep(9, "Results summary ...");
  const bool two_sided = (N_obs > 0);

  auto pr = [&](const char* tag, double ll, double ul, const char* unit_str) {
    if (two_sided)
      std::cout << "  " << std::left << std::setw(36) << tag
                << "[" << ll << ", " << ul << "]  " << unit_str << std::endl;
    else
      std::cout << "  " << std::left << std::setw(36) << tag
                << ul << "  " << unit_str << " (UL)" << std::endl;
  };

  std::cout << "\n--- Inputs ---\n" << std::scientific;
  std::cout << "  PoT                : " << PoT << "\n";
  std::cout << "  Phi                : " << Phi << " /cm^2\n";
  std::cout << "  n_Ar (FV)          : " << n_Ar << "\n";
  std::cout << "  mass_FV            : " << std::fixed << mass_FV/1e3 << " kg\n";
  std::cout << "  N_pred_init        : " << N_pred_init << "\n";
  std::cout << "  N_pred_sel         : " << N_pred_sel  << "  (eff=" << efficiency*100 << "%)\n";
  std::cout << "  N_bkg              : " << N_bkg << "\n";
  std::cout << "  N_obs              : " << N_obs << "\n";
  std::cout << "  sigma_MC_ref       : " << std::scientific << sigma_MC_ref << " cm^2/Ar"
            << "  (" << sigma_MC_ref*unit << " × 10^{-40})\n";

  std::cout << "\n--- Systematics ---\n" << std::fixed;
  std::cout << "  Flux               : " << bk_flux*100 << " %\n";
  std::cout << "  Efficiency         : " << bk_eff*100  << " %\n";
  std::cout << "  Cross-section      : " << bk_xsec*100 << " %\n";
  std::cout << "  Detector           : " << bk_det*100  << " %\n";
  std::cout << "  POT                : " << bk_pot*100  << " %\n";
  std::cout << "  Total (quadrature) : " << syst_total*100 << " %\n";

  std::cout << "\n--- Signal events ---\n" << std::fixed;
  pr("FC      90% CL (stat)",       N_LL_FC,  N_UL_FC,  "events");
  pr("FC      68% CL (stat)",       N_LL_FC68,N_UL_FC68,"events");
  pr("Bayesian 90% CL (stat+syst)", N_LL_B,   N_UL_B,   "events");
  pr("Bayesian 68% CL (stat+syst)", N_LL_B68, N_UL_B68, "events");

  std::cout << "\n--- sigma [10^{-40} cm^2/Ar] ---\n" << std::fixed;
  pr("FC      90% CL (stat)",       sLL_FC *unit, sUL_FC *unit, "");
  pr("FC      68% CL (stat)",       sLL_FC68*unit,sUL_FC68*unit,"");
  pr("Bayesian 90% CL (stat+syst)", sLL_B  *unit, sUL_B  *unit, "");
  pr("Bayesian 68% CL (stat+syst)", sLL_B68*unit, sUL_B68*unit, "");

  std::cout << "\n--- sigma / sigma_MC_ref ---\n" << std::fixed;
  pr("FC      90% CL",       sLL_FC/sigma_MC_ref, sUL_FC/sigma_MC_ref, "");
  pr("Bayesian 90% CL",      sLL_B/sigma_MC_ref,  sUL_B/sigma_MC_ref,  "");
  std::cout << "  Asimov UL (Bayes, 90% CL)   : " << sUL_exp*unit << " × 10^{-40} cm^2/Ar"
            << "  (" << sUL_exp/sigma_MC_ref << " × sigma_MC)\n";

  std::cout << "\n=============================================================\n\n";

  // --------------------------------------------------------------
  // 10. Posterior PDF plot
  // --------------------------------------------------------------
  PrintStep(10, "Posterior PDF plot ...");

  const int    Npts          = 2000;
  const double sigma_max_plt = 4.0 * sUL_B;
  const double dsigma        = sigma_max_plt / Npts;

  std::vector<double> vS(Npts), vP(Npts);
  double norm_post = 0.0;
  for (int i = 0; i < Npts; i++) {
    vS[i]     = (i + 0.5) * dsigma;
    vP[i]     = MarginalLikelihood(vS[i], exposure, syst_frac, N_obs);
    norm_post += vP[i] * dsigma;
  }
  for (auto &p : vP) p /= norm_post;

  // helper: build filled polygon between [lo, hi] for a given vector
  auto MakeFill = [&](double lo, double hi) -> TGraph* {
    std::vector<double> xs, ys;
    xs.push_back(lo); ys.push_back(0.0);
    for (int i = 0; i < Npts; i++) {
      if (vS[i] < lo) continue;
      if (vS[i] > hi) break;
      xs.push_back(vS[i]); ys.push_back(vP[i]);
    }
    xs.push_back(hi); ys.push_back(0.0);
    TGraph *g = new TGraph(xs.size(), xs.data(), ys.data());
    return g;
  };

  const double lo90 = two_sided ? sLL_B   : 0.0;
  const double lo68 = two_sided ? sLL_B68 : 0.0;

  TGraph *gFill90 = MakeFill(lo90, sUL_B);
  TGraph *gFill68 = MakeFill(lo68, sUL_B68);

  // Scale to display units (10^{-40} cm^2/Ar on x; PDF per same unit on y)
  auto ScaleGraph = [&](TGraph *g) {
    for (int i = 0; i < g->GetN(); i++) {
      double x, y; g->GetPoint(i, x, y);
      g->SetPoint(i, x * unit, y / unit);
    }
  };
  TGraph *gCurve = new TGraph(Npts);
  for (int i = 0; i < Npts; i++) gCurve->SetPoint(i, vS[i], vP[i]);
  ScaleGraph(gCurve);  ScaleGraph(gFill90);  ScaleGraph(gFill68);

  gFill90->SetFillColor(kAzure-9); gFill90->SetLineColor(kAzure-9); gFill90->SetFillStyle(1001);
  gFill68->SetFillColor(kAzure-4); gFill68->SetLineColor(kAzure-4); gFill68->SetFillStyle(1001);
  gCurve->SetLineColor(kBlue+1);   gCurve->SetLineWidth(3);

  TCanvas *c = new TCanvas("c_posterior", "Bayesian Posterior", 900, 650);
  c->SetLeftMargin(0.14); c->SetBottomMargin(0.13);

  const double x_max = sigma_max_plt * unit;
  const double y_max = *std::max_element(vP.begin(), vP.end()) / unit * 1.30;

  TH1D *hFrame = new TH1D("hFrame","", 100, 0, x_max);
  hFrame->SetMaximum(y_max); hFrame->SetMinimum(0.0);
  hFrame->GetXaxis()->SetTitle("#sigma(#bar{#nu}_{#mu} Ar #rightarrow #mu^{+}#Sigma^{0})"
                               " Br(#Lambda#rightarrow p#pi^{#minus})  [10^{#minus40} cm^{2}/Ar]");
  hFrame->GetYaxis()->SetTitle("Posterior probability density  [10^{40} cm^{#minus2} Ar^{#minus1}]");
  hFrame->GetXaxis()->SetTitleSize(0.040); hFrame->GetXaxis()->SetLabelSize(0.038);
  hFrame->GetYaxis()->SetTitleSize(0.040); hFrame->GetYaxis()->SetLabelSize(0.038);
  hFrame->GetXaxis()->SetTitleOffset(1.25);
  hFrame->Draw("AXIS");

  gFill90->Draw("F SAME");
  gFill68->Draw("F SAME");
  gCurve->Draw("L SAME");

  // Limit lines
  auto MakeLine = [&](double sigma_val, int col, int wid, int sty) -> TLine* {
    TLine *l = new TLine(sigma_val*unit, 0, sigma_val*unit, y_max*0.88);
    l->SetLineColor(col); l->SetLineWidth(wid); l->SetLineStyle(sty); return l;
  };
  TLine *lUL_B   = MakeLine(sUL_B,        kBlue+1,   3, 1);
  TLine *lUL_B68 = MakeLine(sUL_B68,      kBlue+1,   2, 2);
  TLine *lUL_FC  = MakeLine(sUL_FC,       kRed+1,    3, 1);
  TLine *lMC     = MakeLine(sigma_MC_ref, kOrange+1, 3, 5);

  lUL_B->Draw(); lUL_B68->Draw(); lUL_FC->Draw(); lMC->Draw();

  if (two_sided) {
    MakeLine(sLL_B,   kBlue+1, 3, 1)->Draw();
    MakeLine(sLL_B68, kBlue+1, 2, 2)->Draw();
    MakeLine(sLL_FC,  kRed+1,  3, 1)->Draw();
  }

  // Legend
  TLegend *leg = new TLegend(0.36, 0.18, 0.88, 0.60);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.029);
  leg->AddEntry(gCurve,  "Bayesian posterior (flat prior)",  "L");
  leg->AddEntry(gFill90, "90% CL credible region (Bayes)", "F");
  leg->AddEntry(gFill68, "68% CL credible region (Bayes)", "F");
  if (two_sided) {
    leg->AddEntry(lUL_B,
      Form("#sigma^{Bayes} = [%.2f, %.2f] (90%% CL, stat+syst)",
           sLL_B*unit, sUL_B*unit), "L");
    leg->AddEntry(lUL_B68,
      Form("#sigma^{Bayes} = [%.2f, %.2f] (68%% CL, stat+syst)",
           sLL_B68*unit, sUL_B68*unit), "L");
    leg->AddEntry(lUL_FC,
      Form("#sigma^{FC}    = [%.2f, %.2f] (90%% CL, stat only)",
           sLL_FC*unit, sUL_FC*unit), "L");
  } else {
    leg->AddEntry(lUL_B,
      Form("#sigma_{UL}^{Bayes} = %.2f #times10^{#minus40} cm^{2}/Ar (90%% CL)",
           sUL_B*unit), "L");
    leg->AddEntry(lUL_B68,
      Form("#sigma_{UL}^{Bayes} = %.2f #times10^{#minus40} cm^{2}/Ar (68%% CL)",
           sUL_B68*unit), "L");
    leg->AddEntry(lUL_FC,
      Form("#sigma_{UL}^{FC}    = %.2f #times10^{#minus40} cm^{2}/Ar (90%% CL, stat only)",
           sUL_FC*unit), "L");
  }
  leg->AddEntry(lMC,
    Form("#sigma_{MC}        = %.2f #times10^{#minus40} cm^{2}/Ar (GENIE ref.)",
         sigma_MC_ref*unit), "L");
  leg->Draw();

  TLatex latex;
  latex.SetNDC(); latex.SetTextSize(0.035); latex.SetTextColor(kGray+2);
  latex.DrawLatex(0.37, 0.88, "MicroBooNE BNB run 3b,  4.9#times10^{20} POT");
  latex.DrawLatex(0.37, 0.83, "#bar{#nu}_{#mu} Ar #rightarrow #mu^{+}+#Sigma^{0},"
                               "  #Sigma^{0}#rightarrow#gamma#Lambda,  #Lambda#rightarrow p#pi^{#minus}");
  latex.DrawLatex(0.37, 0.78, Form("N_{obs}=%d,  N_{bkg}=%.3f,  #varepsilon=%.1f%%",
                                    N_obs, N_bkg, efficiency*100));
  latex.DrawLatex(0.37, 0.73, Form("Syst. uncertainty (Bayesian): %.0f%%", syst_frac*100));

  c->Print("Sigma0_CrossSection_UpperLimit.png");
  c->Print("Sigma0_CrossSection_UpperLimit.pdf");

  const double total = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - t_start_cs).count();
  PrintDone("Sigma0_CrossSection_UpperLimit.{png,pdf}  saved");
  std::cout << "\n\033[1;35m  All done in " << std::fixed << std::setprecision(2)
            << total << " s\033[0m\n\n";
}
