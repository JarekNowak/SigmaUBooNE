// ============================================================================
//  PolWeights.C -- per-event Sigma0 polarisation weights for a signal sample.
//
//  GENIE (and NuWro through the same Geant4 chain) produce the Sigma0 unpolarised and
//  decay Sigma0 -> Lambda gamma and Lambda -> N pi isotropically. The production matrix
//  element polarises the Sigma0 (P_L, P_T from HypXsecWeights.py, ~ -0.4 each); the M1
//  decay passes it on as P_Lambda = -(P_Sigma . n) n, n = Lambda direction in the Sigma0
//  frame, and Lambda -> N pi distributes the nucleon as 1 + alpha P_Lambda . q, q its
//  direction in the Lambda frame. Weight:
//        w = 1 - alpha (P_Sigma . n)(n . q),   <w> = 1 over isotropic decays.
//  alpha = 0.747 (p pi-), 0.69 (n pi0).
//
//  Production kinematics (nubar, mu+, Sigma0) from the PeLEE ntuple of the sample; the
//  Lambda and nucleon directions from its hypana ntuples (start -> end points; momentum
//  magnitudes), matched on (run, subrun, event) within the SAME production -- the samples
//  share the run3b EXT stream, so keys are not unique across productions.
//  Checks printed: Lambda momentum in the Sigma0 frame (74.5 MeV/c) and nucleon momentum
//  in the Lambda frame (100.6 / 104.4 MeV/c), which test the rebuilt directions.
//
//  Output: "run sub evt w" (FluxWeights.h format), read by SigmaSelection.C (FileDesc::polW).
//  Run:  root -l -b -q 'PolWeights.C("<pelee.root>", "lists/sig_genie_qe.list", "pol_genie_qe.tsv")'
// ============================================================================
#include "FluxWeights.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

namespace {
const double M_SIG = 1.192642, M_LAM = 1.115683, M_P = 0.938272, M_N = 0.939565, M_MU = 0.105658;
struct Prod { TLorentzVector nu, mu, sig; };

double interp(TH2D* h, double E, double u) {
  const TAxis* ax = h->GetXaxis(); const TAxis* ay = h->GetYaxis();
  double x = std::log10(E);
  x = std::min(std::max(x, ax->GetBinCenter(1)), ax->GetBinCenter(ax->GetNbins()));
  u = std::min(std::max(u, ay->GetBinCenter(1)), ay->GetBinCenter(ay->GetNbins()));
  return h->Interpolate(x, u);
}
// u = (Q2 - Q2min)/(Q2max - Q2min) for nubar p -> mu+ Sigma0 on a free proton at rest
double uOf(double E, double Q2) {
  const double s = M_P * M_P + 2 * M_P * E, rs = std::sqrt(s);
  const double El = (s + M_MU * M_MU - M_SIG * M_SIG) / (2 * rs), pl2 = El * El - M_MU * M_MU;
  if (pl2 <= 0) return 0.5;
  const double pl = std::sqrt(pl2), Ein = (s - M_P * M_P) / (2 * rs);
  const double lo = -M_MU * M_MU + 2 * Ein * (El - pl), hi = -M_MU * M_MU + 2 * Ein * (El + pl);
  return std::min(std::max((Q2 - lo) / (hi - lo), 0.0), 1.0);
}
TVector3 dir(double sx, double sy, double sz, double ex, double ey, double ez) {
  TVector3 d(ex - sx, ey - sy, ez - sz); return d.Mag() > 0 ? d.Unit() : d;
}
}  // namespace

void PolWeights(const char* pelee, const char* hypList, const char* out,
                const char* tables = "hyp_xsec_weights.root")
{
  TFile* ft = TFile::Open(tables);
  if (!ft) { printf("  cannot open %s\n", tables); return; }
  TH2D* hPL = (TH2D*)ft->Get("PL_3212"); TH2D* hPT = (TH2D*)ft->Get("PT_3212");
  if (!hPL || !hPT) { printf("  no PL_3212/PT_3212 in %s (rerun HypXsecWeights.py)\n", tables); return; }

  // ---- production kinematics from the PeLEE ntuple ----
  std::map<unsigned long long, Prod> prod; long nDupP = 0;
  {
    TFile* fp = TFile::Open(pelee);
    TTree* t = fp ? (TTree*)fp->Get("nuselection/NeutrinoSelectionFilter") : nullptr;
    if (!t) { printf("  cannot read %s\n", pelee); return; }
    Int_t run, sub, evt, nu_pdg, ccnc; Float_t nu_e, npx, npy, npz;
    std::vector<int>* pdg = nullptr; std::vector<float> *px = nullptr, *py = nullptr, *pz = nullptr;
    t->SetBranchStatus("*", 0);
    for (auto b : {"run", "sub", "evt", "nu_pdg", "ccnc", "nu_e", "true_nu_px", "true_nu_py", "true_nu_pz",
                   "mc_pdg", "mc_px", "mc_py", "mc_pz"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("run", &run); t->SetBranchAddress("sub", &sub); t->SetBranchAddress("evt", &evt);
    t->SetBranchAddress("nu_pdg", &nu_pdg); t->SetBranchAddress("ccnc", &ccnc); t->SetBranchAddress("nu_e", &nu_e);
    t->SetBranchAddress("true_nu_px", &npx); t->SetBranchAddress("true_nu_py", &npy); t->SetBranchAddress("true_nu_pz", &npz);
    t->SetBranchAddress("mc_pdg", &pdg); t->SetBranchAddress("mc_px", &px); t->SetBranchAddress("mc_py", &py); t->SetBranchAddress("mc_pz", &pz);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (nu_pdg != -14 || ccnc != 0) continue;
      Prod p; bool hm = false, hs = false;
      for (size_t k = 0; k < pdg->size(); ++k) {
        TVector3 v(px->at(k), py->at(k), pz->at(k));
        if (pdg->at(k) == -13 && !hm) { p.mu.SetVectM(v, M_MU); hm = true; }
        if (pdg->at(k) == 3212 && !hs) { p.sig.SetVectM(v, M_SIG); hs = true; }
      }
      if (!hm || !hs) continue;
      p.nu.SetPxPyPzE(npx, npy, npz, nu_e);
      if (!prod.emplace(rseKey(run, sub, evt), p).second) nDupP++;
    }
    delete fp;
  }
  printf("  PeLEE: %zu numubar CC Sigma0 events (%ld duplicate keys, first kept)\n", prod.size(), nDupP);

  // ---- decay chain from hypana, weight ----
  TH1D hLs("hLs", ";p_{#Lambda} in #Sigma^{0} frame [GeV/c]", 100, 0, 0.3);
  TH1D hNs("hNs", ";p_{N} in #Lambda frame [GeV/c]", 100, 0, 0.3);
  TH1D hW("hW", ";w", 100, 0, 2);
  std::ofstream fo(out);
  fo << "# run sub evt w   -- Sigma0 polarisation weight, PolWeights.C (" << pelee << ")\n";
  long nHyp = 0, nMatch = 0, nNoN = 0, nOut = 0; double sw = 0, sw2 = 0;
  std::ifstream lf(hypList); std::string fn;
  while (std::getline(lf, fn)) {
    if (fn.empty()) continue;
    TFile* f = TFile::Open(fn.c_str()); if (!f || f->IsZombie()) { delete f; continue; }
    TTree* t = (TTree*)f->Get("ana/OutputTree"); if (!t) { delete f; continue; }
    UInt_t run, sub, evt; Int_t nupdg, hyp; std::string* ccnc = nullptr;
    double lm, lsx, lsy, lsz, lex, ley, lez;
    std::vector<int>* dp = nullptr;
    std::vector<double> *dm = nullptr, *dsx = nullptr, *dsy = nullptr, *dsz = nullptr, *dex = nullptr, *dey = nullptr, *dez = nullptr;
    t->SetBranchStatus("*", 0);
    auto B = [&](const char* n, void* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
    B("run", &run); B("subrun", &sub); B("event", &evt); B("mc_nu_pdg", &nupdg); B("mc_ccnc", &ccnc); B("mc_hyperon_pdg", &hyp);
    B("mc_sigmazero_lambda_mom", &lm);
    B("mc_sigmazero_lambda_start_x", &lsx); B("mc_sigmazero_lambda_start_y", &lsy); B("mc_sigmazero_lambda_start_z", &lsz);
    B("mc_sigmazero_lambda_end_x", &lex);   B("mc_sigmazero_lambda_end_y", &ley);   B("mc_sigmazero_lambda_end_z", &lez);
    B("mc_decay_pdg", &dp); B("mc_decay_mom", &dm);
    B("mc_decay_start_x", &dsx); B("mc_decay_start_y", &dsy); B("mc_decay_start_z", &dsz);
    B("mc_decay_end_x", &dex);   B("mc_decay_end_y", &dey);   B("mc_decay_end_z", &dez);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (nupdg != -14 || !ccnc || *ccnc != "CC" || hyp != 3212) continue;
      nHyp++;
      auto it = prod.find(rseKey(run, sub, evt)); if (it == prod.end()) continue;
      nMatch++;
      const Prod& p = it->second;
      int kN = -1;
      for (size_t k = 0; k < dp->size(); ++k) if (dp->at(k) == 2212 || dp->at(k) == 2112) { kN = (int)k; break; }
      if (kN < 0 || lm <= 0) { nNoN++; continue; }

      // Sigma0 polarisation in its rest frame (axes from the lab = nucleon rest frame)
      const double Q2 = -(p.nu - p.mu).M2(), E = p.nu.E();
      const double u = uOf(E, Q2);
      const TVector3 eL = p.sig.Vect().Unit();
      TVector3 eT = p.nu.Vect().Unit() - (p.nu.Vect().Unit().Dot(eL)) * eL;
      eT = eT.Mag() > 1e-9 ? eT.Unit() : eL.Orthogonal().Unit();
      const TVector3 P = interp(hPL, E, u) * eL + interp(hPT, E, u) * eT;

      // Lambda in the Sigma0 frame; nucleon in the Lambda frame reached from there
      TLorentzVector lam; lam.SetVectM(lm * dir(lsx, lsy, lsz, lex, ley, lez), M_LAM);
      const double mN = dp->at(kN) == 2212 ? M_P : M_N;
      TLorentzVector nuc; nuc.SetVectM(dm->at(kN) * dir(dsx->at(kN), dsy->at(kN), dsz->at(kN), dex->at(kN), dey->at(kN), dez->at(kN)), mN);
      const TVector3 bS = p.sig.BoostVector();
      lam.Boost(-bS); nuc.Boost(-bS);
      const TVector3 n = lam.Vect().Unit();
      nuc.Boost(-lam.BoostVector());
      hLs.Fill(lam.P()); hNs.Fill(nuc.P());
      const TVector3 q = nuc.Vect().Unit();
      const double alpha = dp->at(kN) == 2212 ? 0.747 : 0.69;
      const double w = 1 - alpha * P.Dot(n) * n.Dot(q);
      hW.Fill(w); sw += w; sw2 += w * w; nOut++;
      fo << run << " " << sub << " " << evt << " " << w << "\n";
    }
    delete f;
  }
  const double mw = nOut ? sw / nOut : 0;
  printf("  hypana: %ld numubar CC Sigma0, %ld matched, %ld without a decay nucleon (w = 1), %ld weighted\n",
         nHyp, nMatch, nNoN, nOut);
  printf("  <w> = %.4f, rms %.4f, range %.3f - %.3f\n", mw, std::sqrt(std::max(0., sw2 / std::max(1L, nOut) - mw * mw)),
         hW.GetBinLowEdge(hW.FindFirstBinAbove(0)), hW.GetBinLowEdge(hW.FindLastBinAbove(0) + 1));
  printf("  check: p_Lambda in Sigma0 frame median-ish mean %.4f rms %.4f (expect 0.0745)\n", hLs.GetMean(), hLs.GetRMS());
  printf("         p_N in Lambda frame mean %.4f rms %.4f (expect 0.1006 p, 0.1044 n)\n", hNs.GetMean(), hNs.GetRMS());
  printf("  wrote %s\n", out);
}
