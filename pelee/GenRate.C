// ============================================================================
//  GenRate.C -- generator-level Sigma0 rate from HepMC inputs (the 15-column
//  format both the NuWro and the GENIE productions feed into the overlay chain).
//
//  For one list of HepMC files and the POT that generated them: numubar CC
//  events with a primary hyperon (Sigma0 by default, hypPdg) and the true vertex in the fiducial volume,
//  scaled to the RHC data exposure, plus the E_nu, Q2, p_mu and p_Sigma shapes.
//  With fluxRW the NuWro events are put on the current NuMI flux with the
//  official new/old map in (E_nu, angle to the beam axis), as SigmaSelection.C
//  does; GENIE is generated on the current flux and takes no weight.
//
//  "Primary" means the Sigma0 listed at the vertex: for NuWro that is only the
//  generated hyperon when the cascade is off (kaskada_on=0), since the HepMC
//  file carries NuWro's post-FSI list.
//
//  Run:  root -l -b -q 'GenRate.C("nuwro.list", 5.37e22, true, "nuwro_bare", "genrate.root")'
// ============================================================================
#include "SigmaFV.h"
#include "HypXsecWeight.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

static const double DATAPOT = 1.1526e21;
static const char* FLUXMAP =
  "/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_oldg4rw_ratiohists/prod_ratios/ratios_g4_10_4_cvonly_fornuwrofakedata.root";

void GenRate(const char* list, double pot, bool fluxRW, const char* label, const char* out = "genrate.root",
             int hypPdg = 3212,   // 3122 for Lambda
             const char* xsecW = "")   // hyp_xsec_weights.root: correct GENIE's Pais cross section (HypXsecWeight.h)
{
  HypXsecWeight::Table xw; const bool useXW = xsecW && *xsecW && xw.load(xsecW);
  TH2D* map = nullptr;
  if (fluxRW) {
    TFile* fm = TFile::Open(FLUXMAP);
    if (fm) { map = (TH2D*)fm->Get("ratio_cv_rhc_enu_angle_numubar"); if (map) map->SetDirectory(nullptr); delete fm; }
    if (!map) { printf("  flux map not readable\n"); return; }
  }
  const TVector3 beam(0.38880857519374290, -0.0584279894529063024, 0.91946400794392302);

  TFile fo(out, "UPDATE"); TDirectory* d = fo.mkdir(label, label, true); d->cd();
  TH1D* hE  = new TH1D("enu",  ";E_{#nu} [GeV];fraction", 40, 0, 6);
  TH1D* hQ2 = new TH1D("q2",   ";Q^{2} [GeV^{2}];fraction", 40, 0, 2);
  TH1D* hPm = new TH1D("pmu",  ";p_{#mu} [GeV/c];fraction", 40, 0, 4);
  TH1D* hPs = new TH1D("psig", ";p_{hyperon} [GeV/c];fraction", 40, 0, 2);
  for (TH1D* h : {hE, hQ2, hPm, hPs}) h->Sumw2();

  double wFV = 0, wAll = 0, w2FV = 0, sumWflux = 0; long nEv = 0, nSig = 0, nFV = 0, nFiles = 0;
  std::ifstream lf(list); std::string fn;
  while (std::getline(lf, fn)) {
    if (fn.empty()) continue;
    std::ifstream in(fn); if (!in) continue; nFiles++;
    std::string line;
    while (std::getline(in, line)) {
      std::istringstream hs(line); long ev; int np;
      if (!(hs >> ev >> np)) continue;
      nEv++;
      int nupdg = 0; TLorentzVector pnu, pmu, psig, pnuc; bool haveNuc = false, haveMu = false, haveSig = false; double vx = 0, vy = 0, vz = 0;
      for (int k = 0; k < np; ++k) {
        if (!std::getline(in, line)) break;
        std::istringstream ps(line);
        int st, pdg, m1, m2, d1, d2; double px, py, pz, E, m, x, y, z, t;
        ps >> st >> pdg >> m1 >> m2 >> d1 >> d2 >> px >> py >> pz >> E >> m >> x >> y >> z >> t;
        if (st == 0 && std::abs(pdg) <= 16 && std::abs(pdg) >= 12) { nupdg = pdg; pnu.SetPxPyPzE(px, py, pz, E); vx = x; vy = y; vz = z; }
        if (st == 11 && (pdg == 2212 || pdg == 2112) && !haveNuc) { pnuc.SetPxPyPzE(px, py, pz, E); haveNuc = true; }
        if (st == 1 && pdg == -13 && !haveMu) { pmu.SetPxPyPzE(px, py, pz, E); haveMu = true; }
        if (st == 1 && pdg == hypPdg && !haveSig) { psig.SetPxPyPzE(px, py, pz, E); haveSig = true; }
      }
      if (nupdg != -14 || !haveMu || !haveSig) continue;   // numubar CC with a primary Sigma0
      double w = 1;
      if (map) {
        const double ang = pnu.Vect().Angle(beam) * 180.0 / M_PI;
        const double v = map->GetBinContent(map->FindBin(pnu.E(), ang));
        if (v > 0) w = v;
      }
      if (useXW) {
        if (!haveNuc) pnuc.SetPxPyPzE(0, 0, 0, HypXsecWeight::nucleonMass(hypPdg));
        w *= xw.weight(hypPdg, pnu, pnuc, pmu);
      }
      nSig++; wAll += w; sumWflux += w;
      if (!SigmaFV::isFiducial(vx, vy, vz)) continue;
      nFV++; wFV += w; w2FV += w * w;
      hE->Fill(pnu.E(), w); hQ2->Fill(-(pnu - pmu).M2(), w); hPm->Fill(pmu.P(), w); hPs->Fill(psig.P(), w);
    }
  }
  const double sc = DATAPOT / pot;
  printf("\n  ===== %s : %ld files, %ld events, POT %.4e, flux weight %s, xsec weight %s =====\n", label, nFiles, nEv, pot,
         map ? "on" : "off", useXW ? "on" : "off");
  printf("  numubar CC primary %d: %ld (mean flux weight %.4f); in FV %ld (%.1f%%)\n", hypPdg,
         nSig, nSig ? sumWflux / nSig : 0., nFV, nSig ? 100. * nFV / nSig : 0.);
  printf("  in-FV hyperon per %.4e POT: %.2f +- %.2f\n", DATAPOT, wFV * sc, std::sqrt(w2FV) * sc);
  printf("  <E_nu> %.3f  <Q2> %.3f  <p_mu> %.3f  <p_hyperon> %.3f\n", hE->GetMean(), hQ2->GetMean(), hPm->GetMean(), hPs->GetMean());
  for (TH1D* h : {hE, hQ2, hPm, hPs}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
  TH1D* hr = new TH1D("rate", "in-FV #Sigma^{0} per data POT", 1, 0, 1);
  hr->SetBinContent(1, wFV * sc); hr->SetBinError(1, std::sqrt(w2FV) * sc);
  d->Write(); fo.Close();
}
