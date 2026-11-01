// ============================================================================
//  ThreshStudy.C -- per signal event: the Sigma0 decay chain in truth and whether each
//  final-state particle was reconstructed, to set the proton, pi- and photon thresholds of
//  the signal definition from this detector's reconstruction (not from other analyses).
//
//  For numubar CC Sigma0 events of a signal production (hypana ntuples):
//    truth: true vertex in the FV, muon momentum, Lambda decay mode (Lambda -> p pi- iff
//           mc_decay_pdg holds 2212 and -211), proton and pi- momenta (mc_decay_mom), Sigma0
//           photon momentum (mc_sigmazero_photon_mom);
//    reco:  a PFP truth-matched to that particle (pfp_true_pdg and pfp_true_origin: 2 Lambda
//           daughter for p / pi-, 5 Sigma0 daughter for the photon), its track score, hit count
//           and completeness. "found" = a matched PFP; "as track" / "as shower" adds the
//           trk_score split (0.5) the selections use.
//  Vetoes (signal definition, 2026-11-01): any meson at the primary vertex (|pdg| 100-999: pi+-, pi0,
//  eta, eta', rho, omega, K, ...; mesonVeto) and any strange particle other than the one Sigma0
//  (kaons are mesons; strange baryons |pdg| 3000-3999 other than a single 3212; strangeVeto).
//  Output tree "th"; run: root -l -b -q 'ThreshStudy.C+("lists/sig_genie_qe.list","th_genie.root")'
// ============================================================================
#include "HypInput.h"
#include "SigmaFV.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <cstdio>

void ThreshStudy(const char* list, const char* out)
{
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf = 0;
  while (std::getline(in, l)) if (l.size()) { ch.Add(l.c_str()); nf++; }
  HypInput H; if (!H.attach(&ch)) { printf("  not a hypana tree\n"); return; }
  // extras after attach() (TChain bind-order rule)
  std::vector<int>* dPdg = nullptr; std::vector<double>* dMom = nullptr; double gMom = -1;
  ch.SetBranchAddress("mc_decay_pdg", &dPdg); ch.SetBranchAddress("mc_decay_mom", &dMom);
  ch.SetBranchAddress("mc_sigmazero_photon_mom", &gMom);

  TFile fo(out, "RECREATE"); TTree T("th", "Sigma0 decay products: truth and reconstruction");
  int run, sub, evt, trueFV, lamPpi, mesonVeto, strangeVeto;
  double pmu, pP, pPi, pG;
  int fP, fPi, fG, tP, tPi, sG, hP, hPi, hG;   // found, as track/shower, hits of the best-matched PFP
  double cP, cPi, cG;                         // completeness of the matched PFP
  T.Branch("run", &run, "run/I"); T.Branch("sub", &sub, "sub/I"); T.Branch("evt", &evt, "evt/I");
  T.Branch("trueFV", &trueFV, "trueFV/I"); T.Branch("mesonVeto", &mesonVeto, "mesonVeto/I");
  T.Branch("strangeVeto", &strangeVeto, "strangeVeto/I");
  T.Branch("pmu", &pmu, "pmu/D"); T.Branch("lamPpi", &lamPpi, "lamPpi/I");
  T.Branch("pP", &pP, "pP/D"); T.Branch("pPi", &pPi, "pPi/D"); T.Branch("pG", &pG, "pG/D");
  T.Branch("fP", &fP, "fP/I"); T.Branch("fPi", &fPi, "fPi/I"); T.Branch("fG", &fG, "fG/I");
  T.Branch("tP", &tP, "tP/I"); T.Branch("tPi", &tPi, "tPi/I"); T.Branch("sG", &sG, "sG/I");
  T.Branch("hP", &hP, "hP/I"); T.Branch("hPi", &hPi, "hPi/I"); T.Branch("hG", &hG, "hG/I");
  T.Branch("cP", &cP, "cP/D"); T.Branch("cPi", &cPi, "cPi/D"); T.Branch("cG", &cG, "cG/D");

  const Long64_t N = ch.GetEntries();
  printf("  %s: %d files, %lld events\n", list, nf, N);
  for (Long64_t ie = 0; ie < N; ++ie) {
    ch.GetEntry(ie); H.convert();
    if (!(H.h_hyperon_pdg == 3212 && H.nu_pdg == -14 && H.ccnc == 0)) continue;
    run = H.run; sub = H.sub; evt = H.evt;
    trueFV = SigmaFV::isFiducial(H.true_vtx_x, H.true_vtx_y, H.true_vtx_z);
    pmu = H.h_lepton_mom;
    mesonVeto = 0; strangeVeto = 0; int nSig0 = 0;
    if (H.h_mc_primary_pdg) for (int pd : *H.h_mc_primary_pdg) {
      const int ap = std::abs(pd);
      if (ap >= 100 && ap <= 999) mesonVeto = 1;                       // every meson, kaons included
      if (pd == 3212) nSig0++;
      else if (ap >= 3000 && ap <= 3999) strangeVeto = 1;             // another strange baryon
    }
    if (nSig0 > 1) strangeVeto = 1;
    pP = pPi = -1; pG = gMom;
    bool hasP = false, hasPi = false;
    if (dPdg) for (size_t k = 0; k < dPdg->size(); ++k) {
      if (dPdg->at(k) == 2212) { hasP = true; pP = (dMom && k < dMom->size()) ? dMom->at(k) : -1; }
      if (dPdg->at(k) == -211) { hasPi = true; pPi = (dMom && k < dMom->size()) ? dMom->at(k) : -1; }
    }
    lamPpi = hasP && hasPi;
    fP = fPi = fG = tP = tPi = sG = 0; hP = hPi = hG = 0; cP = cPi = cG = -1;
    const int n = (int)H.trk_score.size();
    for (int i = 0; i < n; ++i) {
      const int tp = (H.h_pfp_true_pdg && i < (int)H.h_pfp_true_pdg->size()) ? H.h_pfp_true_pdg->at(i) : 0;
      const int to = (H.h_pfp_true_origin && i < (int)H.h_pfp_true_origin->size()) ? H.h_pfp_true_origin->at(i) : 0;
      const double cm = (H.h_pfp_completeness && i < (int)H.h_pfp_completeness->size()) ? H.h_pfp_completeness->at(i) : -1;
      const int hits = (int)H.pfnhits[i];
      const bool trk = H.trk_score[i] > 0.5;
      if (tp == 2212 && to == 2 && hits > hP) { fP = 1; tP = trk; hP = hits; cP = cm; }
      if (tp == -211 && to == 2 && hits > hPi) { fPi = 1; tPi = trk; hPi = hits; cPi = cm; }
      if (tp == 22 && to == 5 && hits > hG) { fG = 1; sG = !trk; hG = hits; cG = cm; }
    }
    T.Fill();
  }
  fo.Write();
  printf("  -> %s\n", out);
}
