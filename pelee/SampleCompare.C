// ============================================================================
//  SampleCompare.C -- truth-level comparison of Sigma0 signal samples (hypana).
//
//  For one file list: generated POT (ana/MetaTree), the in-FV numubar CC Sigma0 rate
//  scaled to the RHC data exposure, split by GENIE/NuWro interaction mode, the
//  primary-hyperon composition, and normalised kinematic histograms of the in-FV
//  signal (Q2, muon and Sigma0 momentum, Sigma0-decay photon momentum, Lambda flight).
//  Histograms go to <out> (one directory per label) for overlaying.
//
//  fluxWeights: the NuWro old->current flux table (FluxWeights.h); empty otherwise.
//  maxFiles <= 0 reads the whole list.
//
//  Run:  root -l -b -q 'SampleCompare.C("lists/sig_genie_qe.list","genie_qe","cmp.root")'
// ============================================================================
#include "SigmaFV.h"
#include "FluxWeights.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

static const double DATAPOT = 1.1526e21;

void SampleCompare(const char* list, const char* label, const char* out = "sample_compare.root",
                   int maxFiles = 0, const char* fluxWeights = "")
{
  const auto FW = loadFluxWeights(fluxWeights);
  std::vector<std::string> files;
  { std::ifstream in(list); std::string l;
    while (std::getline(in, l)) if (l.size() && (maxFiles <= 0 || (int)files.size() < maxFiles)) files.push_back(l); }

  TFile fo(out, "UPDATE");
  TDirectory* d = fo.mkdir(label, label, true); d->cd();
  auto H = [](const char* n, const char* t, int nb, double lo, double hi) {
    TH1D* h = new TH1D(n, t, nb, lo, hi); h->Sumw2(); return h; };
  TH1D* hQ2   = H("q2",     ";Q^{2} [GeV^{2}];fraction",           40, 0, 2);
  TH1D* hPmu  = H("pmu",    ";p_{#mu} [GeV/c];fraction",           40, 0, 4);
  TH1D* hPsig = H("psig",   ";p_{#Sigma^{0}} [GeV/c];fraction",    40, 0, 2);
  TH1D* hPgam = H("pgamma", ";p_{#gamma} [GeV/c];fraction",        40, 0, 0.4);
  TH1D* hFl   = H("lflight",";#Lambda flight [cm];fraction",       40, 0, 40);
  TH1D* hQ2qe = H("q2_qe",  "QEL only;Q^{2} [GeV^{2}];fraction",   40, 0, 2);
  TH1D* hPsqe = H("psig_qe","QEL only;p_{#Sigma^{0}} [GeV/c];fraction", 40, 0, 2);
  TH1D* hPmqe = H("pmu_qe", "QEL only;p_{#mu} [GeV/c];fraction",   40, 0, 4);

  double pot = 0; long nEv = 0, nFiles = 0;
  double wSig = 0, wSigFV = 0, wAll = 0;              // flux-weighted counts
  long   rSigFV = 0;
  std::map<std::string, double> byMode;               // in-FV signal by mode
  std::map<int, double> comp;                          // primary hyperon, numubar CC

  for (auto& fn : files) {
    TFile* f = TFile::Open(fn.c_str(), "READ");
    if (!f || f->IsZombie()) { delete f; continue; }
    TTree* t = (TTree*)f->Get("ana/OutputTree"); TTree* m = (TTree*)f->Get("ana/MetaTree");
    if (!t || !m) { f->Close(); delete f; continue; }
    double p = 0; m->SetBranchAddress("POT", &p);
    for (Long64_t i = 0; i < m->GetEntries(); ++i) { m->GetEntry(i); pot += p; }
    nFiles++;
    t->SetBranchStatus("*", 0);
    UInt_t run = 0, sub = 0, evt = 0; int nupdg = 0, hyp = 0;
    std::string *ccnc = nullptr, *mode = nullptr;
    double vx, vy, vz, q2, pmu, psig, pgam, lsx, lsy, lsz, lex, ley, lez;
    auto B = [&](const char* n, void* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
    B("run", &run); B("subrun", &sub); B("event", &evt); B("mc_nu_pdg", &nupdg);
    B("mc_ccnc", &ccnc); B("mc_mode", &mode); B("mc_hyperon_pdg", &hyp);
    B("mc_nu_pos_x", &vx); B("mc_nu_pos_y", &vy); B("mc_nu_pos_z", &vz);
    B("mc_nu_q2", &q2); B("mc_lepton_mom", &pmu); B("mc_hyperon_mom", &psig);
    B("mc_sigmazero_photon_mom", &pgam);
    B("mc_sigmazero_lambda_start_x", &lsx); B("mc_sigmazero_lambda_start_y", &lsy); B("mc_sigmazero_lambda_start_z", &lsz);
    B("mc_sigmazero_lambda_end_x", &lex);   B("mc_sigmazero_lambda_end_y", &ley);   B("mc_sigmazero_lambda_end_z", &lez);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); nEv++;
      double w = 1;
      if (!FW.empty()) { auto it = FW.find(rseKey(run, sub, evt)); if (it != FW.end()) w = it->second; }
      wAll += w;
      const bool cc = (nupdg == -14 && ccnc && *ccnc == "CC");
      if (cc) comp[hyp] += w;
      if (!(cc && hyp == 3212)) continue;
      wSig += w;
      if (!SigmaFV::isFiducial(vx, vy, vz)) continue;
      wSigFV += w; rSigFV++;
      const std::string md = mode ? *mode : "?";
      byMode[md] += w;
      hQ2->Fill(q2, w); hPmu->Fill(pmu, w); hPsig->Fill(psig, w); hPgam->Fill(pgam, w);
      const double fl = std::sqrt((lex-lsx)*(lex-lsx) + (ley-lsy)*(ley-lsy) + (lez-lsz)*(lez-lsz));
      if (lsx > -900 && lex > -900) hFl->Fill(fl, w);
      if (md == "QEL") { hQ2qe->Fill(q2, w); hPsqe->Fill(psig, w); hPmqe->Fill(pmu, w); }
    }
    f->Close(); delete f;
  }

  const double sc = pot > 0 ? DATAPOT / pot : 0;
  printf("\n  ===== %s : %ld files, %ld events, POT %.4e =====\n", label, nFiles, nEv, pot);
  printf("  numubar CC Sigma0: %.1f (flux-weighted), in FV %.1f (%ld raw, FV fraction %.1f%%)\n",
         wSig, wSigFV, rSigFV, wSig > 0 ? 100 * wSigFV / wSig : 0);
  printf("  in-FV Sigma0 per %.4e POT: %.2f +- %.2f\n", DATAPOT, wSigFV * sc,
         rSigFV > 0 ? wSigFV * sc / std::sqrt((double)rSigFV) : 0);
  printf("  by mode (in FV, per data POT):");
  for (auto& kv : byMode) printf("  %s %.2f (%.1f%%)", kv.first.c_str(), kv.second * sc, 100 * kv.second / wSigFV);
  printf("\n  numubar CC primary hyperon composition:");
  double tot = 0; for (auto& kv : comp) tot += kv.second;
  for (auto& kv : comp) if (kv.second > 0.005 * tot) printf("  %d %.1f%%", kv.first, 100 * kv.second / tot);
  printf("\n");

  // area-normalise the shapes; keep the rate in a 1-bin histogram
  for (TH1D* h : {hQ2, hPmu, hPsig, hPgam, hFl, hQ2qe, hPsqe, hPmqe}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
  TH1D* hr = H("rate", "in-FV #Sigma^{0} per data POT", 1, 0, 1);
  hr->SetBinContent(1, wSigFV * sc); hr->SetBinError(1, rSigFV > 0 ? wSigFV * sc / std::sqrt((double)rSigFV) : 0);
  d->Write(); fo.Close();
}
