// SigmaSelection.C
// ============================================================================
//  Event-loop selection for the MicroBooNE Sigma0 (anti-nu_mu) analysis,
//  rebuilt on standard PeLEE ntuples (tree nuselection/NeutrinoSelectionFilter).
//
//  Structure modeled on ../XSecCCPip/selection/ccpi_selection.C; physics
//  (candidate finding, cut flow) ported from the archived Sigma.C.
//
//  Reads:  NuWro hyperon-enhanced signal, RHC overlay background, dirt, and a
//          fake-data slot (NuWro hyperon sample, reco-only).
//  Writes: SigmaHists.root  — TH1D[var][sample][cut], one TDirectory per cut,
//          histogram name = VarName + SampleName + CutName  (PlotSigma.C layout).
//
//  BDT is intentionally OUT of scope here (see "BDT pass" markers); this build
//  produces the cut-based selection, cut-flow table, and control histograms.
//
//  Build & run:  root -l -b -q 'SigmaSelection.C+'
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TVector3.h"
#include "TString.h"
#include "TSystem.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

R__LOAD_LIBRARY(libTMVA)

#include "SigmaSelConfig.h"
#include "SigmaFV.h"

using namespace SigmaCfg;

// ─── Selection thresholds (ported from archived Sigma.C) ─────────────────────
static const double MUON_TRK_SCORE = 0.8;    // track-score to call a PFP the muon
static const double MUON_TRK_DIST  = 20.0;   // cm, muon start-to-vertex
static const double MUON_TRK_LEN   = 10.0;   // cm
static const double MUON_LLR       = 0.7;    // LLR PID
static const int    NSHOWERS_MAX   = 4;      // require 0 < nshowers < 4
static const int    NTRACKS_MAX    = 5;      // require ntracks < 5
static const double TRK_SHR_SPLIT  = 0.5;    // trk_score_v > 0.5 -> track, else shower

// ─── Added discriminants (optimization pass) ─────────────────────────────────
static const double TOPO_MIN       = 0.20;   // topological_score cosmic/ν veto
static const double BDT_MIN        = -0.90;  // BDTG cut (2-fold CV, full RHC bkg): ~522 sig (≥ cut-based), +20% FoM
static const double M_PROTON       = 0.938272;  // Λ/Σ⁰ mass reconstruction
static const double M_PION         = 0.139570;
// Combined proton/pion discriminant D = LLR/σ_LLR + (χ²_p−χ²_π)/σ_Δχ² (high=pion);
// scales from truth-matched RHC tracks. AUC(p vs π) 0.873, vs 0.838 for χ²_p alone.
static const double LLR_SCALE      = 0.591;
static const double DCHI_SCALE     = 124.5;

// Clamp CV event weight to a sane range, else fall back to 1 (matches ccpi).
static inline double safeWeight(double w) {
  return (std::isfinite(w) && w >= 0.0 && w <= 30.0) ? w : 1.0;
}

void SigmaSelection(const char* outname = "SigmaHists.root", bool writeBDT = true)
{
  // ── Samples ────────────────────────────────────────────────────────────────
  //  Per-file descriptor.  The signal file feeds two samples (0=Sigma signal,
  //  1=other hyperons) split by the truth flag; the other files map to one
  //  sample each.  "isData" toggles the CV-weight reset (none are real data
  //  yet: the fake-data slot is the NuWro hyperon sample, treated reco-only).
  struct FileDesc {
    std::string path;
    int   sample;        // fixed sample index, or -1 when split by truth signal
    bool  truthFV;       // also require true_nu_vtx in FV (signal/hyperon MC)
    bool  isData;        // reset ppfx/weightTune/weightSpline to 1
  };

  const std::string SIG =
    "/data/uboone/prod_hyperon_sigmazero_reco1_prodnuwro_mcc9/pelee/Hyperons_Nuwro_NoCR_pelle_1994.root";

  const std::string NUMI = "/data/uboone/new_numi_flux/";
  std::vector<FileDesc> files = {
    { SIG,                                                                  -1, true,  false }, // 0/1 signal + hyperons
    // 2 background — ALL RHC overlay runs (Run1+2+4a+4b+4c), summed POT below
    { NUMI + "Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root",               2, false, false },
    { NUMI + "Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root",               2, false, false },
    { NUMI + "Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root",              2, false, false },
    { NUMI + "Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root",              2, false, false },
    { NUMI + "Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root",              2, false, false },
    { "/data/uboone/dirt/neutrinoselection_filt_run1_dirt_overlay.root",     3, false, false }, // 3 dirt
    // ── Fake-data slot ───────────────────────────────────────────────────────
    //  NuWro hyperon sample reused as the "data" sample (reco cuts only).  A
    //  true generator-closure test needs a DIFFERENT-generator (genie) hyperon
    //  PeLEE ntuple, which does not yet exist locally — swap this one line when
    //  it does, and flip isData=true for real beam-on data.
    { SIG,                                                                   4, false, false }, // 4 fake data
  };

  // ── POT / normalization (editable block) ─────────────────────────────────────
  //  Default exposure: FULL RHC data POT (Run1+2+4), all samples scaled to it.
  //    signal NuWro generated POT : 2.1764e23  (sum of nuselection/SubRun "pot")
  //    RHC overlay POT (Run1+2+4a+4b+4c) : 9.94482e21  (/data/uboone/pot.txt RunAll_rhc:
  //      8.99715e20+5.78652e21+7.28265e20+2.37975e21+1.50509e20)
  //    RHC Run1 dirt POT          : 4.21274e20 (/data/uboone/dirt/pot.txt)
  //    Full RHC data POT (target) : 1.1526e21  (/data/uboone/pot.txt RunAll_rhc)
  //  Background now uses ALL RHC overlay runs (~11x the Run1 POT), so its statistical
  //  error is no longer Run1-sized.  NOTE: dirt is still Run1-RHC only (scaled up);
  //  add Run3b/Run4 dirt to fix that too.
  const double dataPOT   = 1.1526e21;
  const double signalPOT = 2.1764e23;
  const double bkgPOT    = 9.94482e21;   // sum of all 5 RHC overlay runs
  const double dirtPOT   = 4.21274e20;

  double Scale[Ns];
  Scale[0] = dataPOT / signalPOT;  // Sigma signal
  Scale[1] = dataPOT / signalPOT;  // other hyperons (same file)
  Scale[2] = dataPOT / bkgPOT;     // background overlay
  Scale[3] = dataPOT / dirtPOT;    // dirt
  Scale[4] = dataPOT / signalPOT;  // fake data (NuWro)

  std::cout << "POT scale factors:";
  for (int s = 0; s < Ns; s++) std::cout << "  " << Names[s] << "=" << Scale[s];
  std::cout << "\n";

  // ── Book histograms: TH1D[var][sample][cut] ─────────────────────────────────
  TFile* fout = new TFile(outname, "RECREATE");
  TH1D* Hist[Nv][Ns][Nc];
  for (int c = 0; c < Nc; c++) {
    TDirectory* d = fout->mkdir(CutName[c]);
    d->cd();
    for (int v = 0; v < Nv; v++) {
      int nb; double lo, hi; bool logy;
      Binning(v, nb, lo, hi, logy);
      for (int s = 0; s < Ns; s++) {
        TString hname = VarName[v] + Names[s] + CutName[c];
        Hist[v][s][c] = new TH1D(hname, VarName[v], nb, lo, hi);
        Hist[v][s][c]->Sumw2();
        Hist[v][s][c]->SetDirectory(d);
      }
    }
  }
  fout->cd();

  // ── BDT training-tree dump (flat per-event tree at NTracksCut preselection) ──
  //  One row per event reaching NTracksCut (candidates found, before the
  //  discriminating cuts), so the BDT can learn the final separation from the
  //  full correlated feature set.  TreeS = signal (sample 0), TreeB = all MC
  //  background (samples 1-3); the data slot (4) is excluded.
  static const char* BDTV[] = {
    "MuonPID","MuonTrkLength","MuonTrkDistVtx","MuonTrkScore","MuonNHitsY",
    "PionPID","PionTrkLength","PionTrkDistVtx","PionTrkScore","PionNHitsY","PionCaloEnergyY",
    "ProtonPID","ProtonTrkLength","ProtonTrkDistVtx","ProtonTrkScore","ProtonNHitsY","ProtonCaloEnergyY","ProtonChi2",
    "Shower1NHitsY","Shower1VtxDist","Shower1DeDxY","Shower1Energy",
    "DistMuonPion","DistProtonPion","DistProtonMuon","DistShower1Muon","DistShower1Pion","DistShower1Proton",
    "AngleMuonPion","AngleProtonPion","AngleProtonMuon","AngleShower1Muon","AngleShower1Pion","AngleShower1Proton",
    "LambdaMass","SigmaMass","TopologicalScore","NTracks","NShowers",
    "LambdaPointAngle","LambdaImpact","LambdaDecayDist","PhotonPointAngle","PhotonImpact",
    "LambdaGammaAngle","SigmaMuonAngle",
    "Pi0Mass","Pi0Gammadot","Shower2Energy",
    // Shower Moliere angle: lateral spread, orthogonal to the energy and hit-count
    // variables above and the classic shower-vs-track-stub discriminant. Appended
    // rather than inserted next to the other Shower1 variables so that no existing
    // bdtv[] index shifts. Kept LAST before the Weight spectator.
    "Shower1MoliereAvg","Weight" };
  const int NBV = sizeof(BDTV)/sizeof(BDTV[0]);   // 51
  float bdtv[51];
  float bdtIsOdd = 0;   // event-number parity, spectator for 2-fold CV split
  //  Truth/diagnostic branches (not BDT inputs): let the surviving background be
  //  characterised by interaction type, candidate true identity, π0 content, and
  //  cosmic origin, with the CV BDT score so any threshold can be applied offline.
  float tScore=0, tSample=0, tNuPdg=0, tCcnc=0, tInter=0, tNuE=0, tHasPi0=0,
        tMuPdg=0, tPrPdg=0, tPiPdg=0, tShrPdg=0, tCosmic=0;
  TFile* fbdt = nullptr; TTree* TreeS = nullptr; TTree* TreeB = nullptr;
  if (writeBDT) {
    fbdt = new TFile("BDTInput.root", "RECREATE");
    TreeS = new TTree("TreeS", "Sigma0 signal candidates");
    TreeB = new TTree("TreeB", "background candidates");
    for (auto T : {TreeS, TreeB}) {
      for (int k = 0; k < NBV; k++) T->Branch(BDTV[k], &bdtv[k]);
      T->Branch("IsOdd", &bdtIsOdd);
      T->Branch("BdtScore", &tScore);   T->Branch("Sample", &tSample);
      T->Branch("NuPdg", &tNuPdg);      T->Branch("Ccnc", &tCcnc);
      T->Branch("Interaction", &tInter);T->Branch("NuE", &tNuE);
      T->Branch("HasPi0", &tHasPi0);    T->Branch("CosmicMu", &tCosmic);
      T->Branch("MuBktPdg", &tMuPdg);   T->Branch("ProtonBktPdg", &tPrPdg);
      T->Branch("PionBktPdg", &tPiPdg); T->Branch("Shr1BktPdg", &tShrPdg);
    }
    fout->cd();
  }

  // ── BDT reader: 2-fold cross-application (unbiased) ──────────────────────────
  //  Fold A is trained on even event numbers, fold B on odd (TrainBDT.C).  Each
  //  event is scored by the fold that did NOT train on it: even→B, odd→A.  This
  //  removes the in-sample bias of applying a BDT to its own training MC.
  //  Pass-through if the weights are absent (bootstrap: run → TrainBDT.C → re-run).
  //  Variables registered in training order (BDTV[0..NBV-2]); spectators Weight
  //  and IsOdd match the trained files.
  TMVA::Reader* reader = nullptr;
  const char* wA = "SigmaBDT_A/weights/TMVASigma_BDTG.weights.xml";  // trained on even
  const char* wB = "SigmaBDT_B/weights/TMVASigma_BDTG.weights.xml";  // trained on odd
  const bool bdtReady = (gSystem->AccessPathName(wA) == kFALSE) &&
                        (gSystem->AccessPathName(wB) == kFALSE);
  if (bdtReady) {
    reader = new TMVA::Reader("!Color:!Silent");
    for (int k = 0; k < NBV - 1; k++) reader->AddVariable(BDTV[k], &bdtv[k]);
    reader->AddSpectator("Weight", &bdtv[NBV - 1]);
    reader->AddSpectator("IsOdd", &bdtIsOdd);
    reader->BookMVA("BDTG_A", wA);
    reader->BookMVA("BDTG_B", wB);
    std::cout << "BDTG 2-fold weights loaded -> applying BDTG > " << BDT_MIN
              << " (even→foldB, odd→foldA).\n";
  } else {
    std::cout << "No BDTG weights -> BDTG stage passes through (bootstrap run; "
                 "train with TrainBDT.C then re-run).\n";
  }

  // ── Cut-flow accumulators ───────────────────────────────────────────────────
  //  Selected : POT-weighted yields (physics).  RawSel : unweighted MC event
  //  counts, used to read off the true signal-efficiency per cut (denominator =
  //  RawSel[AllEvents][Sigma], i.e. true Sigma0 in FV) independent of weighting.
  double Selected[Nc][Ns] = {{0}};
  long   RawSel[Nc][Ns]   = {{0}};
  long   nSkippedRagged = 0;

  // ── Branch variables (declared once, reused across files) ───────────────────
  Int_t   run, sub, evt, nu_pdg, ccnc, interaction;
  Float_t nu_e;
  Float_t pi0_energy1_Y, pi0_energy2_Y, pi0_gammadot;   // framework π0 di-photon reco
  Float_t reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
  Float_t true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z;
  Float_t topological_score;
  Float_t ppfx_cv, weightTune, weightSpline;

  std::vector<int>*   mc_pdg              = nullptr;
  std::vector<int>*   backtracked_pdg     = nullptr;
  std::vector<float>* backtracked_overlay_purity = nullptr;  // cosmic tag (genie overlay)
  std::vector<float>* trk_score_v         = nullptr;
  std::vector<float>* trk_len_v           = nullptr;
  std::vector<float>* trk_llr_pid_score_v = nullptr;
  std::vector<float>* trk_sce_start_x_v   = nullptr;
  std::vector<float>* trk_sce_start_y_v   = nullptr;
  std::vector<float>* trk_sce_start_z_v   = nullptr;
  std::vector<float>* trk_dir_x_v         = nullptr;
  std::vector<float>* trk_dir_y_v         = nullptr;
  std::vector<float>* trk_dir_z_v         = nullptr;
  std::vector<float>* trk_calo_energy_y_v = nullptr;
  std::vector<float>* trk_pid_chipr_v     = nullptr;
  std::vector<float>* trk_pid_chipi_v     = nullptr;
  std::vector<float>* trk_energy_proton_v = nullptr;
  std::vector<float>* trk_range_muon_mom_v= nullptr;
  std::vector<float>* shr_start_x_v       = nullptr;
  std::vector<float>* shr_start_y_v       = nullptr;
  std::vector<float>* shr_start_z_v       = nullptr;
  std::vector<float>* shr_dist_v          = nullptr;
  std::vector<float>* shr_energy_y_v      = nullptr;
  std::vector<float>* shr_px_v            = nullptr;
  std::vector<float>* shr_py_v            = nullptr;
  std::vector<float>* shr_pz_v            = nullptr;
  std::vector<float>* shr_tkfit_dedx_y_v  = nullptr;
  std::vector<float>* shr_moliere_avg_v   = nullptr;
  std::vector<int>*   pfnplanehits_Y      = nullptr;

  // ── File loop ───────────────────────────────────────────────────────────────
  for (size_t i_f = 0; i_f < files.size(); i_f++) {
    const FileDesc& fd = files[i_f];

    TFile* f = TFile::Open(fd.path.c_str());
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fd.path << "\n"; if (f) delete f; continue; }
    TTree* t = nullptr;
    f->GetObject("nuselection/NeutrinoSelectionFilter", t);
    if (!t) { std::cerr << "No tree in " << fd.path << "\n"; delete f; continue; }
    std::cout << "Processing [" << i_f << "] " << fd.path << " (" << t->GetEntries() << " events)\n";

    t->SetMakeClass(1);
    t->SetBranchStatus("*", 0);
    const std::vector<const char*> branches = {
      "run","sub","evt","nu_pdg","ccnc","interaction","nu_e","mc_pdg",
      "backtracked_pdg","backtracked_overlay_purity",
      "pi0_energy1_Y","pi0_energy2_Y","pi0_gammadot",
      "reco_nu_vtx_sce_x","reco_nu_vtx_sce_y","reco_nu_vtx_sce_z",
      "true_nu_vtx_x","true_nu_vtx_y","true_nu_vtx_z","topological_score",
      "ppfx_cv","weightTune","weightSpline",
      "trk_score_v","trk_len_v","trk_llr_pid_score_v",
      "trk_sce_start_x_v","trk_sce_start_y_v","trk_sce_start_z_v",
      "trk_dir_x_v","trk_dir_y_v","trk_dir_z_v","trk_calo_energy_y_v",
      "trk_pid_chipr_v","trk_pid_chipi_v","trk_energy_proton_v","trk_range_muon_mom_v",
      "shr_start_x_v","shr_start_y_v","shr_start_z_v","shr_dist_v",
      "shr_energy_y_v","shr_px_v","shr_py_v","shr_pz_v",
      "shr_tkfit_dedx_y_v","shr_moliere_avg_v","pfnplanehits_Y"
    };
    for (auto& b : branches) t->SetBranchStatus(b, 1);

    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("sub", &sub);
    t->SetBranchAddress("evt", &evt);
    t->SetBranchAddress("nu_pdg", &nu_pdg);
    t->SetBranchAddress("ccnc", &ccnc);
    t->SetBranchAddress("interaction", &interaction);
    t->SetBranchAddress("nu_e", &nu_e);
    t->SetBranchAddress("pi0_energy1_Y", &pi0_energy1_Y);
    t->SetBranchAddress("pi0_energy2_Y", &pi0_energy2_Y);
    t->SetBranchAddress("pi0_gammadot", &pi0_gammadot);
    t->SetBranchAddress("mc_pdg", &mc_pdg);
    t->SetBranchAddress("backtracked_pdg", &backtracked_pdg);
    t->SetBranchAddress("backtracked_overlay_purity", &backtracked_overlay_purity);
    t->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x);
    t->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y);
    t->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z);
    t->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x);
    t->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y);
    t->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z);
    t->SetBranchAddress("topological_score", &topological_score);
    t->SetBranchAddress("ppfx_cv", &ppfx_cv);
    t->SetBranchAddress("weightTune", &weightTune);
    t->SetBranchAddress("weightSpline", &weightSpline);
    t->SetBranchAddress("trk_score_v", &trk_score_v);
    t->SetBranchAddress("trk_len_v", &trk_len_v);
    t->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v);
    t->SetBranchAddress("trk_sce_start_x_v", &trk_sce_start_x_v);
    t->SetBranchAddress("trk_sce_start_y_v", &trk_sce_start_y_v);
    t->SetBranchAddress("trk_sce_start_z_v", &trk_sce_start_z_v);
    t->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v);
    t->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v);
    t->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v);
    t->SetBranchAddress("trk_calo_energy_y_v", &trk_calo_energy_y_v);
    t->SetBranchAddress("trk_pid_chipr_v", &trk_pid_chipr_v);
    t->SetBranchAddress("trk_pid_chipi_v", &trk_pid_chipi_v);
    t->SetBranchAddress("trk_energy_proton_v", &trk_energy_proton_v);
    t->SetBranchAddress("trk_range_muon_mom_v", &trk_range_muon_mom_v);
    t->SetBranchAddress("shr_start_x_v", &shr_start_x_v);
    t->SetBranchAddress("shr_start_y_v", &shr_start_y_v);
    t->SetBranchAddress("shr_start_z_v", &shr_start_z_v);
    t->SetBranchAddress("shr_dist_v", &shr_dist_v);
    t->SetBranchAddress("shr_energy_y_v", &shr_energy_y_v);
    t->SetBranchAddress("shr_px_v", &shr_px_v);
    t->SetBranchAddress("shr_py_v", &shr_py_v);
    t->SetBranchAddress("shr_pz_v", &shr_pz_v);
    t->SetBranchAddress("shr_tkfit_dedx_y_v", &shr_tkfit_dedx_y_v);
    t->SetBranchAddress("shr_moliere_avg_v", &shr_moliere_avg_v);
    t->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y);

    const Long64_t nentries = t->GetEntries();
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
      t->GetEntry(jentry);

      if (jentry % 20000 == 0 || jentry == nentries - 1)
        std::cout << "\r  " << (100 * (jentry + 1) / nentries) << "%   " << std::flush;

      // ── Ragged-vector guard: every per-PFP branch must be parallel ───────────
      const std::vector<float>* fvecs[] = {
        trk_score_v, trk_len_v, trk_llr_pid_score_v,
        trk_sce_start_x_v, trk_sce_start_y_v, trk_sce_start_z_v,
        trk_dir_x_v, trk_dir_y_v, trk_dir_z_v, trk_calo_energy_y_v,
        trk_pid_chipr_v, trk_pid_chipi_v, trk_energy_proton_v, trk_range_muon_mom_v,
        shr_start_x_v, shr_start_y_v, shr_start_z_v, shr_dist_v,
        shr_energy_y_v, shr_px_v, shr_py_v, shr_pz_v,
        shr_tkfit_dedx_y_v, shr_moliere_avg_v };
      bool ragged = (trk_score_v == nullptr);
      size_t npfp = ragged ? 0 : trk_score_v->size();
      for (auto* v : fvecs) if (!v || v->size() != npfp) { ragged = true; break; }
      if (!ragged && (!pfnplanehits_Y || pfnplanehits_Y->size() != npfp)) ragged = true;
      if (!ragged && (!backtracked_pdg || backtracked_pdg->size() != npfp)) ragged = true;
      if (!ragged && (!backtracked_overlay_purity || backtracked_overlay_purity->size() != npfp)) ragged = true;
      if (ragged) { nSkippedRagged++; continue; }

      // ── CV weight (MC); reset for the data slot ──────────────────────────────
      double pf = ppfx_cv, wt = weightTune, ws = weightSpline;
      if (fd.isData) { pf = wt = ws = 1.0; }
      const double w_cv = safeWeight(pf * wt * ws);

      // ── Signal truth flag (Sigma0 in anti-nu_mu CC) ──────────────────────────
      //  Photons are NOT stored in mc_pdg, so the Sigma0->Lambda+gamma signature
      //  is identified by the Sigma0 itself (pdg 3212) in the final-state list.
      //  Decay-daughter momenta are unavailable in mc_pdg, so the old truth
      //  proton/pion momentum cuts cannot be ported.
      bool hasSigma0 = false;
      if (mc_pdg) for (int p : *mc_pdg) if (p == 3212) { hasSigma0 = true; break; }
      bool IsSignal = (nu_pdg == -14) && (ccnc == 0) && hasSigma0;

      // ── Sample assignment ────────────────────────────────────────────────────
      int s = (fd.sample >= 0) ? fd.sample : (IsSignal ? 0 : 1);
      const double W = Scale[s] * w_cv;

      // ── Truth fiducial cut (signal/hyperon MC only) ──────────────────────────
      if (fd.truthFV &&
          !SigmaFV::isFiducial(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z)) continue;

      // ── Reco vertex + fiducial flag ──────────────────────────────────────────
      TVector3 vtx(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z);
      bool sel_FV = SigmaFV::isFiducial(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z);

      const int n = (int)trk_score_v->size();

      // ── Per-PFP start position & distance-to-vertex (track vs shower) ────────
      auto isTrack = [&](int i){ return trk_score_v->at(i) > TRK_SHR_SPLIT; };
      std::vector<TVector3> pfpStart(n);
      std::vector<double>   pfpDist(n);
      for (int i = 0; i < n; i++) {
        if (isTrack(i)) {
          pfpStart[i].SetXYZ(trk_sce_start_x_v->at(i), trk_sce_start_y_v->at(i), trk_sce_start_z_v->at(i));
          pfpDist[i] = (pfpStart[i] - vtx).Mag();
        } else {
          pfpStart[i].SetXYZ(shr_start_x_v->at(i), shr_start_y_v->at(i), shr_start_z_v->at(i));
          pfpDist[i] = shr_dist_v->at(i);
        }
      }

      // ── Count tracks and showers ─────────────────────────────────────────────
      int ntracks = 0, nshowers = 0;
      for (int i = 0; i < n; i++) (isTrack(i) ? ntracks : nshowers)++;

      // ── Muon candidate: longest track passing score/dist/len/PID ─────────────
      int muon_index = -1; bool sel_muon = false;
      for (int i = 0; i < n; i++) {
        if (isTrack(i) && trk_score_v->at(i) > MUON_TRK_SCORE && pfpDist[i] < MUON_TRK_DIST &&
            trk_len_v->at(i) > MUON_TRK_LEN && trk_llr_pid_score_v->at(i) > MUON_LLR) {
          sel_muon = true;
          if (muon_index == -1 || trk_len_v->at(i) > trk_len_v->at(muon_index)) muon_index = i;
        }
      }

      // ── Joint proton/pion assignment via a combined PID discriminant ─────────
      //  D = LLR/σ_LLR + (χ²_p − χ²_π)/σ_Δχ²  (high → pion-like, low → proton-like).
      //  Combining the LLR PID with the proton/pion χ² difference separates true
      //  protons from true pions with AUC 0.873, vs 0.838 for χ²_p alone and 0.852
      //  for LLR alone.  LLR is always finite; the χ² term is dropped per-track
      //  when undefined (~30% of tracks).  Proton = argmin D, pion = argmax D.
      auto Dpi = [&](int i){
        double d = trk_llr_pid_score_v->at(i) / LLR_SCALE;
        double dchi = trk_pid_chipr_v->at(i) - trk_pid_chipi_v->at(i);
        if (std::isfinite(dchi)) d += dchi / DCHI_SCALE;
        return d;
      };
      int proton_index = -1; bool sel_proton = false;
      for (int i = 0; i < n; i++) {
        if (i == muon_index || !isTrack(i)) continue;
        sel_proton = true;
        if (proton_index == -1 || Dpi(i) < Dpi(proton_index)) proton_index = i;
      }
      int pion_index = -1; bool sel_pion = false;
      for (int i = 0; i < n; i++) {
        if (i == muon_index || i == proton_index || !isTrack(i)) continue;
        sel_pion = true;
        if (pion_index == -1 || Dpi(i) > Dpi(pion_index)) pion_index = i;
      }

      // ── Shower candidates: the two showers closest to the vertex ─────────────
      int shower1_index = -1, shower2_index = -1;
      double best1 = 1e9, best2 = 1e9;
      for (int i = 0; i < n; i++) {
        if (i == muon_index || i == pion_index || i == proton_index) continue;
        if (isTrack(i)) continue;
        double d = pfpDist[i];
        if (d < best1) { best2 = best1; shower2_index = shower1_index; best1 = d; shower1_index = i; }
        else if (d < best2) { best2 = d; shower2_index = i; }
      }

      // ── Candidate directions ─────────────────────────────────────────────────
      auto trkDir = [&](int i){ return TVector3(trk_dir_x_v->at(i), trk_dir_y_v->at(i), trk_dir_z_v->at(i)); };
      auto shrDir = [&](int i){ TVector3 d = pfpStart[i] - vtx; return (d.Mag() > 0) ? d.Unit() : TVector3(0,0,0); };

      TVector3 muonDir, pionDir, protonDir, shower1Dir;
      if (muon_index   > -1) muonDir    = trkDir(muon_index);
      if (pion_index   > -1) pionDir    = trkDir(pion_index);
      if (proton_index > -1) protonDir  = trkDir(proton_index);
      if (shower1_index> -1) shower1Dir = shrDir(shower1_index);

      // ── Inter-candidate distances ────────────────────────────────────────────
      auto dist = [&](int a, int b){ return (a>-1 && b>-1) ? (pfpStart[a]-pfpStart[b]).Mag() : 99999.; };
      double dMuPi  = dist(muon_index,    pion_index);
      double dPrPi  = dist(proton_index,  pion_index);
      double dPrMu  = dist(proton_index,  muon_index);
      double dS1Mu  = dist(shower1_index, muon_index);
      double dS1Pi  = dist(shower1_index, pion_index);
      double dS1Pr  = dist(shower1_index, proton_index);

      // ── Λ→pπ⁻ and Σ⁰→Λγ reconstruction + vertex-pointing variables ───────────
      //  Proton |p| from range-based KE (proton hyp); pion |p| from range momentum
      //  (muon hyp ≈ MIP).  Photon E from per-PFP shr_energy_y_v (MeV→GeV), dir
      //  from unit shr_px/py/pz_v.  Momentum 3-vectors are kept in the outer scope
      //  so the pointing/collinearity variables below can reuse them.
      double mLambda = -1.0, mSigma = -1.0;
      double lamPointAngle = -1, lamImpact = -1, lamDecayDist = -1;   // Λ traces to vtx?
      double phoPointAngle = -1, phoImpact = -1;                      // γ traces to vtx?
      double lamGamAngle = -1, sigMuAngle = -1;                       // Σ⁰ decay / topology
      TVector3 lambdaP(0,0,0), photonP(0,0,0);
      bool haveLambda = false, havePhoton = false;
      if (proton_index > -1 && pion_index > -1) {
        const double KEp = trk_energy_proton_v->at(proton_index);
        const double Ep  = KEp + M_PROTON;
        const double pp  = (Ep > M_PROTON) ? std::sqrt(Ep*Ep - M_PROTON*M_PROTON) : 0.0;
        const double ppi = trk_range_muon_mom_v->at(pion_index);
        const double Epi = std::sqrt(ppi*ppi + M_PION*M_PION);
        lambdaP = protonDir * pp + pionDir * ppi;
        const double E = Ep + Epi;
        if (E*E > lambdaP.Mag2()) mLambda = std::sqrt(E*E - lambdaP.Mag2());
        haveLambda = true;

        // Λ decay vertex = common origin of the p and π tracks; the Λ flies from
        // the primary vertex to here, so its momentum should be collinear with the
        // flight direction and back-project through the primary vertex.
        TVector3 lamVtx = 0.5 * (pfpStart[proton_index] + pfpStart[pion_index]);
        TVector3 flight = lamVtx - vtx;
        lamDecayDist = flight.Mag();
        if (flight.Mag() > 0 && lambdaP.Mag() > 0) {
          lamPointAngle = lambdaP.Angle(flight);
          lamImpact = (flight.Cross(lambdaP.Unit())).Mag();   // ⟂ dist of vtx to Λ line
        }

        // Photon: should originate at the primary vertex; its reco direction should
        // point from the vertex to the shower start.
        const double shrE = shower1_index > -1 ? shr_energy_y_v->at(shower1_index) : -1;
        if (shower1_index > -1 && shrE > 0 && shrE < 1e4) {
          TVector3 gdir(shr_px_v->at(shower1_index), shr_py_v->at(shower1_index), shr_pz_v->at(shower1_index));
          if (gdir.Mag() > 0) {
            gdir = gdir.Unit();
            photonP = gdir * (shrE / 1000.0);
            havePhoton = true;
            TVector3 gflight = pfpStart[shower1_index] - vtx;
            if (gflight.Mag() > 0) {
              phoPointAngle = gdir.Angle(gflight);
              phoImpact = (gflight.Cross(gdir)).Mag();
            }
            const double E2 = E + shrE / 1000.0;
            TVector3 sigP = lambdaP + photonP;
            if (E2*E2 > sigP.Mag2()) mSigma = std::sqrt(E2*E2 - sigP.Mag2());
            if (lambdaP.Mag() > 0) lamGamAngle = lambdaP.Angle(photonP);
            if (muon_index > -1 && sigP.Mag() > 0) sigMuAngle = sigP.Angle(muonDir);
          }
        }
      }

      // ── π0 rejection: di-photon mass from the framework π0 reco + 2nd shower ──
      //  Background is dominated by RES/DIS π0 (→γγ) events; the Σ⁰ signal has a
      //  single photon.  A reconstructed di-photon mass near mπ0 (135 MeV) tags
      //  the background.  Sentinels (-1 mass / -2 cos / -1 energy) when no second
      //  photon/shower is present.
      double pi0Mass = -1.0, pi0Cos = -2.0, shr2E = -1.0;
      if (pi0_energy1_Y > 0 && pi0_energy2_Y > 0 && pi0_gammadot >= -1 && pi0_gammadot <= 1) {
        pi0Cos = pi0_gammadot;
        const double mm = 2.0 * pi0_energy1_Y * pi0_energy2_Y * (1.0 - pi0_gammadot);
        if (mm > 0) pi0Mass = std::sqrt(mm) / 1000.0;   // MeV → GeV
      }
      if (shower2_index > -1) {
        const double e2 = shr_energy_y_v->at(shower2_index);
        if (e2 > 0 && e2 < 1e4) shr2E = e2 / 1000.0;
      }

      // ── Quality flags (pre-BDT) ──────────────────────────────────────────────
      //  The old hard discriminating cuts (proton-muon distance, photon shower
      //  size, Λ mass window) are gone — the BDTG below subsumes them.  Their
      //  variables are still computed and histogrammed as control plots / BDT
      //  inputs.
      bool sel_nshowers = (nshowers > 0 && nshowers < NSHOWERS_MAX);
      bool sel_ntracks  = (ntracks  < NTRACKS_MAX);
      bool sel_cosmic   = (topological_score > TOPO_MIN);

      // ── Cut chain through the NTracksCut preselection ────────────────────────
      bool Cuts[Nc];
      Cuts[0] = true;                              // AllEvents
      Cuts[1] = Cuts[0] && sel_FV;                 // FiducialV
      Cuts[2] = Cuts[1] && sel_cosmic;             // CosmicVeto (topological_score)
      Cuts[3] = Cuts[2] && sel_nshowers;           // NShowers
      Cuts[4] = Cuts[3] && sel_muon;               // MuonCandidate
      Cuts[5] = Cuts[4] && sel_proton;             // ProtonCandidate (≥1 non-muon trk)
      Cuts[6] = Cuts[5] && sel_pion;               // PionCandidate   (≥2 non-muon trk)
      Cuts[7] = Cuts[6] && sel_ntracks;            // NTracksCut

      // ── BDTG final discriminant (evaluated on the preselected sample) ────────
      //  Fill the feature vector once (also used for the training dump), then
      //  evaluate the BDT.  Sentinels (-1) handled as in training.
      double bdtg = -99.0;
      if (Cuts[7]) {
        const double shrE = shr_energy_y_v->at(shower1_index);
        bdtv[0]=trk_llr_pid_score_v->at(muon_index); bdtv[1]=trk_len_v->at(muon_index);
        bdtv[2]=pfpDist[muon_index]; bdtv[3]=trk_score_v->at(muon_index); bdtv[4]=pfnplanehits_Y->at(muon_index);
        bdtv[5]=trk_llr_pid_score_v->at(pion_index); bdtv[6]=trk_len_v->at(pion_index);
        bdtv[7]=pfpDist[pion_index]; bdtv[8]=trk_score_v->at(pion_index); bdtv[9]=pfnplanehits_Y->at(pion_index);
        bdtv[10]=trk_calo_energy_y_v->at(pion_index);
        bdtv[11]=trk_llr_pid_score_v->at(proton_index); bdtv[12]=trk_len_v->at(proton_index);
        bdtv[13]=pfpDist[proton_index]; bdtv[14]=trk_score_v->at(proton_index); bdtv[15]=pfnplanehits_Y->at(proton_index);
        bdtv[16]=trk_calo_energy_y_v->at(proton_index); bdtv[17]=trk_pid_chipr_v->at(proton_index);
        bdtv[18]=pfnplanehits_Y->at(shower1_index); bdtv[19]=pfpDist[shower1_index];
        bdtv[20]=shr_tkfit_dedx_y_v->at(shower1_index); bdtv[21]=(shrE>0?shrE/1000.0:-1.0);
        bdtv[22]=dMuPi; bdtv[23]=dPrPi; bdtv[24]=dPrMu; bdtv[25]=dS1Mu; bdtv[26]=dS1Pi; bdtv[27]=dS1Pr;
        bdtv[28]=muonDir.Angle(pionDir); bdtv[29]=protonDir.Angle(pionDir); bdtv[30]=protonDir.Angle(muonDir);
        bdtv[31]=shower1Dir.Angle(muonDir); bdtv[32]=shower1Dir.Angle(pionDir); bdtv[33]=shower1Dir.Angle(protonDir);
        bdtv[34]=mLambda; bdtv[35]=mSigma; bdtv[36]=topological_score; bdtv[37]=ntracks; bdtv[38]=nshowers;
        bdtv[39]=lamPointAngle; bdtv[40]=lamImpact; bdtv[41]=lamDecayDist;
        bdtv[42]=phoPointAngle; bdtv[43]=phoImpact; bdtv[44]=lamGamAngle; bdtv[45]=sigMuAngle;
        bdtv[46]=pi0Mass; bdtv[47]=pi0Cos; bdtv[48]=shr2E;
        bdtv[49]=(shower1_index > -1 && shr_moliere_avg_v)
                 ? shr_moliere_avg_v->at(shower1_index) : -1.0f;
        bdtv[50]=W;
        bdtIsOdd = (float)(((evt % 2) + 2) % 2);   // event-number parity (CV fold tag)
        for (int k = 0; k < NBV - 1; k++) if (!std::isfinite(bdtv[k])) bdtv[k] = -1.0f;
        // Score with the fold that did NOT train on this event: even→B, odd→A.
        if (bdtReady) bdtg = reader->EvaluateMVA(bdtIsOdd > 0.5 ? "BDTG_A" : "BDTG_B");
        if (writeBDT && s != 4) {
          bool hasPi0 = false; if (mc_pdg) for (int p : *mc_pdg) if (p == 111) { hasPi0 = true; break; }
          tScore  = bdtg;            tSample = s;
          tNuPdg  = nu_pdg;          tCcnc   = ccnc;       tInter = interaction;
          tNuE    = nu_e;            tHasPi0 = hasPi0 ? 1 : 0;
          tMuPdg  = backtracked_pdg->at(muon_index);
          tPrPdg  = backtracked_pdg->at(proton_index);
          tPiPdg  = backtracked_pdg->at(pion_index);
          tShrPdg = backtracked_pdg->at(shower1_index);
          tCosmic = (backtracked_overlay_purity->at(muon_index) > 0.5) ? 1 : 0;
          if (s == 0) TreeS->Fill(); else TreeB->Fill();
        }
      }
      bool sel_bdt = bdtReady ? (bdtg > BDT_MIN) : true;   // pass-through if no weights
      Cuts[8] = Cuts[7] && sel_bdt;                 // BDTG

      int cmax = -1;
      for (int c = 0; c < Nc; c++) { if (Cuts[c]) cmax = c; else break; }
      if (cmax < 0) continue;

      for (int c = 0; c <= cmax; c++) { Selected[c][s] += W; RawSel[c][s]++; }

      // ── Fill histograms (cumulative: every passing level 0..cmax) ────────────
      auto F = [&](int v, double x){ for (int c = 0; c <= cmax; c++) Hist[v][s][c]->Fill(x, W); };

      F(0, ntracks);
      F(1, nshowers);
      F(2, topological_score);

      if (muon_index > -1) {
        F(3, trk_llr_pid_score_v->at(muon_index));
        F(4, trk_len_v->at(muon_index));
        F(5, (pfpStart[muon_index] - vtx).Mag());
        F(6, trk_score_v->at(muon_index));
        F(7, pfnplanehits_Y->at(muon_index));
      }
      if (pion_index > -1) {
        F(8,  trk_llr_pid_score_v->at(pion_index));
        F(9,  trk_len_v->at(pion_index));
        F(10, (pfpStart[pion_index] - vtx).Mag());
        F(11, trk_score_v->at(pion_index));
        F(12, pfnplanehits_Y->at(pion_index));
        F(13, trk_calo_energy_y_v->at(pion_index));
      }
      if (proton_index > -1) {
        F(14, trk_llr_pid_score_v->at(proton_index));
        F(15, trk_len_v->at(proton_index));
        F(16, (pfpStart[proton_index] - vtx).Mag());
        F(17, trk_score_v->at(proton_index));
        F(18, pfnplanehits_Y->at(proton_index));
        F(19, trk_calo_energy_y_v->at(proton_index));
      }
      if (shower1_index > -1) {
        F(20, pfpDist[shower1_index]);
        F(21, pfnplanehits_Y->at(shower1_index));
        F(22, shr_tkfit_dedx_y_v->at(shower1_index));
        F(23, shr_moliere_avg_v->at(shower1_index));
      }
      if (shower2_index > -1) {
        F(24, pfpDist[shower2_index]);
        F(25, pfnplanehits_Y->at(shower2_index));
        F(26, shr_tkfit_dedx_y_v->at(shower2_index));
      }

      if (proton_index > -1) F(40, trk_pid_chipr_v->at(proton_index));
      if (mLambda > 0)       F(39, mLambda);
      if (mSigma  > 0)       F(41, mSigma);
      if (lamPointAngle >= 0) F(42, lamPointAngle);
      if (lamImpact     >= 0) F(43, lamImpact);
      if (lamDecayDist  >= 0) F(44, lamDecayDist);
      if (phoPointAngle >= 0) F(45, phoPointAngle);
      if (phoImpact     >= 0) F(46, phoImpact);
      if (lamGamAngle   >= 0) F(47, lamGamAngle);
      if (sigMuAngle    >= 0) F(48, sigMuAngle);
      if (bdtg > -90)         F(49, bdtg);   // CV BDTG score (for threshold scans)
      if (pi0Mass >= 0)       F(50, pi0Mass);
      if (pi0Cos  >= -1)      F(51, pi0Cos);
      if (shr2E   >= 0)       F(52, shr2E);

      if (muon_index>-1 && pion_index>-1)     { F(27, muonDir.Angle(pionDir));     F(33, dMuPi); }
      if (proton_index>-1 && pion_index>-1)   { F(28, protonDir.Angle(pionDir));   F(34, dPrPi); }
      if (proton_index>-1 && muon_index>-1)   { F(29, protonDir.Angle(muonDir));   F(35, dPrMu); }
      if (shower1_index>-1 && muon_index>-1)  { F(30, shower1Dir.Angle(muonDir));  F(36, dS1Mu); }
      if (shower1_index>-1 && pion_index>-1)  { F(31, shower1Dir.Angle(pionDir));  F(37, dS1Pi); }
      if (shower1_index>-1 && proton_index>-1){ F(32, shower1Dir.Angle(protonDir));F(38, dS1Pr); }

    } // event loop
    std::cout << "\r  done.        \n";
    delete f;
  } // file loop

  // ── Cut-flow table (POT-weighted) ───────────────────────────────────────────
  std::cout << "\n==================== POT-weighted cut flow ====================\n";
  std::cout << std::left << std::setw(18) << "Cut";
  for (int s = 0; s < Ns; s++) std::cout << std::right << std::setw(12) << Names[s];
  std::cout << "\n";
  for (int c = 0; c < Nc; c++) {
    std::cout << std::left << std::setw(18) << CutName[c];
    for (int s = 0; s < Ns; s++)
      std::cout << std::right << std::setw(12) << std::fixed << std::setprecision(2) << Selected[c][s];
    std::cout << "\n";
  }
  std::cout << "===============================================================\n";
  std::cout << "Skipped ragged entries: " << nSkippedRagged << "\n";

  // ── Signal-efficiency / significance diagnostics ────────────────────────────
  //  eff   = raw signal events at cut c / raw signal at AllEvents (true Sig in FV)
  //  S,B   = POT-weighted signal / (hyperons+background+dirt) yields
  //  S/√(S+B) is the figure of merit we optimize (keep signal, control bkg).
  const long   rawDenom = RawSel[0][0];
  std::cout << "\n=============== Signal efficiency & significance ===============\n";
  std::cout << std::left << std::setw(18) << "Cut"
            << std::right << std::setw(10) << "RawSig"
            << std::setw(9)  << "Eff[%]"
            << std::setw(12) << "S"
            << std::setw(12) << "B"
            << std::setw(11) << "Pur[%]"
            << std::setw(12) << "S/sqrt(S+B)" << "\n";
  for (int c = 0; c < Nc; c++) {
    const double S = Selected[c][0];
    const double B = Selected[c][1] + Selected[c][2] + Selected[c][3];
    const double eff = rawDenom > 0 ? 100.0 * RawSel[c][0] / rawDenom : 0.0;
    const double pur = (S + B) > 0 ? 100.0 * S / (S + B) : 0.0;
    const double sig = (S + B) > 0 ? S / std::sqrt(S + B) : 0.0;
    std::cout << std::left << std::setw(18) << CutName[c] << std::right
              << std::setw(10) << RawSel[c][0]
              << std::setw(9)  << std::fixed << std::setprecision(2) << eff
              << std::setw(12) << std::setprecision(3) << S
              << std::setw(12) << B
              << std::setw(11) << std::setprecision(3) << pur
              << std::setw(12) << std::setprecision(4) << sig << "\n";
  }
  std::cout << "================================================================\n";

  fout->Write();
  fout->Close();
  std::cout << "Wrote " << outname << "\n";

  if (writeBDT) {
    fbdt->cd();
    std::cout << "BDT training trees: TreeS=" << TreeS->GetEntries()
              << " signal, TreeB=" << TreeB->GetEntries() << " background\n";
    TreeS->Write(); TreeB->Write();
    fbdt->Close();
    std::cout << "Wrote BDTInput.root\n";
  }
  delete reader;
}
