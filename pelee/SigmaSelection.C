// SigmaSelection.C
// ============================================================================
//  Event-loop selection for the MicroBooNE Sigma0 (anti-nu_mu) analysis,
//  rebuilt on standard PeLEE ntuples (tree nuselection/NeutrinoSelectionFilter).
//
//  Structure modeled on ../XSecCCPip/selection/ccpi_selection.C; physics
//  (candidate finding, cut flow) ported from the archived Sigma.C.
//
//  Reads:  NuWro hyperon-enhanced signal, RHC overlay background, RHC dirt, the
//          run3b NuMI EXT (beam-off) data, and a fake-data slot (NuWro hyperon
//          sample, reco-only).
//  Writes: SigmaHists.root  — TH1D[var][sample][cut], one TDirectory per cut,
//          histogram name = VarName + SampleName + CutName  (PlotSigma.C layout).
//
//  BDT is intentionally OUT of scope here (see "BDT pass" markers); this build
//  produces the cut-based selection, cut-flow table, and control histograms.
//
//  Build & run:  root -l -b -q 'SigmaSelection.C+'
// ============================================================================

#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
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
#include <map>
#include <unordered_set>
#include <utility>

R__LOAD_LIBRARY(libTMVA)

#include "SigmaSelConfig.h"
#include "SigmaFV.h"
#include "PidAssign.h"

using namespace SigmaCfg;

// ─── Selection thresholds (ported from archived Sigma.C) ─────────────────────
static const double MUON_TRK_SCORE = 0.8;    // track-score to call a PFP the muon
static const double MUON_TRK_DIST  = 20.0;   // cm, muon start-to-vertex
static const double MUON_TRK_LEN   = 10.0;   // cm
static const double MUON_LLR       = 0.7;    // LLR PID
static const int    NSHOWERS_MAX   = 4;      // require 0 < nshowers < 4
static const int    NTRACKS_MAX    = 5;      // require ntracks < 5
static const double TRK_SHR_SPLIT  = 0.5;    // trk_score_v > 0.5 -> track, else shower
// Truth muon-momentum threshold for the signal definition. 0.15 GeV/c matches the
// NuMI CC1mu1pi analysis; it costs 4.2% of the Sigma0 signal.
static const double MUON_MOM_MIN_TRUE = 0.15;   // GeV/c

// ─── Added discriminants (optimization pass) ─────────────────────────────────
static const double TOPO_MIN       = 0.20;   // topological_score cosmic/ν veto
// BDTG operating point: the out-of-sample FoM peak, located by ScanPeleeBDT.C on
// BDTInput.root (scores are cross-validated, even->foldB / odd->foldA).
// Retrained 2026-10-18 with the run3b EXT in the background tree (25% of its weight)
// and on 95% of each fold instead of TMVA's default half (TrainBDT.C testFrac):
//   preselection      S= 6.370  B=26548.8  pur=0.024%  FoM=0.0391   (B incl. EXT)
//   peak     +0.759   S= 2.112  B=  382.5  pur=0.549%  FoM=0.1077   <- this
//   (scanned at every event score; ScanPeleeBDT.C's 0.002 grid gives +0.756, FoM 0.1075 -- flat)
// The 2026-10-15 weights (no EXT, half-fold training) peak on the same events at +0.690
// with S=2.330, B=547.6, FoM=0.0993: EXT in training +5%, the full fold +4% (means over
// three split seeds; this seed gives 0.1049 and 0.1077). 0.759 is
// also the peak with a 2% background systematic (FoM 0.1003).
// History: 2026-10-15 run3b NuWro retrain, +0.616 (S 2.850, B 749.9 before EXT), from
// +0.502 with the earlier weights; the retired Run 1-EXT sample gave S 4.200 only
// through normalisation (x0.677 NuMI flux weight, x0.864 lower Sigma0 rate).
// CAUTION: at the peak a 2% systematic on the background is +-7.6 events against
// 2.1 signal, so the systematic budget -- not S/sqrt(S+B) -- should ultimately set
// this. Rescan after any change to the training or the preselection.
static const double BDT_MIN        = 0.759;
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

// ─── Old → current NuMI flux reweighting (NuWro samples and the dirt) ────────
//  NuWro's flux files were converted from the OLD NuMI gsimple flux (g4numi v6,
//  the bugged beamline geometry); the GENIE overlay and the data use the CURRENT
//  flux (numi_dk2nu_g4_10_4). Every NuWro event therefore has to be weighted by
//  the official MicroBooNE ratio map "New Flux / Old (bugged) Flux", binned in
//  (E_nu [GeV], angle [deg] to the NuMI beam axis). Without it the NuWro yields
//  are too high: <w> = 0.677 over the whole sample.
//  The RHC dirt was generated on the OLD flux too, so it takes the same map. Its
//  ppfx_cv profile in true E_nu is the old-flux overlay's in every bin (anti-numu
//  0.886-0.910 for 0-5 GeV, against 0.904-1.012 on the current-flux overlay).
//  Beam axis = 3rd column of the GENIE microboone-numi-v2 rotation.
//  Source and cross-checks: ~/MicroBooNE/production/nuwro_prod/flux_rw/.
static const char* FLUXRW_MAP =
  "/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_oldg4rw_ratiohists/prod_ratios/"
  "ratios_g4_10_4_cvonly_fornuwrofakedata.root";

struct FluxRW {
  std::map<int, TH2D*> h;
  bool ok = false;
  long nUsed = 0, nEmpty = 0;
  void load(const char* horn = "rhc") {
    TFile* f = TFile::Open(FLUXRW_MAP, "READ");
    if (!f || f->IsZombie()) {
      std::cout << "  WARNING: flux ratio map not readable -> NuWro stays on the OLD flux\n";
      if (f) delete f;
      return;
    }
    const std::pair<int, const char*> m[4] =
      { {14,"numu"}, {-14,"numubar"}, {12,"nue"}, {-12,"nuebar"} };
    for (const auto& p : m) {
      TH2D* t = (TH2D*)f->Get(Form("ratio_cv_%s_enu_angle_%s", horn, p.second));
      if (t) { t->SetDirectory(nullptr); h[p.first] = t; }
    }
    f->Close(); delete f;
    ok = !h.empty();
    if (ok) std::cout << "  NuMI flux reweight map loaded (" << h.size()
                      << " species, horn " << horn << ")\n";
  }
  //  Empty map bins (and anything outside it) fall back to 1, as in
  //  flux_rw/hepmc_flux_weights.C; that is 0.2% of events.
  double w(int pdg, double E, double px, double py, double pz) {
    if (!ok) return 1.0;
    auto it = h.find(pdg);
    if (it == h.end()) return 1.0;
    TVector3 p(px, py, pz);
    if (p.Mag() <= 0.0) return 1.0;
    const double ang = p.Angle(TVector3(0.38880857519374290,
                                        -0.0584279894529063024,
                                        0.91946400794392302)) * 180.0 / TMath::Pi();
    const double v = it->second->GetBinContent(it->second->FindBin(E, ang));
    ++nUsed;
    if (v <= 0.0) { ++nEmpty; return 1.0; }
    return v;
  }
};

//  requireProtonPion=false (default) drops the ProtonCandidate/PionCandidate
//  REQUIREMENTS. Measured on this chain they remove more signal than background
//  (-31.6%/-18.3% and -52.9%/-43.5%) and more than halve the figure of merit
//  (0.0730 -> 0.0346), because roughly half the signal has no reconstructed pion.
//  Candidates are still formed when the tracks exist; when they do not the event
//  is kept, pair-dependent BDT inputs carry sentinels and HasPair flags it. The
//  two cut-flow rows then repeat the MuonCandidate row by construction.
//
//  detvarFile (detector-variation mode): run the chain on ONE GENIE overlay sample
//  instead of the sample list -- a detector-variation file or its CV. Truth signal
//  in the true FV becomes sample 0, everything else sample 2, both scaled by the
//  file's own summed SubRun POT. No BDT dump. The output also holds a per-event
//  tree `dv` (run, sub, evt, nuE, isSig, cmax, bdtg, w) so that DetVar.C can compare
//  each variation with the CV on the events they share: the variations are the
//  same generated events re-simulated, but some lost files.
//    root -l -b -q 'SigmaSelection.C+("SigmaHists_dv_CV.root",false,true,false,
//        "/data/uboone/temp/custom_pelee_ntuples/detvar/run3/numi_detvars_pion_ntuples_run3_rhc_CV_slim.root")'
void SigmaSelection(const char* outname = "SigmaHists.root", bool writeBDT = true,
                    bool requirePion = true, bool requireProtonPion = false,
                    const char* detvarFile = "")
{
  const bool detvar = detvarFile && detvarFile[0];
  if (detvar) writeBDT = false;
  // ── Samples ────────────────────────────────────────────────────────────────
  //  Per-file descriptor.  The signal file feeds two samples (0=Sigma signal,
  //  1=other hyperons) split by the truth flag; the other files map to one
  //  sample each.  "isData" marks real data (the EXT beam-off): CV weights reset,
  //  no backtracking required, duplicate events dropped. The fake-data slot is
  //  the NuWro hyperon sample, treated reco-only but still MC.
  struct FileDesc {
    std::string path;
    int   sample;        // fixed sample index, or -1 when split by truth signal
    bool  truthFV;       // also require true_nu_vtx in FV (signal/hyperon MC)
    bool  isData;        // reset ppfx/weightTune/weightSpline to 1
    bool  fluxRW;        // NuWro: weight onto the current NuMI flux (see FluxRW)
  };

  //  Signal: the run3b NuWro Sigma0 production (2500 overlay jobs on swizzled
  //  NuMI run3b ext_unbiased data, ~/MicroBooNE/production/nuwro_prod/). It
  //  replaced the Run 1-EXT prodnuwro sample on 2026-10-15: 3.3x the raw Sigma0
  //  statistics and cosmics from the same run period as the overlay background.
  //  The hypana form of the same files is lists/sig_nuwro.list.
  const std::string SIG =
    "/data/uboone/prod_nuwro_sigmazero_overlay_numi_rhc_run3b_extunbiased_mcc9/ana/Hyperons_Nuwro_run3b_WithCR_pelee_2500.root";

  const std::string NUMI = "/data/uboone/new_numi_flux/";
  const std::string DIRT = "/data/uboone/temp/custom_pelee_ntuples/dirt_overlay/";
  std::vector<FileDesc> files = {
    { SIG,                                                            -1, true,  false, true }, // 0/1 signal + hyperons
    // 2 background — ALL RHC overlay runs (Run1+2+3+4a+4b+4c), summed POT below
    { NUMI + "Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root",         2, false, false, false },
    { NUMI + "Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root",         2, false, false, false },
    { NUMI + "Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root",        2, false, false, false },
    { NUMI + "Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root",        2, false, false, false },
    { NUMI + "Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root",        2, false, false, false },
    //   Run 3 (5.518e21 POT, 1.44M events, runs 13697-18960) is split over five files
    { NUMI + "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root",      2, false, false, false },
    { NUMI + "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab.root",      2, false, false, false },
    { NUMI + "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac.root",      2, false, false, false },
    { NUMI + "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad.root",      2, false, false, false },
    { NUMI + "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae.root",      2, false, false, false },
    // 3 dirt — RHC, every run that has a sample (old flux, reweighted, see FluxRW)
    { DIRT + "numi_dirt_overlay_pion_ntuples_run1_rhc_slim.root",      3, false, false, true },
    { DIRT + "numi_dirt_overlay_pion_ntuples_run2b_rhc_slim.root",     3, false, false, true },
    { DIRT + "numi_dirt_overlay_pion_ntuples_run3b_rhc_slim.root",     3, false, false, true },
    { DIRT + "numi_dirt_overlay_pion_ntuples_run4a_rhc_slim.root",     3, false, false, true },
    { DIRT + "numi_dirt_overlay_pion_ntuples_run4b_rhc_slim.root",     3, false, false, true },
    // 5 EXT — NuMI beam-off, run3b (cosmic-only data: no truth, no CV weights)
    { "/data/uboone/EXT/neutrinoselection_filt_run3b_beamoff.root",   5, false, true,  false },
    // ── Fake-data slot ───────────────────────────────────────────────────────
    //  NuWro hyperon sample reused as the "data" sample (reco cuts only).  A
    //  true generator-closure test needs a DIFFERENT-generator (genie) hyperon
    //  PeLEE ntuple, which does not yet exist locally — swap this one line when
    //  it does, and flip isData=true for real beam-on data.
    { SIG,                                                             4, false, false, true }, // 4 fake data
  };

  // ── POT / normalization (editable block) ─────────────────────────────────────
  //  Default exposure: FULL RHC data POT (Run1+2+4), all samples scaled to it.
  //    signal NuWro generated POT : 1.3394e24 (sum of nuselection/SubRun
  //      "pot"; Hyperons_Nuwro_run3b_WithCR_pelee_2500.root, 249832 entries,
  //      57841 anti-numu CC Sigma0. Supersedes the Run 1-EXT sample
  //      Hyperons_Nuwro_WithCR_pelee_2982.root, 3.41647e23 / 64760 / 17721, and
  //      before that the NoCR 1994-file sample, 2.1764e23 / 41244 / 11316.)
  //      CAUTION: scaled by this POT the sample predicts 4.3x FEWER Sigma0 CC
  //      events than the plain RHC overlay contains at its own POT (50.4 against
  //      216.8 per 1.1526e21, same truth definition, same fiducial fractions),
  //      and 5.8x fewer than lists/sig_genie.list. See README "Sample bookkeeping".
  //    RHC overlay POT (Run1+2+4a+4b+4c) : 9.94482e21  (/data/uboone/pot.txt RunAll_rhc:
  //      8.99715e20+5.78652e21+7.28265e20+2.37975e21+1.50509e20)
  //      + Run 3 aa-ae 5.518472e21 (summed SubRun pot; added 2026-10-17, it was
  //      missing from this list) = 1.5463292e22
  //    RHC dirt POT (summed SubRun pot): run1 2.304517e20, run2b 1.254292e20,
  //      run3b 1.042393e21, run4a 7.779458e19, run4b 2.648509e20 = 1.7409194e21.
  //      Until 2026-10-17 the dirt was /data/uboone/dirt/neutrinoselection_filt_
  //      run1_dirt_overlay.root (4.21274e20), which is FHC dirt: its true-E_nu
  //      spectra match the FHC dirt samples (chi2/ndf ~1) and not the RHC ones
  //      (chi2/ndf 90-660), <E> numu 1.60 GeV against 2.15 for RHC.
  //    Full RHC data POT (target) : 1.1526e21  (/data/uboone/pot.txt RunAll_rhc)
  //  EXT: scaled by beam-on / beam-off hardware triggers. The beam-on triggers are the
  //  analyser's RHC table (xsec_analyzer configs/ext_trigger_accounting.txt): 23,534,937
  //  for 1.10823e21, taken at dataPOT with that triggers/POT. The run3b beam-off file
  //  covers 32,649,128.65 EXT gates (same table). The table's 1.10823e21 is below
  //  pot.txt's 1.1526e21 for the same runs; adopting it lowers every sample by 3.9%.
  //  The beam-off file repeats 1,655 run/sub/evt (1,619 with identical reco): data
  //  files are de-duplicated in the loop.
  const double dataPOT   = 1.1526e21;
  const double signalPOT = 1.3394e24;
  const double bkgPOT    = 1.5463292e22; // sum of all 6 RHC overlay runs
  const double dirtPOT   = 1.7409194e21;
  const double dataTrig  = dataPOT * 23534937.0 / 1.10823e21;   // 2.4477e7
  const double extGates  = 32649128.65;

  //  Detector-variation mode: the one file replaces the list (sample -2 = split
  //  into signal 0 / background 2), scaled by its own POT.
  double dvPOT = 0;
  if (detvar) {
    files = { { detvarFile, -2, false, false, false } };
    TFile* fp = TFile::Open(detvarFile);
    TTree* sr = nullptr;
    if (fp) fp->GetObject("nuselection/SubRun", sr);
    if (sr) { float p4 = 0; double p8 = 0;
      const bool isF = TString(((TLeaf*)sr->GetBranch("pot")->GetListOfLeaves()->At(0))->GetTypeName()) == "Float_t";
      if (isF) sr->SetBranchAddress("pot", &p4); else sr->SetBranchAddress("pot", &p8);
      for (Long64_t i = 0; i < sr->GetEntries(); i++) { sr->GetEntry(i); dvPOT += isF ? p4 : p8; } }
    delete fp;
    std::cout << "Detector-variation mode: " << detvarFile << "  POT " << dvPOT << "\n";
    if (dvPOT <= 0) { std::cout << "  no POT -> abort\n"; return; }
  }

  double Scale[Ns];
  Scale[0] = dataPOT / signalPOT;  // Sigma signal
  Scale[1] = dataPOT / signalPOT;  // other hyperons (same file)
  Scale[2] = dataPOT / bkgPOT;     // background overlay
  Scale[3] = dataPOT / dirtPOT;    // dirt
  Scale[4] = dataPOT / signalPOT;  // fake data (NuWro)
  Scale[5] = dataTrig / extGates;  // EXT beam-off
  if (detvar) Scale[0] = Scale[2] = dataPOT / dvPOT;

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

  //  Detector-variation per-event tree (see detvarFile above)
  TTree* dvT = nullptr;
  int dvRun = 0, dvSub = 0, dvEvt = 0, dvSig = 0, dvCmax = -1;
  float dvNuE = 0, dvBdtg = -99, dvW = 0;
  if (detvar) {
    dvT = new TTree("dv", "per-event selection, detector-variation mode");
    dvT->Branch("run", &dvRun);   dvT->Branch("sub", &dvSub);   dvT->Branch("evt", &dvEvt);
    dvT->Branch("nuE", &dvNuE);   dvT->Branch("isSig", &dvSig); dvT->Branch("cmax", &dvCmax);
    dvT->Branch("bdtg", &dvBdtg); dvT->Branch("w", &dvW);
  }

  // ── BDT training-tree dump (flat per-event tree at NTracksCut preselection) ──
  //  One row per event reaching NTracksCut (candidates found, before the
  //  discriminating cuts), so the BDT can learn the final separation from the
  //  full correlated feature set.  TreeS = signal (sample 0), TreeB = all
  //  background (samples 1-3 and the EXT, 5; Sample tells them apart); the
  //  fake-data slot (4) is excluded.
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
    "Shower1MoliereAvg",
    // Second-Shower Veto block (gLEE). SliceUnassocFrac is gLEE's "fraction of
    // unassociated hits in slice": hits in the neutrino slice belonging to no
    // reconstructed object, which is what a missed low-energy photon leaves behind.
    // SliceUnassocFrac and SliceUnassocHits are DELIBERATELY EXCLUDED. They separate
    // almost perfectly -- TMVA separation 0.86 and 0.83, six times any physics
    // variable -- but the separation comes from the PRODUCTION CHAIN, not physics.
    // slnunhits is not populated in the hyperon samples: 99.3% of signal events carry
    // zero unassociated slice hits (mean 7.6) against a mean of 900 in the standard
    // RHC overlay. This is NOT a cosmic effect: the hyperon sample does contain
    // cosmics and, by the independent measure of overlay-matched PFPs, carries MORE
    // of them than the background (15.3% against 9.5%). _closestNuCosmicDist is
    // identically zero in the same sample, so a whole class of slice/cosmic
    // bookkeeping variables differs by production chain and will fake discrimination.
    // Check any future addition for this before trusting it.
    "SliceClustFrac","SecondShowerNHit","SecondShowerDot","ShrSubClusters",
    // Recomputed replacements for the two excluded branches: the slice hits that
    // belong to no reconstructed PFP, which is the quantity gLEE's SSV actually
    // wants. Built from slnhits and pfnhits, both of which ARE populated in the
    // hyperon production. Between the hyperon sample and the RHC overlay these
    // overlap 93.4% and 95.2% with separations of 0.007 and 0.005, against 12.5%
    // and 0.774 for the slnunhits branch -- i.e. they carry the physics without
    // the production-chain artefact.
    "UnassocHitsCalc","UnassocFracCalc",
    // Whether a proton/pion pair was formed at all. With requireProtonPion=false
    // the pair is optional, every pair-dependent input above carries its sentinel
    // when it is missing, and this flags those events for the BDT.
    "HasPair",
    "Weight" };
  const int NBV = sizeof(BDTV)/sizeof(BDTV[0]);   // 58
  float bdtv[58];
  float bdtIsOdd = 0;   // event-number parity, spectator for 2-fold CV split
  //  Truth/diagnostic branches (not BDT inputs): let the surviving background be
  //  characterised by interaction type, candidate true identity, π0 content, and
  //  cosmic origin, with the CV BDT score so any threshold can be applied offline.
  float tScore=0, tSample=0, tNuPdg=0, tCcnc=0, tInter=0, tNuE=0, tHasPi0=0,
        tMuPdg=0, tPrPdg=0, tPiPdg=0, tShrPdg=0, tCosmic=0;
  //  ── All-PFP block ─────────────────────────────────────────────────────────
  //  Every reconstructed quantity for every particle in the event, rather than
  //  the four hand-picked candidate roles. PFPs are ordered by collection-plane
  //  hit count (descending) -- a stable, activity-based ordering that does not
  //  presume which particle is which -- and slots beyond the event's multiplicity
  //  carry the -1 sentinel. The preselection caps tracks at 4 and showers at 3,
  //  so 8 slots hold every surviving event.
  static const int NSLOT = 8, NPV = 15;
  static const char* PFPV[NPV] = {"score","len","llr","chipr","chipi","dvtx","nhity",
                                  "caloy","eproton","rmumom","dedxy","shrE","moliere",
                                  "cosmu","costh"};
  float pfpv[NSLOT][NPV];
  float nPFPtot = 0;
  TFile* fbdt = nullptr; TTree* TreeS = nullptr; TTree* TreeB = nullptr;
  if (writeBDT) {
    fbdt = new TFile(requirePion ? "BDTInput.root" : "BDTInput_nopion.root", "RECREATE");
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
      T->Branch("nPFPtot", &nPFPtot);
      for (int k = 0; k < NSLOT; k++)
        for (int j = 0; j < NPV; j++)
          T->Branch(Form("p%d_%s", k, PFPV[j]), &pfpv[k][j]);
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
  long   nSkippedRagged = 0, nDupData = 0, nOverlaySig = 0;
  //  The overlay's own (GENIE) signal events, counted apart from every sample:
  //  the GENIE prediction for the same signal definition and selection.
  double OvlSig[Nc] = {0}, OvlSigW2[Nc] = {0};
  long   OvlSigRaw[Nc] = {0};
  std::unordered_set<unsigned long long> seenData;   // run/sub/evt of data events read

  // ── Branch variables (declared once, reused across files) ───────────────────
  Int_t   run, sub, evt, nu_pdg, ccnc, interaction;
  Float_t nu_e;
  Float_t pi0_energy1_Y, pi0_energy2_Y, pi0_gammadot;   // framework π0 di-photon reco
  Float_t reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
  Float_t true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z;
  Float_t true_nu_px = 0.f, true_nu_py = 0.f, true_nu_pz = 0.f;  // flux reweighting
  Float_t topological_score;
  // Second-Shower-Veto inputs, after the gLEE single-photon analysis: evidence that
  // a second photon EXISTS even where Pandora never reconstructed a shower from it.
  // Our pi0 mass variable needs both photons reconstructed and therefore fires on
  // only a minority of events; these do not.
  Int_t   slnhits = 0, slnunhits = 0, secondshower_Y_nhit = 0;
  Float_t slclustfrac = -1.f, secondshower_Y_dot = -2.f;
  UInt_t  shrsubclusters0 = 0, shrsubclusters1 = 0, shrsubclusters2 = 0;
  Float_t ppfx_cv, weightTune, weightSpline;

  std::vector<int>*   mc_pdg              = nullptr;
  std::vector<float>* mc_px               = nullptr;   // truth momenta, for the
  std::vector<float>* mc_py               = nullptr;   // muon threshold in the
  std::vector<float>* mc_pz               = nullptr;   // signal definition
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
  std::vector<unsigned int>* pfnhits      = nullptr;  // hits per PFP, for the
                                                      // recomputed unassociated count

  //  Flux reweight map, loaded once; only files with fd.fluxRW use it.
  FluxRW frw; frw.load("rhc");

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
      "mc_px","mc_py","mc_pz",
      "backtracked_pdg","backtracked_overlay_purity",
      "pi0_energy1_Y","pi0_energy2_Y","pi0_gammadot",
      "reco_nu_vtx_sce_x","reco_nu_vtx_sce_y","reco_nu_vtx_sce_z",
      "true_nu_vtx_x","true_nu_vtx_y","true_nu_vtx_z","topological_score",
      "true_nu_px","true_nu_py","true_nu_pz",
      "ppfx_cv","weightTune","weightSpline",
      "trk_score_v","trk_len_v","trk_llr_pid_score_v",
      "trk_sce_start_x_v","trk_sce_start_y_v","trk_sce_start_z_v",
      "trk_dir_x_v","trk_dir_y_v","trk_dir_z_v","trk_calo_energy_y_v",
      "trk_pid_chipr_v","trk_pid_chipi_v","trk_energy_proton_v","trk_range_muon_mom_v",
      "shr_start_x_v","shr_start_y_v","shr_start_z_v","shr_dist_v",
      "shr_energy_y_v","shr_px_v","shr_py_v","shr_pz_v",
      "shr_tkfit_dedx_y_v","shr_moliere_avg_v","pfnplanehits_Y","pfnhits",
      "slnhits","slnunhits","slclustfrac","secondshower_Y_nhit","secondshower_Y_dot",
      "shrsubclusters0","shrsubclusters1","shrsubclusters2"
    };
    for (auto& b : branches) if (t->GetBranch(b)) t->SetBranchStatus(b, 1);

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
    if (t->GetBranch("mc_px")) t->SetBranchAddress("mc_px", &mc_px);
    if (t->GetBranch("mc_py")) t->SetBranchAddress("mc_py", &mc_py);
    if (t->GetBranch("mc_pz")) t->SetBranchAddress("mc_pz", &mc_pz);
    t->SetBranchAddress("backtracked_pdg", &backtracked_pdg);
    t->SetBranchAddress("backtracked_overlay_purity", &backtracked_overlay_purity);
    t->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x);
    t->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y);
    t->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z);
    t->SetBranchAddress("true_nu_px", &true_nu_px);
    t->SetBranchAddress("true_nu_py", &true_nu_py);
    t->SetBranchAddress("true_nu_pz", &true_nu_pz);
    t->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x);
    t->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y);
    t->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z);
    t->SetBranchAddress("topological_score", &topological_score);
    // CV weights: absent from the data files, where they are reset anyway
    ppfx_cv = weightTune = weightSpline = 1.f;
    if (t->GetBranch("ppfx_cv"))      t->SetBranchAddress("ppfx_cv", &ppfx_cv);
    if (t->GetBranch("weightTune"))   t->SetBranchAddress("weightTune", &weightTune);
    if (t->GetBranch("weightSpline")) t->SetBranchAddress("weightSpline", &weightSpline);
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
    // Guarded: absent from some ntuple generations; the sentinels above then stand.
    if (t->GetBranch("slnhits"))             t->SetBranchAddress("slnhits", &slnhits);
    if (t->GetBranch("slnunhits"))           t->SetBranchAddress("slnunhits", &slnunhits);
    if (t->GetBranch("slclustfrac"))         t->SetBranchAddress("slclustfrac", &slclustfrac);
    if (t->GetBranch("secondshower_Y_nhit")) t->SetBranchAddress("secondshower_Y_nhit", &secondshower_Y_nhit);
    if (t->GetBranch("secondshower_Y_dot"))  t->SetBranchAddress("secondshower_Y_dot", &secondshower_Y_dot);
    if (t->GetBranch("shrsubclusters0"))     t->SetBranchAddress("shrsubclusters0", &shrsubclusters0);
    if (t->GetBranch("shrsubclusters1"))     t->SetBranchAddress("shrsubclusters1", &shrsubclusters1);
    if (t->GetBranch("shrsubclusters2"))     t->SetBranchAddress("shrsubclusters2", &shrsubclusters2);
    t->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y);
    if (t->GetBranch("pfnhits")) t->SetBranchAddress("pfnhits", &pfnhits);

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
      //  Data carry no backtracking: both vectors are EMPTY there whenever the event
      //  has PFPs (all 374,958 neutrino-slice events of the run3b beam-off file), so
      //  requiring them parallel would silently drop every candidate EXT event.
      const bool needBT = !fd.isData;
      if (!ragged && needBT && (!backtracked_pdg || backtracked_pdg->size() != npfp)) ragged = true;
      if (!ragged && needBT && (!backtracked_overlay_purity || backtracked_overlay_purity->size() != npfp)) ragged = true;
      if (ragged) { nSkippedRagged++; continue; }
      if (fd.isData && !seenData.insert(((unsigned long long)run << 40) |
                                        ((unsigned long long)sub << 20) |
                                        (unsigned long long)evt).second) { nDupData++; continue; }

      // ── CV weight (MC); reset for the data slot ──────────────────────────────
      //  NuWro samples and the dirt carry an extra factor putting them on the
      //  current NuMI flux (fd.fluxRW); it is applied to the fake-data slot too,
      //  since that slot stands in for data and data is on the current flux.
      double pf = ppfx_cv, wt = weightTune, ws = weightSpline;
      if (fd.isData) { pf = wt = ws = 1.0; }
      const double w_flux = fd.fluxRW
        ? frw.w(nu_pdg, nu_e, true_nu_px, true_nu_py, true_nu_pz) : 1.0;
      const double w_cv = safeWeight(pf * wt * ws) * w_flux;

      // ── Signal truth flag (Sigma0 in anti-nu_mu CC) ──────────────────────────
      //  Photons are NOT stored in mc_pdg, so the Sigma0->Lambda+gamma signature
      //  is identified by the Sigma0 itself (pdg 3212) in the final-state list.
      //  Decay-daughter momenta are unavailable in mc_pdg, so the old truth
      //  proton/pion momentum cuts cannot be ported.
      //  SIGNAL DEFINITION
      //
      //  anti-nu_mu CC, a Sigma0 in the final state, the muon above threshold, and no
      //  meson other than the Lambda decay products.
      //
      //  What can and cannot be required here is set by the truth record, which in
      //  these PeLEE ntuples stops at the PRIMARY vertex: mc_pdg lists the Sigma0
      //  itself (1.00 per event) and the muon (1.00), but NOT the Sigma0 -> Lambda
      //  gamma -> p pi- gamma chain. Measured on the signal sample: nproton averages
      //  0.06 and is zero in 94.4% of events, npion is zero in 100%, PDG 22 never
      //  appears at all, and zero of 17721 events have a backtracked pi- PFP -- yet
      //  every one of them contains a proton, a pi- and a photon from the decay.
      //  Proton, pion and photon thresholds are therefore IMPOSSIBLE to impose from
      //  this record, and the momentum cuts of the archived OutputTree analysis cannot
      //  be ported. Doing so needs a truth record that follows the decay chain, which
      //  is what the HyperonProduction analyzer provides.
      bool hasSigma0 = false;
      double pmu_true = -1.0;
      bool hasOtherMeson = false;
      if (mc_pdg) {
        for (size_t k = 0; k < mc_pdg->size(); ++k) {
          const int pd = mc_pdg->at(k), ap = std::abs(pd);
          if (pd == 3212) hasSigma0 = true;
          if (ap == 13 && mc_px && k < mc_px->size())
            pmu_true = std::sqrt( mc_px->at(k)*mc_px->at(k)
                                + mc_py->at(k)*mc_py->at(k)
                                + mc_pz->at(k)*mc_pz->at(k) );
          //  Any meson at the primary vertex. The Lambda decay pi- is not in this
          //  list (see above), so this vetoes genuinely EXTRA mesons and does not
          //  reject the signal's own decay products. pi0, eta, and the kaons.
          if (ap == 111 || ap == 221 || ap == 331 ||
              ap == 130 || ap == 310 || ap == 311 || ap == 321) hasOtherMeson = true;
        }
      }
      //  Muon threshold. The muon IS in the record, mean momentum 1.005 GeV/c, and
      //  only 4.2% of signal muons fall below 0.15 GeV/c -- the same threshold the
      //  CC1mu1pi NuMI analysis uses, chosen there for reconstruction rather than
      //  physics reasons and adopted here for consistency.
      //  If mc_px is absent from a sample the momentum cannot be formed, and
      //  pmu_true stays at its -1 sentinel. Requiring the threshold then would silently
      //  reject EVERY event as non-signal, so fall back to not applying it and say so
      //  once. Preferring a wrong-but-loud failure to a wrong-and-silent one.
      static bool warnedNoMuP = false;
      const bool haveMuP = ( pmu_true >= 0.0 );
      if ( !haveMuP && !warnedNoMuP ) {
        std::cout << "  WARNING: mc_px/py/pz unavailable -- the truth muon-momentum "
                     "threshold in the signal definition is NOT being applied.\n";
        warnedNoMuP = true;
      }
      const bool muAboveThreshold = haveMuP ? ( pmu_true > MUON_MOM_MIN_TRUE ) : true;

      //  The meson veto currently removes NOTHING: no signal event in this NuWro
      //  hyperon sample contains a pi0, kaon or eta. It is imposed anyway so that the
      //  measured quantity is defined rather than inherited from whatever the
      //  generator happened to produce, and because a different generator -- the GENIE
      //  overlay samples in preparation -- may well produce Sigma0 with an associated
      //  pi0.
      bool IsSignal = (nu_pdg == -14) && (ccnc == 0) && hasSigma0
                      && muAboveThreshold && !hasOtherMeson;

      // ── Sample assignment ────────────────────────────────────────────────────
      int s = (fd.sample >= 0) ? fd.sample : (IsSignal ? 0 : 1);
      if (fd.sample == -2)    // detector variation: signal in the true FV, or background
        s = (IsSignal && SigmaFV::isFiducial(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z)) ? 0 : 2;
      const double W = Scale[s] * w_cv;

      // ── Truth fiducial cut (signal/hyperon MC only) ──────────────────────────
      if (fd.truthFV &&
          !SigmaFV::isFiducial(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z)) continue;

      //  The inclusive overlay contains GENIE signal events too. The signal is
      //  modelled by the dedicated sample, so they are not background: keep them
      //  out of B, the BDT trees and the histograms, as the hypana chain does
      //  (IslandProcess.C, isSig), and count them in OvlSig as the GENIE signal
      //  prediction. The BDT never trains on them, so their score is unbiased.
      //  Signal outside the true FV stays, as background.
      const bool ovlSig = (s == 2 && fd.sample == 2 && IsSignal &&
          SigmaFV::isFiducial(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z));
      if (ovlSig) nOverlaySig++;

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

      // ── Joint proton/pion assignment (see PidAssign.h) ──────────────────────
      //  The old rule ranked non-muon tracks by a single monotonic discriminant
      //  D = LLR/σ_LLR + (χ²_p − χ²_π)/σ_Δχ² and took its two extremes. Measured
      //  against HyperonProduction truth, the pion slot that produces is a true
      //  π⁻ only 33.5% of the time: photons (20.4%) and muons (15.4%) sit at the
      //  same high-D end as a real pion. Scoring each track separately for the
      //  two roles and taking the best ORDERED PAIR lifts that to 42.5% with the
      //  proton unchanged, at no cost in efficiency -- the surviving events are
      //  identical, only the indices change.
      std::vector<int> cand;
      for (int i = 0; i < n; i++) {
        if (i == muon_index || !isTrack(i)) continue;
        cand.push_back(i);
      }
      const bool sel_proton = !cand.empty();       // ProtonCandidate: ≥1 non-muon trk
      const bool sel_pion   = cand.size() >= 2;    // PionCandidate:   ≥2 non-muon trk
      int proton_index = -1, pion_index = -1;
      if (sel_pion) {
        std::vector<double> dvv(n, 0.0), nhv(n, 0.0);
        for (int i = 0; i < n; i++) {
          dvv[i] = pfpDist[i];
          nhv[i] = (pfnplanehits_Y && i < (int)pfnplanehits_Y->size())
                     ? (double)pfnplanehits_Y->at(i) : 0.0;
        }
        PidAssign::assign(cand, *trk_score_v, *trk_llr_pid_score_v,
                          *trk_pid_chipr_v, *trk_pid_chipi_v, *trk_len_v,
                          dvv, nhv, proton_index, pion_index);
      } else if (sel_proton) {
        proton_index = cand[0];                    // only one track can fill the slot
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
      Cuts[5] = Cuts[4] && (requireProtonPion ? sel_proton : true);  // ProtonCandidate
      Cuts[6] = Cuts[5] && (requireProtonPion ? sel_pion   : true);  // PionCandidate
      Cuts[7] = Cuts[6] && sel_ntracks;            // NTracksCut

      // Preselection gate for the BDT dump and evaluation. With requirePion=false
      // the third-track requirement is dropped: events with no pion candidate are
      // KEPT and their pion-derived inputs carry sentinels, so the comparison
      // isolates the effect of the CUT rather than of a different feature set.
      const bool presel = requirePion ? Cuts[7] : (Cuts[5] && sel_ntracks);

      // ── BDTG final discriminant (evaluated on the preselected sample) ────────
      //  Fill the feature vector once (also used for the training dump), then
      //  evaluate the BDT.  Sentinels (-1) handled as in training.
      double bdtg = -99.0;
      if (presel) {
        const double shrE = shr_energy_y_v->at(shower1_index);
        bdtv[0]=trk_llr_pid_score_v->at(muon_index); bdtv[1]=trk_len_v->at(muon_index);
        bdtv[2]=pfpDist[muon_index]; bdtv[3]=trk_score_v->at(muon_index); bdtv[4]=pfnplanehits_Y->at(muon_index);
        // Pion inputs go through a guarded accessor: with requirePion=false the
        // pion candidate may not exist, and the sentinel is what the BDT sees.
        auto PATf=[&](std::vector<float>* v){
          return (pion_index>-1 && v && pion_index<(int)v->size()) ? (double)v->at(pion_index) : -1.0; };
        auto PATi=[&](std::vector<int>* v){
          return (pion_index>-1 && v && pion_index<(int)v->size()) ? (double)v->at(pion_index) : -1.0; };
        bdtv[5]=PATf(trk_llr_pid_score_v); bdtv[6]=PATf(trk_len_v);
        bdtv[7]=(pion_index>-1?pfpDist[pion_index]:-1.0); bdtv[8]=PATf(trk_score_v); bdtv[9]=PATi(pfnplanehits_Y);
        bdtv[10]=PATf(trk_calo_energy_y_v);
        // Proton inputs guarded the same way: the slot can be empty once the
        // candidate requirements are optional.
        auto RATf=[&](std::vector<float>* v){
          return (proton_index>-1 && v && proton_index<(int)v->size()) ? (double)v->at(proton_index) : -1.0; };
        auto RATi=[&](std::vector<int>* v){
          return (proton_index>-1 && v && proton_index<(int)v->size()) ? (double)v->at(proton_index) : -1.0; };
        bdtv[11]=RATf(trk_llr_pid_score_v); bdtv[12]=RATf(trk_len_v);
        bdtv[13]=(proton_index>-1?pfpDist[proton_index]:-1.0); bdtv[14]=RATf(trk_score_v); bdtv[15]=RATi(pfnplanehits_Y);
        bdtv[16]=RATf(trk_calo_energy_y_v); bdtv[17]=RATf(trk_pid_chipr_v);
        bdtv[18]=pfnplanehits_Y->at(shower1_index); bdtv[19]=pfpDist[shower1_index];
        bdtv[20]=shr_tkfit_dedx_y_v->at(shower1_index); bdtv[21]=(shrE>0?shrE/1000.0:-1.0);
        bdtv[22]=dMuPi; bdtv[23]=dPrPi; bdtv[24]=dPrMu; bdtv[25]=dS1Mu; bdtv[26]=dS1Pi; bdtv[27]=dS1Pr;
        bdtv[28]=muonDir.Angle(pionDir); bdtv[29]=protonDir.Angle(pionDir); bdtv[30]=protonDir.Angle(muonDir);
        bdtv[31]=shower1Dir.Angle(muonDir); bdtv[32]=shower1Dir.Angle(pionDir); bdtv[33]=shower1Dir.Angle(protonDir);
        bdtv[34]=mLambda; bdtv[35]=mSigma; bdtv[36]=topological_score; bdtv[37]=ntracks; bdtv[38]=nshowers;
        bdtv[39]=lamPointAngle; bdtv[40]=lamImpact; bdtv[41]=lamDecayDist;
        bdtv[42]=phoPointAngle; bdtv[43]=phoImpact; bdtv[44]=lamGamAngle; bdtv[45]=sigMuAngle;
        bdtv[46]=pi0Mass; bdtv[47]=pi0Cos; bdtv[48]=shr2E;
        // Shower Moliere angle, in DEGREES (range ~0-180). The branch carries
        // FLT_MAX (3.4e38) when undefined -- 0.03% of background events -- which
        // would swamp any BDT split, so anything outside the physical range is
        // mapped onto the same -1 sentinel used for a missing shower.
        {
          float mol = (shower1_index > -1 && shr_moliere_avg_v)
                      ? shr_moliere_avg_v->at(shower1_index) : -1.0f;
          if ( !std::isfinite(mol) || mol < 0.f || mol > 180.f ) mol = -1.0f;
          bdtv[49] = mol;
        }
        // --- Second-Shower Veto block (slice unassociated-hit variables excluded,
        //     see the comment on BDTV above) ---
        bdtv[50] = std::isfinite(slclustfrac) ? slclustfrac : -1.0f;
        bdtv[51] = float(secondshower_Y_nhit);
        bdtv[52] = (secondshower_Y_dot >= -1.f && secondshower_Y_dot <= 1.f)
                   ? secondshower_Y_dot : -2.0f;
        bdtv[53] = float(shrsubclusters0 + shrsubclusters1 + shrsubclusters2);
        // slice hits belonging to no reconstructed PFP (see the BDTV comment)
        {
          long pfsum = 0;
          if ( pfnhits ) for ( auto h : *pfnhits ) pfsum += h;
          float una = ( slnhits > 0 ) ? float(slnhits - pfsum) : -1.0f;
          if ( una < 0.f && slnhits > 0 ) una = 0.f;
          bdtv[54] = una;
          bdtv[55] = ( slnhits > 0 && una >= 0.f ) ? una / float(slnhits) : -1.0f;
        }
        bdtv[56] = (proton_index>-1 && pion_index>-1) ? 1.0f : 0.0f;   // HasPair
        bdtv[57] = W;
        bdtIsOdd = (float)(((evt % 2) + 2) % 2);   // event-number parity (CV fold tag)
        for (int k = 0; k < NBV - 1; k++) if (!std::isfinite(bdtv[k])) bdtv[k] = -1.0f;
        // Score with the fold that did NOT train on this event: even→B, odd→A.
        if (bdtReady) bdtg = reader->EvaluateMVA(bdtIsOdd > 0.5 ? "BDTG_A" : "BDTG_B");
        if (writeBDT && s != 4 && !ovlSig) {
          bool hasPi0 = false; if (mc_pdg) for (int p : *mc_pdg) if (p == 111) { hasPi0 = true; break; }
          tScore  = bdtg;            tSample = s;
          tNuPdg  = nu_pdg;          tCcnc   = ccnc;       tInter = interaction;
          tNuE    = nu_e;            tHasPi0 = hasPi0 ? 1 : 0;
          if (fd.isData) {               // EXT: no backtracking, and all of it is cosmic
            tMuPdg = tPrPdg = tPiPdg = tShrPdg = 0; tCosmic = 1;
          } else {
            tMuPdg  = backtracked_pdg->at(muon_index);
            tPrPdg  = (proton_index>-1) ? backtracked_pdg->at(proton_index) : 0;
            tPiPdg  = (pion_index>-1) ? backtracked_pdg->at(pion_index) : 0;
            tShrPdg = backtracked_pdg->at(shower1_index);
            tCosmic = (backtracked_overlay_purity->at(muon_index) > 0.5) ? 1 : 0;
          }

          // ── all-PFP block ───────────────────────────────────────────────────
          for (int k = 0; k < NSLOT; k++) for (int j = 0; j < NPV; j++) pfpv[k][j] = -1.0f;
          nPFPtot = (float)n;
          {
            std::vector<int> ord;
            for (int i = 0; i < n; i++) ord.push_back(i);
            auto hitsOf = [&](int i){
              return (pfnplanehits_Y && i < (int)pfnplanehits_Y->size())
                     ? pfnplanehits_Y->at(i) : 0; };
            std::sort(ord.begin(), ord.end(),
                      [&](int a, int b){ return hitsOf(a) > hitsOf(b); });
            auto GF = [&](std::vector<float>* v, int i){
              return (v && i < (int)v->size() && std::isfinite(v->at(i)))
                     ? (float)v->at(i) : -1.0f; };
            for (int k = 0; k < NSLOT && k < (int)ord.size(); k++) {
              const int i = ord[k];
              TVector3 d(trk_dir_x_v->at(i), trk_dir_y_v->at(i), trk_dir_z_v->at(i));
              const bool haveDir = d.Mag() > 0;
              if (haveDir) d = d.Unit();
              pfpv[k][0]  = GF(trk_score_v, i);
              pfpv[k][1]  = GF(trk_len_v, i);
              pfpv[k][2]  = GF(trk_llr_pid_score_v, i);
              pfpv[k][3]  = GF(trk_pid_chipr_v, i);
              pfpv[k][4]  = GF(trk_pid_chipi_v, i);
              pfpv[k][5]  = (float)pfpDist[i];
              pfpv[k][6]  = (pfnplanehits_Y && i < (int)pfnplanehits_Y->size())
                            ? (float)pfnplanehits_Y->at(i) : -1.0f;
              pfpv[k][7]  = GF(trk_calo_energy_y_v, i);
              pfpv[k][8]  = GF(trk_energy_proton_v, i);
              pfpv[k][9]  = GF(trk_range_muon_mom_v, i);
              pfpv[k][10] = GF(shr_tkfit_dedx_y_v, i);
              pfpv[k][11] = GF(shr_energy_y_v, i);
              pfpv[k][12] = GF(shr_moliere_avg_v, i);
              pfpv[k][13] = (haveDir && muonDir.Mag() > 0) ? (float)d.Dot(muonDir) : -2.0f;
              pfpv[k][14] = haveDir ? (float)d.Z() : -2.0f;
            }
            // pfpDist can be non-finite when a start position is unset, and TMVA
            // aborts outright on NaN/inf rather than ignoring the event. Sweep the
            // whole block rather than guarding each entry.
            for (int k = 0; k < NSLOT; k++)
              for (int j = 0; j < NPV; j++)
                if (!std::isfinite(pfpv[k][j])) pfpv[k][j] = -1.0f;
          }
          if (s == 0) TreeS->Fill(); else TreeB->Fill();
        }
      }
      bool sel_bdt = bdtReady ? (bdtg > BDT_MIN) : true;   // pass-through if no weights
      Cuts[8] = Cuts[7] && sel_bdt;                 // BDTG

      int cmax = -1;
      for (int c = 0; c < Nc; c++) { if (Cuts[c]) cmax = c; else break; }
      if (ovlSig) { for (int c = 0; c <= cmax; c++) { OvlSig[c] += W; OvlSigW2[c] += W * W; OvlSigRaw[c]++; } continue; }
      if (dvT) {
        dvRun = run; dvSub = sub; dvEvt = evt; dvNuE = nu_e; dvSig = (s == 0);
        dvCmax = cmax; dvBdtg = bdtg; dvW = W; dvT->Fill();
      }
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
  std::cout << "Skipped ragged entries: " << nSkippedRagged
            << ", duplicate data events: " << nDupData
            << ", overlay signal events removed from the background: " << nOverlaySig << "\n";

  // ── GENIE signal: the overlay's own signal events (not in any sample) ───────
  std::cout << "\n========== GENIE overlay signal (same definition, not in B) ==========\n";
  std::cout << std::left << std::setw(18) << "Cut" << std::right << std::setw(10) << "Raw"
            << std::setw(12) << "S_GENIE" << std::setw(9) << "Eff[%]" << std::setw(13) << "GENIE/NuWro" << "\n";
  for (int c = 0; c < Nc; c++)
    std::cout << std::left << std::setw(18) << CutName[c] << std::right
              << std::setw(10) << OvlSigRaw[c]
              << std::setw(12) << std::fixed << std::setprecision(3) << OvlSig[c]
              << std::setw(9) << std::setprecision(2) << (OvlSig[0] > 0 ? 100.0 * OvlSig[c] / OvlSig[0] : 0.0)
              << std::setw(13) << std::setprecision(3) << (Selected[c][0] > 0 ? OvlSig[c] / Selected[c][0] : 0.0) << "\n";
  std::cout << "======================================================================\n";
  fout->cd();
  TH1D* hGenie = new TH1D("GenieSignal", "GENIE overlay signal per cut (not in B);cut;events", Nc, 0, Nc);
  for (int c = 0; c < Nc; c++) {
    hGenie->GetXaxis()->SetBinLabel(c + 1, CutName[c]);
    hGenie->SetBinContent(c + 1, OvlSig[c]);
    hGenie->SetBinError(c + 1, std::sqrt(OvlSigW2[c]));
  }

  // ── Signal-efficiency / significance diagnostics ────────────────────────────
  //  eff   = raw signal events at cut c / raw signal at AllEvents (true Sig in FV)
  //  S,B   = POT-weighted signal / (hyperons+background+dirt+EXT) yields
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
    const double B = Selected[c][1] + Selected[c][2] + Selected[c][3] + Selected[c][5];
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
  if (frw.ok)
    std::cout << "NuMI flux reweight applied to " << frw.nUsed << " NuWro events ("
              << frw.nEmpty << " outside the map -> weight 1)\n";
  std::cout << "Wrote " << outname << "\n";

  if (writeBDT) {
    fbdt->cd();
    std::cout << "BDT training trees: TreeS=" << TreeS->GetEntries()
              << " signal, TreeB=" << TreeB->GetEntries() << " background\n";
    TreeS->Write(); TreeB->Write();
    fbdt->Close();
    std::cout << "Wrote " << (requirePion ? "BDTInput.root" : "BDTInput_nopion.root") << "\n";
  }
  delete reader;
}
