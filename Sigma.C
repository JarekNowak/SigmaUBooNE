#define Sigma_cxx
#include "Sigma.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include "SigmaConfig.h"
#include "TMVA/Reader.h"

using namespace std;
using namespace SigmaCfg;   // Ns, Nc, Nv, Names[], CutName[], VarName[], Binning()

// ============================================================
//  Helper: safely retrieve from a vector, returning a sentinel
//  if the index is out of range or the pointer is null.
// ============================================================
template<typename T>
static inline float SafeGet(const std::vector<T>* v, int idx, float sentinel = -999.f)
{
  if (!v || idx < 0 || idx >= (int)v->size()) return sentinel;
  return (float)v->at(idx);
}

template<typename T>
static inline int SafeGetI(const std::vector<T>* v, int idx, int sentinel = -999)
{
  if (!v || idx < 0 || idx >= (int)v->size()) return sentinel;
  return (int)v->at(idx);
}

// ============================================================
void Sigma::Loop()
// ============================================================
{
  gStyle->SetOptStat(1001111);

  // -----------------------------------------------------------
  // Sample labels and bookkeeping
  // -----------------------------------------------------------
  // Ns, Names[], Nc, CutName[], Nv, VarName[] and the per-variable Binning()
  // all come from SigmaConfig.h (shared with the plotting macro).
  //
  // Sample is assigned per event from the tree index within the chain
  // (see the main loop below), so no hardcoded entry boundaries are needed.

  // Protons on Target — used for luminosity scaling to run-3b data
  double PoT[Ns] = {0};
  PoT[0] = 552.0;   // *1e20  Sigma MC
  PoT[1] = PoT[0];  //        Hyperons MC  (same sample)
  PoT[2] = 9.280;   // *1e20  Background MC   (FIX: was swapped with Dirt)
  PoT[3] = 17.30;   // *1e20  Dirt MC
  PoT[4] = 4.90;    // *1e20  Data (run 3b)

  double scale[Ns] = {0};
  for (int i = 0; i < Ns - 1; i++) scale[i] = PoT[4] / PoT[i];
  scale[4] = 1.0;

  for (int i = 0; i < Ns; i++) cout << Names[i] << "  scale=" << scale[i] << endl;

  // -----------------------------------------------------------
  // Histogram booking — explicit per-variable binning from SigmaConfig.h
  // (shared with the plotting macro, so all samples share identical axes
  //  and can be stacked/overlaid).
  // -----------------------------------------------------------
  TH1F *Hist[Nv][Ns][Nc];
  for (int v = 0; v < Nv; v++) {
    int nb; double lo, hi; bool logy;
    Binning(v, nb, lo, hi, logy);
    for (int s = 0; s < Ns; s++)
      for (int c = 0; c < Nc; c++)
        Hist[v][s][c] = new TH1F(VarName[v] + Names[s] + CutName[c], "", nb, lo, hi);
  }

  // -----------------------------------------------------------
  // TMVA output file and trees
  // -----------------------------------------------------------
  TFile *tmvaFile = new TFile("TMVA_input.root", "RECREATE");

  TTree *TreeS = new TTree("TreeS", "Signal tree for TMVA (Sigma events after all cuts)");
  TTree *TreeB = new TTree("TreeB", "Background tree for TMVA (Bkg+Dirt after all cuts)");

  // --- TMVA branch variables (Float_t required by TMVA) ---
  // Track-level PID and kinematics
  Float_t t_MuonPID,        t_MuonTrkLength,    t_MuonTrkDistVtx,  t_MuonTrkScore;
  Float_t t_MuonDeDxP0,     t_MuonDeDxP1,       t_MuonDeDxP2;
  Float_t t_MuonNHitsY,     t_MuonPhotonPID;

  Float_t t_PionPID,        t_PionTrkLength,     t_PionTrkDistVtx,  t_PionTrkScore;
  Float_t t_PionDeDxP0,     t_PionDeDxP1,        t_PionDeDxP2;
  Float_t t_PionNHitsY,     t_PionPhotonPID;

  Float_t t_ProtonPID,      t_ProtonTrkLength,   t_ProtonTrkDistVtx,t_ProtonTrkScore;
  Float_t t_ProtonDeDxP0,   t_ProtonDeDxP1,      t_ProtonDeDxP2;
  Float_t t_ProtonNHitsY,   t_ProtonPhotonPID;

  // Shower variables
  Float_t t_Shr1NHitsY,     t_Shr1VtxDist,      t_Shr1DeDxP0,      t_Shr1DeDxP1,    t_Shr1DeDxP2;
  Float_t t_Shr1OpenAngle;
  Float_t t_Shr2NHitsY,     t_Shr2VtxDist,      t_Shr2DeDxP0,      t_Shr2DeDxP1,    t_Shr2DeDxP2;
  Float_t t_Shr2OpenAngle;

  // Opening angles between candidates
  Float_t t_AngleMuonPion,  t_AngleProtonPion,   t_AngleProtonMuon;
  Float_t t_AngleShr1Muon,  t_AngleShr1Pion,     t_AngleShr1Proton;
  Float_t t_AngleShr2Muon,  t_AngleShr2Pion,     t_AngleShr2Proton;

  // Distances between candidates
  Float_t t_DistMuonPion,   t_DistProtonPion,    t_DistProtonMuon;
  Float_t t_DistShr1Muon,   t_DistShr1Pion,      t_DistShr1Proton;
  Float_t t_DistShr2Muon,   t_DistShr2Pion,      t_DistShr2Proton;

  // Topological
  Float_t t_NTracks,        t_NShowers;

  // Event weight for TMVA
  Float_t t_Weight;

  // Helper lambda to add identical branches to both trees
  auto AddBranch = [&](const char* name, Float_t& var) {
    TreeS->Branch(name, &var, (std::string(name) + "/F").c_str());
    TreeB->Branch(name, &var, (std::string(name) + "/F").c_str());
  };

  AddBranch("MuonPID",           t_MuonPID);
  AddBranch("MuonTrkLength",     t_MuonTrkLength);
  AddBranch("MuonTrkDistVtx",    t_MuonTrkDistVtx);
  AddBranch("MuonTrkScore",      t_MuonTrkScore);
  AddBranch("MuonDeDxPlane0",    t_MuonDeDxP0);
  AddBranch("MuonDeDxPlane1",    t_MuonDeDxP1);
  AddBranch("MuonDeDxPlane2",    t_MuonDeDxP2);
  AddBranch("MuonNHitsY",        t_MuonNHitsY);
  AddBranch("MuonPhotonPID",     t_MuonPhotonPID);

  AddBranch("PionPID",           t_PionPID);
  AddBranch("PionTrkLength",     t_PionTrkLength);
  AddBranch("PionTrkDistVtx",    t_PionTrkDistVtx);
  AddBranch("PionTrkScore",      t_PionTrkScore);
  AddBranch("PionDeDxPlane0",    t_PionDeDxP0);
  AddBranch("PionDeDxPlane1",    t_PionDeDxP1);
  AddBranch("PionDeDxPlane2",    t_PionDeDxP2);
  AddBranch("PionNHitsY",        t_PionNHitsY);
  AddBranch("PionPhotonPID",     t_PionPhotonPID);

  AddBranch("ProtonPID",         t_ProtonPID);
  AddBranch("ProtonTrkLength",   t_ProtonTrkLength);
  AddBranch("ProtonTrkDistVtx",  t_ProtonTrkDistVtx);
  AddBranch("ProtonTrkScore",    t_ProtonTrkScore);
  AddBranch("ProtonDeDxPlane0",  t_ProtonDeDxP0);
  AddBranch("ProtonDeDxPlane1",  t_ProtonDeDxP1);
  AddBranch("ProtonDeDxPlane2",  t_ProtonDeDxP2);
  AddBranch("ProtonNHitsY",      t_ProtonNHitsY);
  AddBranch("ProtonPhotonPID",   t_ProtonPhotonPID);

  AddBranch("Shr1NHitsY",        t_Shr1NHitsY);
  AddBranch("Shr1VtxDist",       t_Shr1VtxDist);
  AddBranch("Shr1DeDxPlane0",    t_Shr1DeDxP0);
  AddBranch("Shr1DeDxPlane1",    t_Shr1DeDxP1);
  AddBranch("Shr1DeDxPlane2",    t_Shr1DeDxP2);
  AddBranch("Shr1OpenAngle",     t_Shr1OpenAngle);
  AddBranch("Shr2NHitsY",        t_Shr2NHitsY);
  AddBranch("Shr2VtxDist",       t_Shr2VtxDist);
  AddBranch("Shr2DeDxPlane0",    t_Shr2DeDxP0);
  AddBranch("Shr2DeDxPlane1",    t_Shr2DeDxP1);
  AddBranch("Shr2DeDxPlane2",    t_Shr2DeDxP2);
  AddBranch("Shr2OpenAngle",     t_Shr2OpenAngle);

  AddBranch("AngleMuonPion",     t_AngleMuonPion);
  AddBranch("AngleProtonPion",   t_AngleProtonPion);
  AddBranch("AngleProtonMuon",   t_AngleProtonMuon);
  AddBranch("AngleShr1Muon",     t_AngleShr1Muon);
  AddBranch("AngleShr1Pion",     t_AngleShr1Pion);
  AddBranch("AngleShr1Proton",   t_AngleShr1Proton);
  AddBranch("AngleShr2Muon",     t_AngleShr2Muon);
  AddBranch("AngleShr2Pion",     t_AngleShr2Pion);
  AddBranch("AngleShr2Proton",   t_AngleShr2Proton);

  AddBranch("DistMuonPion",      t_DistMuonPion);
  AddBranch("DistProtonPion",    t_DistProtonPion);
  AddBranch("DistProtonMuon",    t_DistProtonMuon);
  AddBranch("DistShr1Muon",      t_DistShr1Muon);
  AddBranch("DistShr1Pion",      t_DistShr1Pion);
  AddBranch("DistShr1Proton",    t_DistShr1Proton);
  AddBranch("DistShr2Muon",      t_DistShr2Muon);
  AddBranch("DistShr2Pion",      t_DistShr2Pion);
  AddBranch("DistShr2Proton",    t_DistShr2Proton);

  AddBranch("NTracks",           t_NTracks);
  AddBranch("NShowers",          t_NShowers);
  AddBranch("Weight",            t_Weight);

  // -----------------------------------------------------------
  // TMVA Reader — apply trained classifiers as the last selection cuts
  // -----------------------------------------------------------
  Float_t r_MuonTrkLength,    r_MuonTrkDistVtx,  r_MuonTrkScore;
  Float_t r_MuonNHitsY,       r_MuonPhotonPID;
  Float_t r_PionPID,          r_PionTrkLength,   r_PionTrkDistVtx,  r_PionTrkScore;
  Float_t r_PionNHitsY,       r_PionDeDxP2,      r_PionPhotonPID;
  Float_t r_ProtonPID,        r_ProtonTrkLength, r_ProtonTrkDistVtx,r_ProtonTrkScore;
  Float_t r_ProtonNHitsY,     r_ProtonDeDxP2,    r_ProtonPhotonPID;
  Float_t r_Shr1NHitsY,       r_Shr1VtxDist,     r_Shr1OpenAngle;
  Float_t r_Shr2NHitsY,       r_Shr2VtxDist;
  Float_t r_AngleMuonPion,    r_AngleProtonPion,  r_AngleProtonMuon;
  Float_t r_AngleShr1Muon,    r_AngleShr1Pion,    r_AngleShr1Proton;
  Float_t r_DistMuonPion,     r_DistProtonPion,   r_DistProtonMuon;
  Float_t r_Weight;

  TMVA::Reader *mvaReader = new TMVA::Reader("!Color:!Silent");

  mvaReader->AddVariable("MuonTrkLength",    &r_MuonTrkLength);
  mvaReader->AddVariable("MuonTrkDistVtx",   &r_MuonTrkDistVtx);
  mvaReader->AddVariable("MuonTrkScore",     &r_MuonTrkScore);
  mvaReader->AddVariable("MuonNHitsY",       &r_MuonNHitsY);
  mvaReader->AddVariable("MuonPhotonPID",    &r_MuonPhotonPID);
  mvaReader->AddVariable("PionPID",          &r_PionPID);
  mvaReader->AddVariable("PionTrkLength",    &r_PionTrkLength);
  mvaReader->AddVariable("PionTrkDistVtx",   &r_PionTrkDistVtx);
  mvaReader->AddVariable("PionTrkScore",     &r_PionTrkScore);
  mvaReader->AddVariable("PionNHitsY",       &r_PionNHitsY);
  mvaReader->AddVariable("PionDeDxPlane2",   &r_PionDeDxP2);
  mvaReader->AddVariable("PionPhotonPID",    &r_PionPhotonPID);
  mvaReader->AddVariable("ProtonPID",        &r_ProtonPID);
  mvaReader->AddVariable("ProtonTrkLength",  &r_ProtonTrkLength);
  mvaReader->AddVariable("ProtonTrkDistVtx", &r_ProtonTrkDistVtx);
  mvaReader->AddVariable("ProtonTrkScore",   &r_ProtonTrkScore);
  mvaReader->AddVariable("ProtonNHitsY",     &r_ProtonNHitsY);
  mvaReader->AddVariable("ProtonDeDxPlane2", &r_ProtonDeDxP2);
  mvaReader->AddVariable("ProtonPhotonPID",  &r_ProtonPhotonPID);
  mvaReader->AddVariable("Shr1NHitsY",       &r_Shr1NHitsY);
  mvaReader->AddVariable("Shr1VtxDist",      &r_Shr1VtxDist);
  mvaReader->AddVariable("Shr1OpenAngle",    &r_Shr1OpenAngle);
  mvaReader->AddVariable("Shr2NHitsY",       &r_Shr2NHitsY);
  mvaReader->AddVariable("Shr2VtxDist",      &r_Shr2VtxDist);
  mvaReader->AddVariable("AngleMuonPion",    &r_AngleMuonPion);
  mvaReader->AddVariable("AngleProtonPion",  &r_AngleProtonPion);
  mvaReader->AddVariable("AngleProtonMuon",  &r_AngleProtonMuon);
  mvaReader->AddVariable("AngleShr1Muon",    &r_AngleShr1Muon);
  mvaReader->AddVariable("AngleShr1Pion",    &r_AngleShr1Pion);
  mvaReader->AddVariable("AngleShr1Proton",  &r_AngleShr1Proton);
  mvaReader->AddVariable("DistMuonPion",     &r_DistMuonPion);
  mvaReader->AddVariable("DistProtonPion",   &r_DistProtonPion);
  mvaReader->AddVariable("DistProtonMuon",   &r_DistProtonMuon);

  mvaReader->AddSpectator("Weight", &r_Weight);

  mvaReader->BookMVA("BDT",    "SigmaDataset/weights/TMVASigma_BDT.weights.xml");
  mvaReader->BookMVA("BDTAda", "SigmaDataset/weights/TMVASigma_BDTAda.weights.xml");
  mvaReader->BookMVA("BDTg",   "SigmaDataset/weights/TMVASigma_BDTg.weights.xml");

  // -----------------------------------------------------------
  // Guard against empty chain
  // -----------------------------------------------------------
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();   // true total (GetEntriesFast() returns a bogus sentinel for a TChain)

  // -----------------------------------------------------------
  // Main event loop
  // -----------------------------------------------------------
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    // ---- Progress status line (updated in place) ----
    if (jentry % 10000 == 0 || jentry == nentries - 1)
      cout << "\rProcessing entry " /*<< (jentry + 1) << " / " << nentries << "  ("*/ << (100 * (jentry + 1) / nentries) <<"%" << flush;

    // ---- Reset per-event variables ----
    bool IsSignal         = false;
    int  mc_decay         = 0;
    int  s                = -1;
    int  muon_index       = -1;
    int  pion_index       = -1;
    int  proton_index     = -1;
    int  shower1_index    = -1;
    int  shower2_index    = -1;
    bool Cuts[Nc]         = {false};

    bool sel_FV              = false;
    bool sel_ntracks         = false;
    bool sel_nshowers        = false;
    bool sel_muon_candidate  = false;
    bool sel_pion_candidate  = false;
    bool sel_proton_candidate= false;

      double proton_mom=0;
      double pion_mom=0;
    // ---- Truth-level decay product counting ----
    for (int i = 0; i < (int)mc_decay_pdg->size(); i++) {
        if (mc_decay_pdg->at(i) ==  2212) {mc_decay++; proton_mom = SafeGet(mc_decay_mom, i);}  // proton
        if (mc_decay_pdg->at(i) == -211)  {mc_decay++; pion_mom   = SafeGet(mc_decay_mom, i);} // pi-
    }

    // ---- Truth fiducial-volume flags (applied per-sample below) ----
    bool inFV_truth = (mc_nu_pos_x >   3.00  && mc_nu_pos_x < 253.35 &&
                       mc_nu_pos_y > -112.53  && mc_nu_pos_y < 114.47 &&
                       mc_nu_pos_z >   3.1    && mc_nu_pos_z < 1033.0);
    bool inDeadRegion = (mc_nu_pos_z > 675.1  && mc_nu_pos_z < 775.1);

    // ---- Signal definition: anti-nu_mu + Sigma0 with proton + pi- decay products ----
    if ((mc_nu_pdg == -14) && (mc_hyperon_pdg == 3212) && (mc_decay == 2)
        && proton_mom>0.3 && pion_mom>0.1) IsSignal = true;

    // ---- Sample assignment from the tree index within the chain ----
    //   tree 0 = Hyperon/Sigma file, 1 = Background, 2 = Dirt (order of chain->Add in Sigma.h)
    int treenum = fChain->GetTreeNumber();
    if      (treenum == 0) s = IsSignal ? 0 : 1;  // Sigma signal / other hyperons
    else if (treenum == 1) s = 2;                 // Background
    else if (treenum == 2) s = 3;                 // Dirt
    // Data (s=4) not yet available per original code

    if (s < 0) continue;  // Safety: skip unclassified events

    // ---- Truth FV cut applies ONLY to samples with an in-FV true vertex by design
    //      (signal & other hyperons). Background uses overlay neutrinos and Dirt is
    //      out-of-FV by construction, so both rely on the reco FV cut (Cuts[1]). ----
    if ((s == 0 || s == 1) && (!inFV_truth || inDeadRegion)) continue;

    // ---- Reconstructed vertex ----
    TVector3 vtx(reco_primary_vtx_x, reco_primary_vtx_y, reco_primary_vtx_z);
    if (reco_primary_vtx_inFV) sel_FV = true;

    // ---- Count tracks and showers; check for large showers ----
    int  ntracks      = 0;
    int  nshowers     = 0;
    bool NoBigShowers = true;

    for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
      if (pfp_trk_shr_score->at(i) < 0.5) {
        nshowers++;
        if (pfpnplanehits_Y->at(i) > 600) NoBigShowers = false;
      } else {
        ntracks++;
      }
    }

    // ---- Muon candidate: longest track, PID > 0.7, within 20 cm of vertex ----
    for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
      TVector3 start(trk_start_x->at(i), trk_start_y->at(i), trk_start_z->at(i));
      double dist = (start - vtx).Mag();
      if (pfp_trk_shr_score->at(i) > 0.8 && dist < 20 &&
          trk_length->at(i) > 10          && trk_llrpid->at(i) > 0.7) {
        sel_muon_candidate = true;
        if (muon_index == -1 || trk_length->at(i) > trk_length->at(muon_index))
          muon_index = i;
      }
    }

    // ---- Pion candidate: highest LLR PID (excluding muon) ----
    for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
      if (i == muon_index) continue;
      if (pfp_trk_shr_score->at(i) > 0.5) {
        sel_pion_candidate = true;
        if (pion_index == -1 || trk_llrpid->at(i) > trk_llrpid->at(pion_index))
          pion_index = i;
      }
    }

    // ---- Proton candidate: lowest LLR PID (excluding muon + pion) ----
    for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
      if (i == muon_index || i == pion_index) continue;
      if (pfp_trk_shr_score->at(i) > 0.5)  {
        sel_proton_candidate = true;
        if (proton_index == -1 || trk_llrpid->at(i) < trk_llrpid->at(proton_index))
          proton_index = i;
      }
    }

    // ---- Shower candidates: the TWO showers closest to the vertex ----
    double shr_best1 = 1e9, shr_best2 = 1e9;
    for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
      if (i == muon_index || i == pion_index || i == proton_index) continue;
      if (pfp_trk_shr_score->at(i) >= 0.5) continue;  // showers only
      TVector3 start(trk_start_x->at(i), trk_start_y->at(i), trk_start_z->at(i));
      double dist = (start - vtx).Mag();
      if (dist < shr_best1) {
        shr_best2 = shr_best1;  shower2_index = shower1_index;
        shr_best1 = dist;       shower1_index = i;
      } else if (dist < shr_best2) {
        shr_best2 = dist;       shower2_index = i;
      }
    }

    // ---- Compute inter-candidate positions & distances (computed ONCE per event) ----
    TVector3 muonStart,  pionStart,   protonStart;
    TVector3 shower1Start, shower2Start;
    TVector3 muonDir,    pionDir,     protonDir;
    TVector3 shower1Dir, shower2Dir;

    if (muon_index    > -1) {
      muonStart.SetXYZ(trk_start_x->at(muon_index),    trk_start_y->at(muon_index),    trk_start_z->at(muon_index));
      muonDir.SetXYZ  (trk_dir_x->at(muon_index),      trk_dir_y->at(muon_index),      trk_dir_z->at(muon_index));
    }
    if (pion_index    > -1) {
      pionStart.SetXYZ(trk_start_x->at(pion_index),    trk_start_y->at(pion_index),    trk_start_z->at(pion_index));
      pionDir.SetXYZ  (trk_dir_x->at(pion_index),      trk_dir_y->at(pion_index),      trk_dir_z->at(pion_index));
    }
    if (proton_index  > -1) {
      protonStart.SetXYZ(trk_start_x->at(proton_index),trk_start_y->at(proton_index),  trk_start_z->at(proton_index));
      protonDir.SetXYZ  (trk_dir_x->at(proton_index),  trk_dir_y->at(proton_index),    trk_dir_z->at(proton_index));
    }
    if (shower1_index > -1) {
      shower1Start.SetXYZ(trk_start_x->at(shower1_index),trk_start_y->at(shower1_index),trk_start_z->at(shower1_index));
      shower1Dir.SetXYZ  (trk_dir_x->at(shower1_index),  trk_dir_y->at(shower1_index),  trk_dir_z->at(shower1_index));
    }
    if (shower2_index > -1) {
      shower2Start.SetXYZ(trk_start_x->at(shower2_index),trk_start_y->at(shower2_index),trk_start_z->at(shower2_index));
      shower2Dir.SetXYZ  (trk_dir_x->at(shower2_index),  trk_dir_y->at(shower2_index),  trk_dir_z->at(shower2_index));
    }

    double distMuonPion    = (muon_index>-1  && pion_index>-1)   ? (muonStart  -pionStart).Mag()   : 99999;
    double distProtonPion  = (proton_index>-1&& pion_index>-1)   ? (protonStart-pionStart).Mag()   : 99999;
    double distProtonMuon  = (proton_index>-1&& muon_index>-1)   ? (protonStart-muonStart).Mag()   : 99999;
    double distShower1Muon = (shower1_index>-1&&muon_index>-1)   ? (shower1Start-muonStart).Mag()  : 99999;
    double distShower1Pion = (shower1_index>-1&&pion_index>-1)   ? (shower1Start-pionStart).Mag()  : 99999;
    double distShower1Prot = (shower1_index>-1&&proton_index>-1) ? (shower1Start-protonStart).Mag(): 99999;
    double distShower2Muon = (shower2_index>-1&&muon_index>-1)   ? (shower2Start-muonStart).Mag()  : 99999;
    double distShower2Pion = (shower2_index>-1&&pion_index>-1)   ? (shower2Start-pionStart).Mag()  : 99999;
    double distShower2Prot = (shower2_index>-1&&proton_index>-1) ? (shower2Start-protonStart).Mag(): 99999;

    // ---- Boolean flags for shower quality ----
    bool sel_more_showers = true;
    if (shower2_index > -1 && pfpnplanehits_Y->at(shower2_index) > 100) sel_more_showers = false;

    bool sel_1st_shower = false;
    if (shower1_index > -1) {
      int hits1 = pfpnplanehits_Y->at(shower1_index);
      if (hits1 > 10 && hits1 < 300) sel_1st_shower = true;
    }

    if (nshowers > 0 && nshowers < 4) sel_nshowers = true;  // loosened: allow up to 3 showers
    if (ntracks < 5)                  sel_ntracks   = true;

    // ---- FIX: corrected bitwise & -> && in distance cut ----
    bool sel_distances =
      (distMuonPion < 30 && distProtonMuon < 20 && distProtonPion < 20 &&
       distShower1Muon < 100 && distShower1Pion < 100 && distShower1Prot < 100);

    // ---- Build cumulative cut flags ----
    Cuts[0]  = true;
    Cuts[1]  = Cuts[0]  && sel_FV;
    Cuts[2]  = Cuts[1]  && sel_nshowers;
    Cuts[3]  = Cuts[2];  // NoBigShowers removed: fully redundant after BDTg cut
    Cuts[4]  = Cuts[3]  && sel_muon_candidate;
    Cuts[5]  = Cuts[4]  && sel_pion_candidate;
    Cuts[6]  = Cuts[5]  && sel_proton_candidate;
    Cuts[7]  = Cuts[6]  && sel_ntracks;
    Cuts[8]  = Cuts[7];  // 2nd shower veto removed: recovers +2% signal with same BDTg purity
    Cuts[9]  = Cuts[8]  && sel_1st_shower;
    Cuts[10] = Cuts[9]  && sel_distances;

    // ---- MVA classifier evaluation ----
    float mva_BDT = -999.f, mva_BDTAda = -999.f, mva_BDTg = -999.f;
    if (Cuts[10]) {
      r_MuonTrkLength  = (muon_index>-1)    ? SafeGet(trk_length,           muon_index)             : -999.f;
      r_MuonTrkDistVtx = (muon_index>-1)    ? (Float_t)(muonStart  - vtx).Mag()                    : -999.f;
      r_MuonTrkScore   = (muon_index>-1)    ? SafeGet(pfp_trk_shr_score,   muon_index)             : -999.f;
      r_MuonNHitsY     = (muon_index>-1)    ? (Float_t)SafeGetI(pfpnplanehits_Y, muon_index)       : -999.f;
      r_MuonPhotonPID  = (muon_index>-1)    ? SafeGet(pfp_photon_bdt_score, muon_index)            : -999.f;

      r_PionPID        = (pion_index>-1)    ? SafeGet(trk_llrpid,             pion_index)            : -999.f;
      r_PionTrkLength  = (pion_index>-1)    ? SafeGet(trk_length,              pion_index)            : -999.f;
      r_PionTrkDistVtx = (pion_index>-1)    ? (Float_t)(pionStart  - vtx).Mag()                     : -999.f;
      r_PionTrkScore   = (pion_index>-1)    ? SafeGet(pfp_trk_shr_score,      pion_index)            : -999.f;
      r_PionNHitsY     = (pion_index>-1)    ? (Float_t)SafeGetI(pfpnplanehits_Y, pion_index)         : -999.f;
      r_PionDeDxP2     = (pion_index>-1)    ? SafeGet(trk_mean_dedx_plane2,   pion_index)            : -999.f;
      r_PionPhotonPID  = (pion_index>-1)    ? SafeGet(pfp_photon_bdt_score,    pion_index)            : -999.f;

      r_ProtonPID        = (proton_index>-1) ? SafeGet(trk_llrpid,              proton_index)         : -999.f;
      r_ProtonTrkLength  = (proton_index>-1) ? SafeGet(trk_length,               proton_index)         : -999.f;
      r_ProtonTrkDistVtx = (proton_index>-1) ? (Float_t)(protonStart - vtx).Mag()                    : -999.f;
      r_ProtonTrkScore   = (proton_index>-1) ? SafeGet(pfp_trk_shr_score,       proton_index)         : -999.f;
      r_ProtonNHitsY     = (proton_index>-1) ? (Float_t)SafeGetI(pfpnplanehits_Y, proton_index)       : -999.f;
      r_ProtonDeDxP2     = (proton_index>-1) ? SafeGet(trk_mean_dedx_plane2,    proton_index)         : -999.f;
      r_ProtonPhotonPID  = (proton_index>-1) ? SafeGet(pfp_photon_bdt_score,     proton_index)         : -999.f;

      r_Shr1NHitsY    = (shower1_index>-1) ? (Float_t)SafeGetI(pfpnplanehits_Y, shower1_index)      : -999.f;
      r_Shr1VtxDist   = (shower1_index>-1) ? (Float_t)(shower1Start - vtx).Mag()                   : -999.f;
      r_Shr1OpenAngle = (shower1_index>-1) ? SafeGet(shr_open_angle, shower1_index)                 : -999.f;
      r_Shr2NHitsY    = (shower2_index>-1) ? (Float_t)SafeGetI(pfpnplanehits_Y, shower2_index)      : -999.f;
      r_Shr2VtxDist   = (shower2_index>-1) ? (Float_t)(shower2Start - vtx).Mag()                   : -999.f;

      r_AngleMuonPion    = (muon_index>-1   && pion_index>-1)    ? (Float_t)muonDir.Angle(pionDir)     : -999.f;
      r_AngleProtonPion  = (proton_index>-1 && pion_index>-1)    ? (Float_t)protonDir.Angle(pionDir)   : -999.f;
      r_AngleProtonMuon  = (proton_index>-1 && muon_index>-1)    ? (Float_t)protonDir.Angle(muonDir)   : -999.f;
      r_AngleShr1Muon    = (shower1_index>-1&& muon_index>-1)    ? (Float_t)shower1Dir.Angle(muonDir)  : -999.f;
      r_AngleShr1Pion    = (shower1_index>-1&& pion_index>-1)    ? (Float_t)shower1Dir.Angle(pionDir)  : -999.f;
      r_AngleShr1Proton  = (shower1_index>-1&& proton_index>-1)  ? (Float_t)shower1Dir.Angle(protonDir): -999.f;

      r_DistMuonPion   = (Float_t)distMuonPion;
      r_DistProtonPion = (Float_t)distProtonPion;
      r_DistProtonMuon = (Float_t)distProtonMuon;
      r_Weight         = (Float_t)scale[s];

      mva_BDT    = mvaReader->EvaluateMVA("BDT");
      mva_BDTAda = mvaReader->EvaluateMVA("BDTAda");
      mva_BDTg   = mvaReader->EvaluateMVA("BDTg");
    }

    Cuts[11] = Cuts[10];                            // BDT — covered by BDTg, pass-through
    Cuts[12] = Cuts[11];                            // BDTAda — covered by BDTg, pass-through
    Cuts[13] = Cuts[12] && (mva_BDTg   > 0.96f);

    // ======================================================
    // Fill histograms.  Cuts are cumulative, so the set of passing
    // levels is the prefix 0..cmax. Fill() writes each value into all
    // passing levels in one shot, so the PFP vector is traversed ONCE
    // per event instead of once per cut.
    // ======================================================
    int cmax = -1;
    for (int c = 0; c < Nc; c++) { if (Cuts[c]) cmax = c; else break; }

    auto Fill = [&](int v, double x) {
      for (int c = 0; c <= cmax; c++) Hist[v][s][c]->Fill(x, scale[s]);
    };

    if (cmax >= 0) {
      // Topological
      Fill(0, ntracks);
      Fill(1, nshowers);

      // Shower hit counts in Y plane
      if (nshowers > 0)
        for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++)
          if (pfp_trk_shr_score->at(i) < 0.5)
            Fill(2, pfpnplanehits_Y->at(i));

      // Opening angles and distances (pre-computed above)
      if (muon_index>-1  && pion_index>-1)    { Fill(58, muonDir.Angle(pionDir));      Fill(71, distMuonPion);    }
      if (proton_index>-1&& pion_index>-1)    { Fill(59, protonDir.Angle(pionDir));    Fill(72, distProtonPion);  }
      if (proton_index>-1&& muon_index>-1)    { Fill(60, protonDir.Angle(muonDir));    Fill(73, distProtonMuon);  }
      if (shower1_index>-1&&muon_index>-1)    { Fill(61, shower1Dir.Angle(muonDir));   Fill(74, distShower1Muon); }
      if (shower1_index>-1&&pion_index>-1)    { Fill(62, shower1Dir.Angle(pionDir));   Fill(75, distShower1Pion); }
      if (shower1_index>-1&&proton_index>-1)  { Fill(63, shower1Dir.Angle(protonDir)); Fill(76, distShower1Prot); }
      if (shower2_index>-1&&muon_index>-1)    { Fill(64, shower2Dir.Angle(muonDir));   Fill(77, distShower2Muon); }
      if (shower2_index>-1&&pion_index>-1)    { Fill(65, shower2Dir.Angle(pionDir));   Fill(78, distShower2Pion); }
      if (shower2_index>-1&&proton_index>-1)  { Fill(66, shower2Dir.Angle(protonDir)); Fill(79, distShower2Prot); }

      // Shower open angles and dEdx
      if (muon_index    > -1) Fill(47, SafeGet(shr_open_angle, muon_index));
      if (pion_index    > -1) Fill(48, SafeGet(shr_open_angle, pion_index));
      if (proton_index  > -1) Fill(49, SafeGet(shr_open_angle, proton_index));
      if (shower1_index > -1) {
        Fill(50, SafeGet(shr_open_angle,  shower1_index));
        Fill(52, SafeGet(shr_dedx_plane0, shower1_index));
        Fill(53, SafeGet(shr_dedx_plane1, shower1_index));
        Fill(54, SafeGet(shr_dedx_plane2, shower1_index));
        Fill(67, SafeGetI(pfpnplanehits_Y, shower1_index));
        Fill(69, (shower1Start - vtx).Mag());
      }
      if (shower2_index > -1) {
        Fill(51, SafeGet(shr_open_angle,  shower2_index));
        Fill(55, SafeGet(shr_dedx_plane0, shower2_index));
        Fill(56, SafeGet(shr_dedx_plane1, shower2_index));
        Fill(57, SafeGet(shr_dedx_plane2, shower2_index));
        Fill(68, SafeGetI(pfpnplanehits_Y, shower2_index));
        Fill(70, (shower2Start - vtx).Mag());
      }

      // Per-PFP particle-type histograms (species block laid out by index arithmetic)
      for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
        if (trk_length->at(i) <= 0) continue;
        TVector3 start(trk_start_x->at(i), trk_start_y->at(i), trk_start_z->at(i));
        double dist = (start - vtx).Mag();
        int pdg = pfp_true_pdg->at(i);

        int base = -1;
        if      (pdg == -13)  base = 0;   // muon
        else if (pdg == -211) base = 1;   // pion
        else if (pdg ==  2212)base = 2;   // proton
        else if (pdg ==  22)  base = 3;   // gamma
        if (base < 0) continue;

        int off = base * 4;  // stride for PID/Length/Dist/Score block (3,4,5,6)
        Fill(3  + off,    trk_llrpid->at(i));
        Fill(4  + off,    trk_length->at(i));
        Fill(5  + off,    dist);
        Fill(6  + off,    pfp_trk_shr_score->at(i));
        Fill(19 + base,   pfpnplanehits_Y->at(i));
        Fill(23 + base,   pfpnplanehits_U->at(i));
        Fill(27 + base,   pfpnplanehits_V->at(i));
        Fill(31 + base*3, trk_mean_dedx_plane0->at(i));
        Fill(32 + base*3, trk_mean_dedx_plane1->at(i));
        Fill(33 + base*3, trk_mean_dedx_plane2->at(i));
        Fill(43 + base,   pfp_photon_bdt_score->at(i));
      }
      // MVA scores (only valid for events passing Cuts[10])
      if (cmax >= 10) {
        Fill(80, mva_BDT);
        Fill(81, mva_BDTAda);
        Fill(82, mva_BDTg);
      }
    } // end fill block

    // ======================================================
    // TMVA tree filling — events passing ALL cuts (Cuts[10])
    // Signal:     s == 0        (Sigma)
    // Background: s == 1,2,3     (other hyperons + Background MC + Dirt)
    // ======================================================
    if (Cuts[10] && (s == 0 || s == 1 || s == 2 || s == 3)) {

      // Reset to sentinels
      t_MuonPID=t_MuonTrkLength=t_MuonTrkDistVtx=t_MuonTrkScore=-999;
      t_MuonDeDxP0=t_MuonDeDxP1=t_MuonDeDxP2=-999;
      t_MuonNHitsY=t_MuonPhotonPID=-999;

      t_PionPID=t_PionTrkLength=t_PionTrkDistVtx=t_PionTrkScore=-999;
      t_PionDeDxP0=t_PionDeDxP1=t_PionDeDxP2=-999;
      t_PionNHitsY=t_PionPhotonPID=-999;

      t_ProtonPID=t_ProtonTrkLength=t_ProtonTrkDistVtx=t_ProtonTrkScore=-999;
      t_ProtonDeDxP0=t_ProtonDeDxP1=t_ProtonDeDxP2=-999;
      t_ProtonNHitsY=t_ProtonPhotonPID=-999;

      t_Shr1NHitsY=t_Shr1VtxDist=t_Shr1DeDxP0=t_Shr1DeDxP1=t_Shr1DeDxP2=t_Shr1OpenAngle=-999;
      t_Shr2NHitsY=t_Shr2VtxDist=t_Shr2DeDxP0=t_Shr2DeDxP1=t_Shr2DeDxP2=t_Shr2OpenAngle=-999;

      t_AngleMuonPion=t_AngleProtonPion=t_AngleProtonMuon=-999;
      t_AngleShr1Muon=t_AngleShr1Pion=t_AngleShr1Proton=-999;
      t_AngleShr2Muon=t_AngleShr2Pion=t_AngleShr2Proton=-999;

      t_DistMuonPion=t_DistProtonPion=t_DistProtonMuon=-999;
      t_DistShr1Muon=t_DistShr1Pion=t_DistShr1Proton=-999;
      t_DistShr2Muon=t_DistShr2Pion=t_DistShr2Proton=-999;

      t_NTracks  = (Float_t)ntracks;
      t_NShowers = (Float_t)nshowers;
      t_Weight   = (Float_t)scale[s];

      // Muon
      if (muon_index > -1) {
        TVector3 st(trk_start_x->at(muon_index), trk_start_y->at(muon_index), trk_start_z->at(muon_index));
        t_MuonPID        = SafeGet(trk_llrpid,          muon_index);
        t_MuonTrkLength  = SafeGet(trk_length,           muon_index);
        t_MuonTrkDistVtx = (float)(st - vtx).Mag();
        t_MuonTrkScore   = SafeGet(pfp_trk_shr_score,   muon_index);
        t_MuonDeDxP0     = SafeGet(trk_mean_dedx_plane0,muon_index);
        t_MuonDeDxP1     = SafeGet(trk_mean_dedx_plane1,muon_index);
        t_MuonDeDxP2     = SafeGet(trk_mean_dedx_plane2,muon_index);
        t_MuonNHitsY     = (float)SafeGetI(pfpnplanehits_Y, muon_index);
        t_MuonPhotonPID  = SafeGet(pfp_photon_bdt_score, muon_index);
      }
      // Pion
      if (pion_index > -1) {
        TVector3 st(trk_start_x->at(pion_index), trk_start_y->at(pion_index), trk_start_z->at(pion_index));
        t_PionPID        = SafeGet(trk_llrpid,           pion_index);
        t_PionTrkLength  = SafeGet(trk_length,            pion_index);
        t_PionTrkDistVtx = (float)(st - vtx).Mag();
        t_PionTrkScore   = SafeGet(pfp_trk_shr_score,    pion_index);
        t_PionDeDxP0     = SafeGet(trk_mean_dedx_plane0, pion_index);
        t_PionDeDxP1     = SafeGet(trk_mean_dedx_plane1, pion_index);
        t_PionDeDxP2     = SafeGet(trk_mean_dedx_plane2, pion_index);
        t_PionNHitsY     = (float)SafeGetI(pfpnplanehits_Y, pion_index);
        t_PionPhotonPID  = SafeGet(pfp_photon_bdt_score,  pion_index);
      }
      // Proton
      if (proton_index > -1) {
        TVector3 st(trk_start_x->at(proton_index), trk_start_y->at(proton_index), trk_start_z->at(proton_index));
        t_ProtonPID        = SafeGet(trk_llrpid,           proton_index);
        t_ProtonTrkLength  = SafeGet(trk_length,            proton_index);
        t_ProtonTrkDistVtx = (float)(st - vtx).Mag();
        t_ProtonTrkScore   = SafeGet(pfp_trk_shr_score,    proton_index);
        t_ProtonDeDxP0     = SafeGet(trk_mean_dedx_plane0, proton_index);
        t_ProtonDeDxP1     = SafeGet(trk_mean_dedx_plane1, proton_index);
        t_ProtonDeDxP2     = SafeGet(trk_mean_dedx_plane2, proton_index);
        t_ProtonNHitsY     = (float)SafeGetI(pfpnplanehits_Y, proton_index);
        t_ProtonPhotonPID  = SafeGet(pfp_photon_bdt_score,  proton_index);
      }
      // Shower 1
      if (shower1_index > -1) {
        t_Shr1NHitsY    = (float)SafeGetI(pfpnplanehits_Y, shower1_index);
        t_Shr1VtxDist   = (float)(shower1Start - vtx).Mag();
        t_Shr1DeDxP0    = SafeGet(shr_dedx_plane0,  shower1_index);
        t_Shr1DeDxP1    = SafeGet(shr_dedx_plane1,  shower1_index);
        t_Shr1DeDxP2    = SafeGet(shr_dedx_plane2,  shower1_index);
        t_Shr1OpenAngle = SafeGet(shr_open_angle,    shower1_index);
      }
      // Shower 2
      if (shower2_index > -1) {
        t_Shr2NHitsY    = (float)SafeGetI(pfpnplanehits_Y, shower2_index);
        t_Shr2VtxDist   = (float)(shower2Start - vtx).Mag();
        t_Shr2DeDxP0    = SafeGet(shr_dedx_plane0,  shower2_index);
        t_Shr2DeDxP1    = SafeGet(shr_dedx_plane1,  shower2_index);
        t_Shr2DeDxP2    = SafeGet(shr_dedx_plane2,  shower2_index);
        t_Shr2OpenAngle = SafeGet(shr_open_angle,    shower2_index);
      }

      // Opening angles
      if (muon_index>-1  && pion_index>-1)    t_AngleMuonPion   = (float)muonDir.Angle(pionDir);
      if (proton_index>-1&& pion_index>-1)    t_AngleProtonPion = (float)protonDir.Angle(pionDir);
      if (proton_index>-1&& muon_index>-1)    t_AngleProtonMuon = (float)protonDir.Angle(muonDir);
      if (shower1_index>-1&&muon_index>-1)    t_AngleShr1Muon   = (float)shower1Dir.Angle(muonDir);
      if (shower1_index>-1&&pion_index>-1)    t_AngleShr1Pion   = (float)shower1Dir.Angle(pionDir);
      if (shower1_index>-1&&proton_index>-1)  t_AngleShr1Proton = (float)shower1Dir.Angle(protonDir);
      if (shower2_index>-1&&muon_index>-1)    t_AngleShr2Muon   = (float)shower2Dir.Angle(muonDir);
      if (shower2_index>-1&&pion_index>-1)    t_AngleShr2Pion   = (float)shower2Dir.Angle(pionDir);
      if (shower2_index>-1&&proton_index>-1)  t_AngleShr2Proton = (float)shower2Dir.Angle(protonDir);

      // Distances
      t_DistMuonPion   = (float)distMuonPion;
      t_DistProtonPion = (float)distProtonPion;
      t_DistProtonMuon = (float)distProtonMuon;
      t_DistShr1Muon   = (float)distShower1Muon;
      t_DistShr1Pion   = (float)distShower1Pion;
      t_DistShr1Proton = (float)distShower1Prot;
      t_DistShr2Muon   = (float)distShower2Muon;
      t_DistShr2Pion   = (float)distShower2Pion;
      t_DistShr2Proton = (float)distShower2Prot;

      // Fill the appropriate TMVA tree
      if (s == 0)            TreeS->Fill();  // Signal
      else /* s==1,2,3 */    TreeB->Fill();  // Background (hyperons + bkg + dirt)
    }

  } // end main event loop

  // -----------------------------------------------------------
  // Cut-flow summary table  (FIX: includes Cuts[10])
  // -----------------------------------------------------------
  cout << "\ncut \t\t\t Sig \t Hyp \t Bkg \t Dirt \t Eff[%] \t Purity[%]" << endl;
  for (int c = 0; c < Nc; c++) {
    double sig  = Hist[0][0][c]->Integral();
    double back = Hist[0][1][c]->Integral() + Hist[0][2][c]->Integral() + Hist[0][3][c]->Integral();
    double sig0 = Hist[0][0][0]->Integral();
    cout << CutName[c]
         << "\t\t"  << sig
         << "\t"    << Hist[0][1][c]->Integral()
         << "\t"    << Hist[0][2][c]->Integral()
         << "\t"    << Hist[0][3][c]->Integral()
         << "\t"    << (sig0 > 0 ? sig / sig0 * 100 : 0)
         << "\t"    << ((sig + back) > 0 ? sig / (sig + back) * 100 : 0)
         << endl;
  }

  cout << "\nTMVA TreeS entries: " << TreeS->GetEntries()
       << "   TreeB entries: "     << TreeB->GetEntries() << endl;

  // -----------------------------------------------------------
  // Write TMVA trees and close the output file
  // -----------------------------------------------------------
  tmvaFile->cd();
  TreeS->Write();
  TreeB->Write();
  tmvaFile->Close();
  cout << "TMVA input written to TMVA_input.root" << endl;

  // -----------------------------------------------------------
  // Write all histograms to SigmaHists.root (one TDirectory per cut).
  // Plotting is now done separately by PlotSigma.C, which reads this file.
  // -----------------------------------------------------------
  TFile *histFile = new TFile("SigmaHists.root", "RECREATE");
  for (int c = 0; c < Nc; c++) {
    TDirectory *d = histFile->mkdir(CutName[c]);
    d->cd();
    for (int v = 0; v < Nv; v++)
      for (int s = 0; s < Ns; s++)
        Hist[v][s][c]->Write();
  }
  histFile->Close();
  cout << "Histograms written to SigmaHists.root" << endl;
}
