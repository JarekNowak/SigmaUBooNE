#define Sigma_cxx
#include "Sigma.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

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
  const int Ns = 5;
  TString Names[Ns] = {"Sigma", "Hyperons", "Background", "Dirt", "Data"};

  // Cumulative event boundaries (used to assign sample index via jentry)
  int NEvents[Ns] = {0};
  NEvents[0] = 262576;                       // Sigma
  NEvents[1] = NEvents[0];                   // Hyperons share the same file range
  NEvents[2] = 651830  + NEvents[0];         // Background (cumulative)
  NEvents[3] = 349879  + NEvents[2];         // Dirt       (cumulative)
  NEvents[4] = 0       + NEvents[3];         // Data       (cumulative)

  for (int i = 0; i < Ns; i++) cout << Names[i] << "  NEvents=" << NEvents[i] << endl;

  // Protons on Target — used for luminosity scaling to run-3b data
  double PoT[Ns] = {0};
  PoT[0] = 552.0;   // *1e20  Sigma MC
  PoT[1] = PoT[0];  //        Hyperons MC  (same sample)
  PoT[2] = 9.280;   // *1e20  Background MC
  PoT[3] = 17.30;   // *1e20  Dirt MC
  PoT[4] = 4.90;    // *1e20  Data (run 3b)

  double scale[Ns] = {0};
  for (int i = 0; i < Ns - 1; i++) scale[i] = PoT[4] / PoT[i];
  scale[4] = 1.0;

  for (int i = 0; i < Ns; i++) cout << Names[i] << "  scale=" << scale[i] << endl;

  // -----------------------------------------------------------
  // Cut names  (Nc = 11)
  // -----------------------------------------------------------
  const int Nc = 11;
  TString CutName[Nc] = {
    "AllEvents",
    "FiducialV",          // FIX: was "FducialV"
    "NShowers",
    "NoBigShowers",
    "MuonCandidate",
    "PionCandidate",
    "ProtonCandidate",
    "NTracksCut",
    "2ndShower",
    "1stShower",
    "Distances"
  };

  // -----------------------------------------------------------
  // Variable names  (Nv = 80)
  // -----------------------------------------------------------
  const int Nv = 80;
  TString VarName[Nv];
  VarName[0]  = "NTracks";
  VarName[1]  = "NShowers";
  VarName[2]  = "NHitsY";
  VarName[3]  = "MuonPID";
  VarName[4]  = "MuonTrkLength";       // FIX: was "MuonTrkLenght"
  VarName[5]  = "MuonTrkDistVtx";
  VarName[6]  = "MuonTrkScore";
  VarName[7]  = "PionPID";
  VarName[8]  = "PionTrkLength";       // FIX: was "PionTrkLenght"
  VarName[9]  = "PionTrkDistVtx";
  VarName[10] = "PionTrkScore";
  VarName[11] = "ProtonPID";
  VarName[12] = "ProtonTrkLength";     // FIX: was "ProtonTrkLenght"
  VarName[13] = "ProtonTrkDistVtx";
  VarName[14] = "ProtonTrkScore";
  VarName[15] = "GammaPID";
  VarName[16] = "GammaTrkLength";      // FIX: was "GammaTrkLenght"
  VarName[17] = "GammaTrkDistVtx";
  VarName[18] = "GammaTrkScore";
  VarName[19] = "MuonNHitsY";
  VarName[20] = "PionNHitsY";
  VarName[21] = "ProtonNHitsY";
  VarName[22] = "GammaNHitsY";
  VarName[23] = "MuonNHitsU";
  VarName[24] = "PionNHitsU";
  VarName[25] = "ProtonNHitsU";
  VarName[26] = "GammaNHitsU";
  VarName[27] = "MuonNHitsV";
  VarName[28] = "PionNHitsV";
  VarName[29] = "ProtonNHitsV";
  VarName[30] = "GammaNHitsV";
  VarName[31] = "MuonMeanDeDxPlane0";
  VarName[32] = "MuonMeanDeDxPlane1";
  VarName[33] = "MuonMeanDeDxPlane2";
  VarName[34] = "PionMeanDeDxPlane0";
  VarName[35] = "PionMeanDeDxPlane1";
  VarName[36] = "PionMeanDeDxPlane2";
  VarName[37] = "ProtonMeanDeDxPlane0";
  VarName[38] = "ProtonMeanDeDxPlane1";
  VarName[39] = "ProtonMeanDeDxPlane2";
  VarName[40] = "GammaMeanDeDxPlane0";
  VarName[41] = "GammaMeanDeDxPlane1";
  VarName[42] = "GammaMeanDeDxPlane2";
  VarName[43] = "MuonPhotonPID";
  VarName[44] = "PionPhotonPID";
  VarName[45] = "ProtonPhotonPID";
  VarName[46] = "GammaPhotonPID";
  VarName[47] = "MuonCandShrOpenAngle";
  VarName[48] = "PionCandShrOpenAngle";
  VarName[49] = "ProtonCandShrOpenAngle";
  VarName[50] = "Shower1CandShrOpenAngle";
  VarName[51] = "Shower2CandShrOpenAngle";
  VarName[52] = "Shower1CandDeDxPlane0";  // FIX: was "Plan"
  VarName[53] = "Shower1CandDeDxPlane1";
  VarName[54] = "Shower1CandDeDxPlane2";
  VarName[55] = "Shower2CandDeDxPlane0";
  VarName[56] = "Shower2CandDeDxPlane1";  // FIX: was "Shower3"
  VarName[57] = "Shower2CandDeDxPlane2";  // FIX: was "Shower4"
  VarName[58] = "OpeningAngleMuonPion";
  VarName[59] = "OpeningAngleProtonPion";
  VarName[60] = "OpeningAngleProtonMuon";
  VarName[61] = "OpeningAngleShower1Muon";
  VarName[62] = "OpeningAngleShower1Pion";
  VarName[63] = "OpeningAngleShower1Proton";
  VarName[64] = "OpeningAngleShower2Muon";
  VarName[65] = "OpeningAngleShower2Pion";
  VarName[66] = "OpeningAngleShower2Proton";
  VarName[67] = "Shower1NHitsY";
  VarName[68] = "Shower2NHitsY";
  VarName[69] = "Shower1VtxDist";
  VarName[70] = "Shower2VtxDist";
  VarName[71] = "DistMuonPion";
  VarName[72] = "DistProtonPion";
  VarName[73] = "DistProtonMuon";
  VarName[74] = "DistShower1Muon";
  VarName[75] = "DistShower1Pion";
  VarName[76] = "DistShower1Proton";
  VarName[77] = "DistShower2Muon";
  VarName[78] = "DistShower2Pion";
  VarName[79] = "DistShower2Proton";

  // -----------------------------------------------------------
  // Histogram booking — use (0,0) range for ROOT auto-range
  // -----------------------------------------------------------
  TH1F *Hist[Nv][Ns][Nc];
  for (int s = 0; s < Ns; s++)
    for (int v = 0; v < Nv; v++)
      for (int c = 0; c < Nc; c++)
        Hist[v][s][c] = new TH1F(VarName[v] + Names[s] + CutName[c], "", 100, 0, 0);

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
  // Guard against empty chain
  // -----------------------------------------------------------
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // -----------------------------------------------------------
  // Main event loop
  // -----------------------------------------------------------
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

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

    // ---- Truth-level decay product counting ----
    for (int i = 0; i < (int)mc_decay_pdg->size(); i++) {
      if (mc_decay_pdg->at(i) ==  2212) mc_decay++;  // proton
      if (mc_decay_pdg->at(i) == -211)  mc_decay++;  // pi-
    }

    // ---- FIX: corrected fiducial volume cut (all three axes, no impossible condition) ----
    // Skip events whose true neutrino vertex lies OUTSIDE the FV
    bool inFV_truth = (mc_nu_pos_x >   3.00  && mc_nu_pos_x < 253.35 &&
                       mc_nu_pos_y > -112.53  && mc_nu_pos_y < 114.47 &&
                       mc_nu_pos_z >   3.1    && mc_nu_pos_z < 1033.0);
    bool inDeadRegion = (mc_nu_pos_z > 675.1  && mc_nu_pos_z < 775.1);
    if (!inFV_truth || inDeadRegion) continue;

    // ---- Signal definition: anti-nu_mu + Sigma0 with proton + pi- decay products ----
    if ((mc_nu_pdg == -14) && (mc_hyperon_pdg == 3212) && (mc_decay == 2)) IsSignal = true;

    // ---- Sample assignment ----
    if (IsSignal  && jentry <= NEvents[0]) s = 0;  // Sigma signal
    if (!IsSignal && jentry <= NEvents[0]) s = 1;  // Hyperons (same file, not signal)
    if (jentry > NEvents[0] && jentry <= NEvents[2]) s = 2;  // Background
    if (jentry > NEvents[2] && jentry <= NEvents[3]) s = 3;  // Dirt
    // Data (s=4) not yet available per original code

    if (s < 0) continue;  // Safety: skip unclassified events

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

    // ---- Shower candidates: closest to vertex = shower1, second = shower2 ----
    double shr_dist1 = 1e9;
    for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
      if (i == muon_index || i == pion_index || i == proton_index) continue;
      if (pfp_trk_shr_score->at(i) < 0.5) {
        TVector3 start(trk_start_x->at(i), trk_start_y->at(i), trk_start_z->at(i));
        double dist = (start - vtx).Mag();
        if (shower1_index == -1) {
          shower1_index = i;
          shr_dist1     = dist;
        } else if (dist < shr_dist1) {
          shower2_index = shower1_index;
          shower1_index = i;
          shr_dist1     = dist;
        }
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

    if (nshowers > 0 && nshowers < 3) sel_nshowers = true;
    if (ntracks < 5)                  sel_ntracks   = true;

    // ---- FIX: corrected bitwise & -> && in distance cut ----
    bool sel_distances =
      (distMuonPion < 30 && distProtonMuon < 20 && distProtonPion < 20 &&
       distShower1Muon < 100 && distShower1Pion < 100 && distShower1Prot < 100);

    // ---- Build cumulative cut flags ----
    Cuts[0]  = true;
    Cuts[1]  = Cuts[0]  && sel_FV;
    Cuts[2]  = Cuts[1]  && sel_nshowers;
    Cuts[3]  = Cuts[2]  && NoBigShowers;
    Cuts[4]  = Cuts[3]  && sel_muon_candidate;
    Cuts[5]  = Cuts[4]  && sel_pion_candidate;
    Cuts[6]  = Cuts[5]  && sel_proton_candidate;
    Cuts[7]  = Cuts[6]  && sel_ntracks;
    Cuts[8]  = Cuts[7]  && sel_more_showers;
    Cuts[9]  = Cuts[8]  && sel_1st_shower;
    Cuts[10] = Cuts[9] && sel_distances;

    // ======================================================
    // Fill histograms (loop over cuts, then over PFPs)
    // ======================================================
    for (int c = 0; c < Nc; c++) {
      if (!Cuts[c]) continue;

      // Topological
      Hist[0][s][c]->Fill(ntracks,  scale[s]);
      Hist[1][s][c]->Fill(nshowers, scale[s]);

      // Shower hit counts in Y plane
      if (nshowers > 0)
        for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++)
          if (pfp_trk_shr_score->at(i) < 0.5)
            Hist[2][s][c]->Fill(pfpnplanehits_Y->at(i), scale[s]);

      // Opening angles and distances (pre-computed above)
      if (muon_index>-1  && pion_index>-1)    { Hist[58][s][c]->Fill(muonDir.Angle(pionDir),     scale[s]); Hist[71][s][c]->Fill(distMuonPion,    scale[s]); }
      if (proton_index>-1&& pion_index>-1)    { Hist[59][s][c]->Fill(protonDir.Angle(pionDir),   scale[s]); Hist[72][s][c]->Fill(distProtonPion,  scale[s]); }
      if (proton_index>-1&& muon_index>-1)    { Hist[60][s][c]->Fill(protonDir.Angle(muonDir),   scale[s]); Hist[73][s][c]->Fill(distProtonMuon,  scale[s]); }
      if (shower1_index>-1&&muon_index>-1)    { Hist[61][s][c]->Fill(shower1Dir.Angle(muonDir),  scale[s]); Hist[74][s][c]->Fill(distShower1Muon, scale[s]); }
      if (shower1_index>-1&&pion_index>-1)    { Hist[62][s][c]->Fill(shower1Dir.Angle(pionDir),  scale[s]); Hist[75][s][c]->Fill(distShower1Pion, scale[s]); }
      if (shower1_index>-1&&proton_index>-1)  { Hist[63][s][c]->Fill(shower1Dir.Angle(protonDir),scale[s]); Hist[76][s][c]->Fill(distShower1Prot, scale[s]); }
      if (shower2_index>-1&&muon_index>-1)    { Hist[64][s][c]->Fill(shower2Dir.Angle(muonDir),  scale[s]); Hist[77][s][c]->Fill(distShower2Muon, scale[s]); }
      if (shower2_index>-1&&pion_index>-1)    { Hist[65][s][c]->Fill(shower2Dir.Angle(pionDir),  scale[s]); Hist[78][s][c]->Fill(distShower2Pion, scale[s]); }
      if (shower2_index>-1&&proton_index>-1)  { Hist[66][s][c]->Fill(shower2Dir.Angle(protonDir),scale[s]); Hist[79][s][c]->Fill(distShower2Prot, scale[s]); }

      // Shower open angles and dEdx
      if (muon_index    > -1) Hist[47][s][c]->Fill(SafeGet(shr_open_angle, muon_index),    scale[s]);
      if (pion_index    > -1) Hist[48][s][c]->Fill(SafeGet(shr_open_angle, pion_index),    scale[s]);
      if (proton_index  > -1) Hist[49][s][c]->Fill(SafeGet(shr_open_angle, proton_index),  scale[s]);
      if (shower1_index > -1) {
        Hist[50][s][c]->Fill(SafeGet(shr_open_angle,  shower1_index), scale[s]);
        Hist[52][s][c]->Fill(SafeGet(shr_dedx_plane0, shower1_index), scale[s]);
        Hist[53][s][c]->Fill(SafeGet(shr_dedx_plane1, shower1_index), scale[s]);
        Hist[54][s][c]->Fill(SafeGet(shr_dedx_plane2, shower1_index), scale[s]);
        Hist[67][s][c]->Fill(SafeGetI(pfpnplanehits_Y, shower1_index),scale[s]);
        Hist[69][s][c]->Fill((shower1Start - vtx).Mag(),               scale[s]);
      }
      if (shower2_index > -1) {
        Hist[51][s][c]->Fill(SafeGet(shr_open_angle,  shower2_index), scale[s]);
        Hist[55][s][c]->Fill(SafeGet(shr_dedx_plane0, shower2_index), scale[s]);
        Hist[56][s][c]->Fill(SafeGet(shr_dedx_plane1, shower2_index), scale[s]);
        Hist[57][s][c]->Fill(SafeGet(shr_dedx_plane2, shower2_index), scale[s]);
        Hist[68][s][c]->Fill(SafeGetI(pfpnplanehits_Y, shower2_index),scale[s]);
        Hist[70][s][c]->Fill((shower2Start - vtx).Mag(),               scale[s]);
      }

      // Per-PFP particle type histograms
      for (int i = 0; i < (int)pfp_trk_shr_score->size(); i++) {
        if (trk_length->at(i) <= 0) continue;
        TVector3 start(trk_start_x->at(i), trk_start_y->at(i), trk_start_z->at(i));
        double dist = (start - vtx).Mag();
        int pdg = pfp_true_pdg->at(i);

        int base = -1;
        if      (pdg == -13)  base = 0;   // muon     -> indices 3-6, 19-21, 31-33, 43
        else if (pdg == -211) base = 1;   // pion     -> indices 7-10, 20-22, 34-36, 44
        else if (pdg ==  2212)base = 2;   // proton   -> indices 11-14, 21-23, 37-39, 45
        else if (pdg ==  22)  base = 3;   // gamma    -> indices 15-18, 22-24, 40-42, 46
        if (base < 0) continue;

        int off = base * 4;  // stride for PID/Length/Dist/Score block
        Hist[3  + off][s][c]->Fill(trk_llrpid->at(i),             scale[s]);
        Hist[4  + off][s][c]->Fill(trk_length->at(i),             scale[s]);
        Hist[5  + off][s][c]->Fill(dist,                           scale[s]);
        Hist[6  + off][s][c]->Fill(pfp_trk_shr_score->at(i),      scale[s]);
        Hist[19 + base][s][c]->Fill(pfpnplanehits_Y->at(i),        scale[s]);
        Hist[23 + base][s][c]->Fill(pfpnplanehits_U->at(i),        scale[s]);
        Hist[27 + base][s][c]->Fill(pfpnplanehits_V->at(i),        scale[s]);
        Hist[31 + base*3][s][c]->Fill(trk_mean_dedx_plane0->at(i), scale[s]);
        Hist[32 + base*3][s][c]->Fill(trk_mean_dedx_plane1->at(i), scale[s]);
        Hist[33 + base*3][s][c]->Fill(trk_mean_dedx_plane2->at(i), scale[s]);
        Hist[43 + base][s][c]->Fill(pfp_photon_bdt_score->at(i),   scale[s]);
      }
    } // end cut loop

    // ======================================================
    // TMVA tree filling — events passing ALL cuts (Cuts[10])
    // Signal:     s == 0  (Sigma)
    // Background: s == 2 or s == 3  (Background MC + Dirt)
    // ======================================================
    if (Cuts[10] && (s == 0 || s == 2 || s == 3)) {

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
      if (s == 0)          TreeS->Fill();  // Signal
      else /* s==2||s==3*/ TreeB->Fill();  // Background
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
         << "\t"    << (back > 0 ? sig / back * 100 : 0)
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
  // Plot all variables across all samples and cuts
  // -----------------------------------------------------------
 /* TCanvas *c2 = new TCanvas("c2", "", 5000, 4000);
  for (int v = 0; v < Nv; v++) {
    c2->Divide(Nc, 4);
    for (int c = 0; c < Nc; c++)
      for (int s = 0; s < 4; s++) {
        c2->cd(c + 1 + s * Nc);
        Hist[v][s][c]->Draw();
      }
    c2->Print("figs/" + VarName[v] + ".png");
    c2->Clear();
  }
    */
}
