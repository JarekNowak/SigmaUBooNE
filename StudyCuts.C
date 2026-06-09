// ============================================================
//  StudyCuts.C
//  Single-pass study of pre-MVA cut optimisation.
//
//  For each scenario (skip a cut / loosen a threshold) the macro
//  evaluates what the final BDTg > 0.96 yield would be, so we
//  can identify redundant cuts and efficiency-recovery opportunities.
//
//  Usage:
//    root -b -l -q StudyCuts.C
// ============================================================

#define Sigma_cxx
#include "Sigma.h"
void Sigma::Loop() {}   // stub — Loop() is defined in Sigma.C; not used here
#include "TMVA/Reader.h"
#include "SigmaConfig.h"

#include <TVector3.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

// ---- SafeGet helpers (same as Sigma.C) ----
template<typename T>
static inline float SafeGet(const vector<T>* v, int i, float s=-999.f)
{ return (!v||i<0||i>=(int)v->size()) ? s : (float)v->at(i); }

template<typename T>
static inline int SafeGetI(const vector<T>* v, int i, int s=-999)
{ return (!v||i<0||i>=(int)v->size()) ? s : (int)v->at(i); }

// ---- Per-scenario yield container ----
struct Scenario {
    string name;
    double S, Hyp, Bkg, Dirt;   // weighted counts after final MVA cut
    Scenario(const string& n) : name(n), S(0), Hyp(0), Bkg(0), Dirt(0) {}
    double B()   const { return Hyp+Bkg+Dirt; }
    double eff(double S0) const { return S0>0 ? S/S0*100 : 0; }
    double pur() const { double T=S+B(); return T>0 ? S/T*100 : 0; }
    double fom() const { return B()>0 ? S/sqrt(B()) : 0; }
};

// ---- Scenario indices ----
enum {
    kBase=0,
    // Cut-removal studies
    kNoNShowers, kNoNoBigShr, kNoNTracks, kNo2ndShr, kNo1stShr, kNoDist,
    kNoMuon, kNoPion, kNoProton,
    // Threshold loosening
    kMuPID05, kMuPID06,          // muon LLR PID threshold
    kMuScore05, kMuScore06,      // muon track score threshold
    kMuDist30,                   // muon distance to vtx
    kShr3,                       // allow up to 3 showers
    kNTracks5, kNTracks6,        // ntracks upper bound
    kShr1HitsMax400,             // widen 1st shower hit window
    kShr1HitsMin5,
    kDistMuPi40, kDistMuPi50,   // muon-pion distance cut
    kDistProt30,                 // proton-X distance cut
    kDistShr150,                 // shower-candidate distance cut
    kNScenarios
};

void StudyCuts()
{
    // -----------------------------------------------------------
    // PoT scales (same as Sigma.C)
    // -----------------------------------------------------------
    const int Ns = 5;
    double PoT[Ns] = { 552.0, 552.0, 9.280, 17.30, 4.90 };
    double scale[Ns];
    for (int i = 0; i < Ns-1; i++) scale[i] = PoT[4] / PoT[i];
    scale[4] = 1.0;

    // -----------------------------------------------------------
    // TMVA Reader — identical variable order to RunTMVA.C
    // -----------------------------------------------------------
    Float_t r_MuonTrkLength, r_MuonTrkDistVtx, r_MuonTrkScore;
    Float_t r_MuonNHitsY, r_MuonPhotonPID;
    Float_t r_PionPID, r_PionTrkLength, r_PionTrkDistVtx, r_PionTrkScore;
    Float_t r_PionNHitsY, r_PionDeDxP2, r_PionPhotonPID;
    Float_t r_ProtonPID, r_ProtonTrkLength, r_ProtonTrkDistVtx, r_ProtonTrkScore;
    Float_t r_ProtonNHitsY, r_ProtonDeDxP2, r_ProtonPhotonPID;
    Float_t r_Shr1NHitsY, r_Shr1VtxDist, r_Shr1OpenAngle;
    Float_t r_Shr2NHitsY, r_Shr2VtxDist;
    Float_t r_AngleMuonPion, r_AngleProtonPion, r_AngleProtonMuon;
    Float_t r_AngleShr1Muon, r_AngleShr1Pion, r_AngleShr1Proton;
    Float_t r_DistMuonPion, r_DistProtonPion, r_DistProtonMuon;
    Float_t r_Weight;

    TMVA::Reader *rd = new TMVA::Reader("!Color:!Silent");
    rd->AddVariable("MuonTrkLength",    &r_MuonTrkLength);
    rd->AddVariable("MuonTrkDistVtx",   &r_MuonTrkDistVtx);
    rd->AddVariable("MuonTrkScore",     &r_MuonTrkScore);
    rd->AddVariable("MuonNHitsY",       &r_MuonNHitsY);
    rd->AddVariable("MuonPhotonPID",    &r_MuonPhotonPID);
    rd->AddVariable("PionPID",          &r_PionPID);
    rd->AddVariable("PionTrkLength",    &r_PionTrkLength);
    rd->AddVariable("PionTrkDistVtx",   &r_PionTrkDistVtx);
    rd->AddVariable("PionTrkScore",     &r_PionTrkScore);
    rd->AddVariable("PionNHitsY",       &r_PionNHitsY);
    rd->AddVariable("PionDeDxPlane2",   &r_PionDeDxP2);
    rd->AddVariable("PionPhotonPID",    &r_PionPhotonPID);
    rd->AddVariable("ProtonPID",        &r_ProtonPID);
    rd->AddVariable("ProtonTrkLength",  &r_ProtonTrkLength);
    rd->AddVariable("ProtonTrkDistVtx", &r_ProtonTrkDistVtx);
    rd->AddVariable("ProtonTrkScore",   &r_ProtonTrkScore);
    rd->AddVariable("ProtonNHitsY",     &r_ProtonNHitsY);
    rd->AddVariable("ProtonDeDxPlane2", &r_ProtonDeDxP2);
    rd->AddVariable("ProtonPhotonPID",  &r_ProtonPhotonPID);
    rd->AddVariable("Shr1NHitsY",       &r_Shr1NHitsY);
    rd->AddVariable("Shr1VtxDist",      &r_Shr1VtxDist);
    rd->AddVariable("Shr1OpenAngle",    &r_Shr1OpenAngle);
    rd->AddVariable("Shr2NHitsY",       &r_Shr2NHitsY);
    rd->AddVariable("Shr2VtxDist",      &r_Shr2VtxDist);
    rd->AddVariable("AngleMuonPion",    &r_AngleMuonPion);
    rd->AddVariable("AngleProtonPion",  &r_AngleProtonPion);
    rd->AddVariable("AngleProtonMuon",  &r_AngleProtonMuon);
    rd->AddVariable("AngleShr1Muon",    &r_AngleShr1Muon);
    rd->AddVariable("AngleShr1Pion",    &r_AngleShr1Pion);
    rd->AddVariable("AngleShr1Proton",  &r_AngleShr1Proton);
    rd->AddVariable("DistMuonPion",     &r_DistMuonPion);
    rd->AddVariable("DistProtonPion",   &r_DistProtonPion);
    rd->AddVariable("DistProtonMuon",   &r_DistProtonMuon);
    rd->AddSpectator("Weight",          &r_Weight);

    rd->BookMVA("BDT",    "SigmaDataset/weights/TMVASigma_BDT.weights.xml");
    rd->BookMVA("BDTAda", "SigmaDataset/weights/TMVASigma_BDTAda.weights.xml");
    rd->BookMVA("BDTg",   "SigmaDataset/weights/TMVASigma_BDTg.weights.xml");

    // -----------------------------------------------------------
    // Scenario containers
    // -----------------------------------------------------------
    vector<Scenario> sc;
    sc.emplace_back("Baseline (current)");
    sc.emplace_back("Skip NShowers");
    sc.emplace_back("Skip NoBigShowers");
    sc.emplace_back("Skip NTracksCut");
    sc.emplace_back("Skip 2ndShower veto");
    sc.emplace_back("Skip 1stShower");
    sc.emplace_back("Skip Distances");
    sc.emplace_back("Skip MuonCandidate");
    sc.emplace_back("Skip PionCandidate");
    sc.emplace_back("Skip ProtonCandidate");
    sc.emplace_back("MuonPID > 0.5");
    sc.emplace_back("MuonPID > 0.6");
    sc.emplace_back("MuonTrkScore > 0.5");
    sc.emplace_back("MuonTrkScore > 0.6");
    sc.emplace_back("MuonDistVtx < 30");
    sc.emplace_back("NShowers <= 3");
    sc.emplace_back("NTracks < 5");
    sc.emplace_back("NTracks < 6");
    sc.emplace_back("Shr1 NHitsY < 400");
    sc.emplace_back("Shr1 NHitsY > 5");
    sc.emplace_back("DistMuonPion < 40");
    sc.emplace_back("DistMuonPion < 50");
    sc.emplace_back("DistProton-X < 30");
    sc.emplace_back("DistShr1-X < 150");

    // -----------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------
    Sigma a;
    if (!a.fChain) { cerr << "No chain\n"; return; }
    Long64_t nentries = a.fChain->GetEntries();

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        if (a.LoadTree(jentry) < 0) break;
        a.fChain->GetEntry(jentry);

        if (jentry % 50000 == 0)
            cout << "\r" << 100*(jentry+1)/nentries << "%" << flush;

        // ---- Sample assignment ----
        int mc_decay = 0;
        double proton_mom = 0, pion_mom = 0;
        for (int i = 0; i < (int)a.mc_decay_pdg->size(); i++) {
            if (a.mc_decay_pdg->at(i) ==  2212) { mc_decay++; proton_mom = SafeGet(a.mc_decay_mom,i); }
            if (a.mc_decay_pdg->at(i) == -211)  { mc_decay++; pion_mom   = SafeGet(a.mc_decay_mom,i); }
        }
        bool IsSignal = (a.mc_nu_pdg==-14 && a.mc_hyperon_pdg==3212 && mc_decay==2
                         && proton_mom>0.3 && pion_mom>0.1);

        int treenum = a.fChain->GetTreeNumber();
        int s = -1;
        if      (treenum==0) s = IsSignal ? 0 : 1;
        else if (treenum==1) s = 2;
        else if (treenum==2) s = 3;
        if (s < 0) continue;

        bool inFV_truth = (a.mc_nu_pos_x>3.0 && a.mc_nu_pos_x<253.35 &&
                           a.mc_nu_pos_y>-112.53 && a.mc_nu_pos_y<114.47 &&
                           a.mc_nu_pos_z>3.1 && a.mc_nu_pos_z<1033.0);
        bool inDead = (a.mc_nu_pos_z>675.1 && a.mc_nu_pos_z<775.1);
        if ((s==0||s==1) && (!inFV_truth||inDead)) continue;

        double w = scale[s];

        // ---- Reco quantities (same as Sigma.C) ----
        TVector3 vtx(a.reco_primary_vtx_x, a.reco_primary_vtx_y, a.reco_primary_vtx_z);
        bool sel_FV = a.reco_primary_vtx_inFV;
        if (!sel_FV) continue;   // FV is kept in all scenarios (safety cut)

        int ntracks=0, nshowers=0;
        bool NoBigShowers = true;
        for (int i=0; i<(int)a.pfp_trk_shr_score->size(); i++) {
            if (a.pfp_trk_shr_score->at(i)<0.5) {
                nshowers++;
                if (a.pfpnplanehits_Y->at(i)>600) NoBigShowers=false;
            } else ntracks++;
        }

        // ---- Muon candidate (baseline: score>0.8, dist<20, len>10, pid>0.7) ----
        int muon_idx=-1;
        bool sel_muon=false, sel_muon_pid05=false, sel_muon_pid06=false;
        bool sel_muon_score05=false, sel_muon_score06=false, sel_muon_dist30=false;

        for (int i=0; i<(int)a.pfp_trk_shr_score->size(); i++) {
            TVector3 st(a.trk_start_x->at(i),a.trk_start_y->at(i),a.trk_start_z->at(i));
            double dist=(st-vtx).Mag();
            double score=a.pfp_trk_shr_score->at(i);
            double pid=a.trk_llrpid->at(i);
            double len=a.trk_length->at(i);
            if (score>0.8 && dist<20 && len>10 && pid>0.7) {
                sel_muon=true;
                if (muon_idx==-1||len>a.trk_length->at(muon_idx)) muon_idx=i;
            }
            if (score>0.8 && dist<20 && len>10 && pid>0.5) sel_muon_pid05=true;
            if (score>0.8 && dist<20 && len>10 && pid>0.6) sel_muon_pid06=true;
            if (score>0.5 && dist<20 && len>10 && pid>0.7) sel_muon_score05=true;
            if (score>0.6 && dist<20 && len>10 && pid>0.7) sel_muon_score06=true;
            if (score>0.8 && dist<30 && len>10 && pid>0.7) sel_muon_dist30=true;
        }
        // For loosened muon scenarios, re-pick the muon index with the looser cut
        // (for simplicity, re-run the index for each scenario separately below)

        // ---- Helper: pick candidates given muon_idx ----
        auto PickCandidates = [&](int mu_idx,
                                  int &pi_idx, int &pr_idx,
                                  int &sh1_idx, int &sh2_idx)
        {
            pi_idx=-1; pr_idx=-1; sh1_idx=-1; sh2_idx=-1;
            for (int i=0; i<(int)a.pfp_trk_shr_score->size(); i++) {
                if (i==mu_idx) continue;
                if (a.pfp_trk_shr_score->at(i)>0.5) {
                    if (pi_idx==-1 || a.trk_llrpid->at(i)>a.trk_llrpid->at(pi_idx))
                        pi_idx=i;
                }
            }
            for (int i=0; i<(int)a.pfp_trk_shr_score->size(); i++) {
                if (i==mu_idx||i==pi_idx) continue;
                if (a.pfp_trk_shr_score->at(i)>0.5) {
                    if (pr_idx==-1 || a.trk_llrpid->at(i)<a.trk_llrpid->at(pr_idx))
                        pr_idx=i;
                }
            }
            double b1=1e9, b2=1e9;
            for (int i=0; i<(int)a.pfp_trk_shr_score->size(); i++) {
                if (i==mu_idx||i==pi_idx||i==pr_idx) continue;
                if (a.pfp_trk_shr_score->at(i)>=0.5) continue;
                TVector3 st(a.trk_start_x->at(i),a.trk_start_y->at(i),a.trk_start_z->at(i));
                double d=(st-vtx).Mag();
                if (d<b1) { b2=b1; sh2_idx=sh1_idx; b1=d; sh1_idx=i; }
                else if (d<b2) { b2=d; sh2_idx=i; }
            }
        };

        // Baseline candidate indices
        int pion_idx, proton_idx, shr1_idx, shr2_idx;
        PickCandidates(muon_idx, pion_idx, proton_idx, shr1_idx, shr2_idx);

        bool sel_pion   = (pion_idx   >= 0);
        bool sel_proton = (proton_idx >= 0);

        // ---- Compute directions & distances ----
        TVector3 muDir, piDir, prDir, sh1Dir, sh2Dir;
        TVector3 muSt, piSt, prSt, sh1St, sh2St;
        if (muon_idx  >-1){ muSt.SetXYZ(a.trk_start_x->at(muon_idx),a.trk_start_y->at(muon_idx),a.trk_start_z->at(muon_idx));
                             muDir.SetXYZ(a.trk_dir_x->at(muon_idx),a.trk_dir_y->at(muon_idx),a.trk_dir_z->at(muon_idx)); }
        if (pion_idx  >-1){ piSt.SetXYZ(a.trk_start_x->at(pion_idx),a.trk_start_y->at(pion_idx),a.trk_start_z->at(pion_idx));
                             piDir.SetXYZ(a.trk_dir_x->at(pion_idx),a.trk_dir_y->at(pion_idx),a.trk_dir_z->at(pion_idx)); }
        if (proton_idx>-1){ prSt.SetXYZ(a.trk_start_x->at(proton_idx),a.trk_start_y->at(proton_idx),a.trk_start_z->at(proton_idx));
                             prDir.SetXYZ(a.trk_dir_x->at(proton_idx),a.trk_dir_y->at(proton_idx),a.trk_dir_z->at(proton_idx)); }
        if (shr1_idx  >-1){ sh1St.SetXYZ(a.trk_start_x->at(shr1_idx),a.trk_start_y->at(shr1_idx),a.trk_start_z->at(shr1_idx));
        if (shr2_idx  >-1)  sh2St.SetXYZ(a.trk_start_x->at(shr2_idx),a.trk_start_y->at(shr2_idx),a.trk_start_z->at(shr2_idx));
                             sh1Dir.SetXYZ(a.trk_dir_x->at(shr1_idx),a.trk_dir_y->at(shr1_idx),a.trk_dir_z->at(shr1_idx)); }

        double dMuPi  = (muon_idx>-1&&pion_idx>-1)   ? (muSt-piSt).Mag()  : 99999;
        double dPrPi  = (proton_idx>-1&&pion_idx>-1) ? (prSt-piSt).Mag()  : 99999;
        double dPrMu  = (proton_idx>-1&&muon_idx>-1) ? (prSt-muSt).Mag()  : 99999;
        double dS1Mu  = (shr1_idx>-1&&muon_idx>-1)   ? (sh1St-muSt).Mag() : 99999;
        double dS1Pi  = (shr1_idx>-1&&pion_idx>-1)   ? (sh1St-piSt).Mag() : 99999;
        double dS1Pr  = (shr1_idx>-1&&proton_idx>-1) ? (sh1St-prSt).Mag() : 99999;

        // ---- Shower quality flags ----
        bool sel_2nd_shr = true;
        if (shr2_idx>-1 && a.pfpnplanehits_Y->at(shr2_idx)>100) sel_2nd_shr=false;

        bool sel_1st_shr = false;
        if (shr1_idx>-1) {
            int h=a.pfpnplanehits_Y->at(shr1_idx);
            if (h>10 && h<300) sel_1st_shr=true;
        }

        // ---- Standard distance cut ----
        bool sel_dist = (dMuPi<30 && dPrMu<20 && dPrPi<20 &&
                         dS1Mu<100 && dS1Pi<100 && dS1Pr<100);

        // ---- Fill reader variables (for BDTg evaluation) ----
        r_MuonTrkLength  = (muon_idx>-1)    ? SafeGet(a.trk_length,            muon_idx)  : -999.f;
        r_MuonTrkDistVtx = (muon_idx>-1)    ? (Float_t)(muSt-vtx).Mag()                  : -999.f;
        r_MuonTrkScore   = (muon_idx>-1)    ? SafeGet(a.pfp_trk_shr_score,    muon_idx)  : -999.f;
        r_MuonNHitsY     = (muon_idx>-1)    ? (Float_t)SafeGetI(a.pfpnplanehits_Y,muon_idx) : -999.f;
        r_MuonPhotonPID  = (muon_idx>-1)    ? SafeGet(a.pfp_photon_bdt_score,  muon_idx)  : -999.f;

        r_PionPID        = (pion_idx>-1)    ? SafeGet(a.trk_llrpid,              pion_idx)  : -999.f;
        r_PionTrkLength  = (pion_idx>-1)    ? SafeGet(a.trk_length,               pion_idx)  : -999.f;
        r_PionTrkDistVtx = (pion_idx>-1)    ? (Float_t)(piSt-vtx).Mag()                      : -999.f;
        r_PionTrkScore   = (pion_idx>-1)    ? SafeGet(a.pfp_trk_shr_score,       pion_idx)  : -999.f;
        r_PionNHitsY     = (pion_idx>-1)    ? (Float_t)SafeGetI(a.pfpnplanehits_Y, pion_idx) : -999.f;
        r_PionDeDxP2     = (pion_idx>-1)    ? SafeGet(a.trk_mean_dedx_plane2,    pion_idx)  : -999.f;
        r_PionPhotonPID  = (pion_idx>-1)    ? SafeGet(a.pfp_photon_bdt_score,     pion_idx)  : -999.f;

        r_ProtonPID        = (proton_idx>-1) ? SafeGet(a.trk_llrpid,               proton_idx) : -999.f;
        r_ProtonTrkLength  = (proton_idx>-1) ? SafeGet(a.trk_length,                proton_idx) : -999.f;
        r_ProtonTrkDistVtx = (proton_idx>-1) ? (Float_t)(prSt-vtx).Mag()                        : -999.f;
        r_ProtonTrkScore   = (proton_idx>-1) ? SafeGet(a.pfp_trk_shr_score,        proton_idx) : -999.f;
        r_ProtonNHitsY     = (proton_idx>-1) ? (Float_t)SafeGetI(a.pfpnplanehits_Y, proton_idx) : -999.f;
        r_ProtonDeDxP2     = (proton_idx>-1) ? SafeGet(a.trk_mean_dedx_plane2,     proton_idx) : -999.f;
        r_ProtonPhotonPID  = (proton_idx>-1) ? SafeGet(a.pfp_photon_bdt_score,      proton_idx) : -999.f;

        r_Shr1NHitsY    = (shr1_idx>-1) ? (Float_t)SafeGetI(a.pfpnplanehits_Y,shr1_idx)  : -999.f;
        r_Shr1VtxDist   = (shr1_idx>-1) ? (Float_t)(sh1St-vtx).Mag()                     : -999.f;
        r_Shr1OpenAngle = (shr1_idx>-1) ? SafeGet(a.shr_open_angle,           shr1_idx)   : -999.f;
        r_Shr2NHitsY    = (shr2_idx>-1) ? (Float_t)SafeGetI(a.pfpnplanehits_Y,shr2_idx)  : -999.f;
        r_Shr2VtxDist   = (shr2_idx>-1) ? (Float_t)(sh2St-vtx).Mag()                     : -999.f;

        r_AngleMuonPion    = (muon_idx>-1   &&pion_idx>-1)   ? (Float_t)muDir.Angle(piDir)  : -999.f;
        r_AngleProtonPion  = (proton_idx>-1 &&pion_idx>-1)   ? (Float_t)prDir.Angle(piDir)  : -999.f;
        r_AngleProtonMuon  = (proton_idx>-1 &&muon_idx>-1)   ? (Float_t)prDir.Angle(muDir)  : -999.f;
        r_AngleShr1Muon    = (shr1_idx>-1  &&muon_idx>-1)    ? (Float_t)sh1Dir.Angle(muDir) : -999.f;
        r_AngleShr1Pion    = (shr1_idx>-1  &&pion_idx>-1)    ? (Float_t)sh1Dir.Angle(piDir) : -999.f;
        r_AngleShr1Proton  = (shr1_idx>-1  &&proton_idx>-1)  ? (Float_t)sh1Dir.Angle(prDir) : -999.f;

        r_DistMuonPion   = (Float_t)dMuPi;
        r_DistProtonPion = (Float_t)dPrPi;
        r_DistProtonMuon = (Float_t)dPrMu;
        r_Weight         = (Float_t)w;

        // ---- MVA scores (only meaningful when all candidates exist) ----
        float bdtg = -999.f;
        float bdt  = -999.f;
        float bdta = -999.f;
        bool have_candidates = (muon_idx>-1 && pion_idx>-1 && proton_idx>-1 && shr1_idx>-1);
        if (have_candidates) {
            bdt  = rd->EvaluateMVA("BDT");
            bdta = rd->EvaluateMVA("BDTAda");
            bdtg = rd->EvaluateMVA("BDTg");
        }
        bool passMVA = (bdt>0.84f && bdta>0.0f && bdtg>0.96f);

        // ---- Accumulate for each scenario ----
        auto add = [&](int idx, bool presel) {
            if (!presel) return;
            if (!passMVA) return;
            switch(s) {
                case 0: sc[idx].S   += w; break;
                case 1: sc[idx].Hyp += w; break;
                case 2: sc[idx].Bkg += w; break;
                case 3: sc[idx].Dirt+= w; break;
            }
        };

        // Baseline — exact reproduction of Sigma.C
        bool base_pre = (nshowers>0 && nshowers<3) && NoBigShowers
                     && sel_muon && sel_pion && sel_proton
                     && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist;
        add(kBase, base_pre);

        // Skip NShowers cut
        add(kNoNShowers, NoBigShowers && sel_muon && sel_pion && sel_proton
                      && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Skip NoBigShowers cut
        add(kNoNoBigShr, (nshowers>0&&nshowers<3) && sel_muon && sel_pion && sel_proton
                      && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Skip NTracksCut
        add(kNoNTracks, (nshowers>0&&nshowers<3) && NoBigShowers
                     && sel_muon && sel_pion && sel_proton
                     && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Skip 2ndShower veto
        add(kNo2ndShr, (nshowers>0&&nshowers<3) && NoBigShowers
                    && sel_muon && sel_pion && sel_proton
                    && (ntracks<5) && sel_1st_shr && sel_dist);

        // Skip 1stShower cut
        add(kNo1stShr, (nshowers>0&&nshowers<3) && NoBigShowers
                    && sel_muon && sel_pion && sel_proton
                    && (ntracks<5) && sel_2nd_shr && sel_dist);

        // Skip Distances cut
        add(kNoDist, (nshowers>0&&nshowers<3) && NoBigShowers
                  && sel_muon && sel_pion && sel_proton
                  && (ntracks<5) && sel_2nd_shr && sel_1st_shr);

        // Skip MuonCandidate
        add(kNoMuon, (nshowers>0&&nshowers<3) && NoBigShowers
                  && sel_pion && sel_proton
                  && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Skip PionCandidate
        add(kNoPion, (nshowers>0&&nshowers<3) && NoBigShowers
                  && sel_muon && sel_proton
                  && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Skip ProtonCandidate
        add(kNoProton, (nshowers>0&&nshowers<3) && NoBigShowers
                    && sel_muon && sel_pion
                    && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Looser muon PID thresholds
        add(kMuPID05, (nshowers>0&&nshowers<3) && NoBigShowers
                   && sel_muon_pid05 && sel_pion && sel_proton
                   && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);
        add(kMuPID06, (nshowers>0&&nshowers<3) && NoBigShowers
                   && sel_muon_pid06 && sel_pion && sel_proton
                   && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Looser muon track score
        add(kMuScore05, (nshowers>0&&nshowers<3) && NoBigShowers
                     && sel_muon_score05 && sel_pion && sel_proton
                     && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);
        add(kMuScore06, (nshowers>0&&nshowers<3) && NoBigShowers
                     && sel_muon_score06 && sel_pion && sel_proton
                     && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Looser muon distance to vertex
        add(kMuDist30, (nshowers>0&&nshowers<3) && NoBigShowers
                    && sel_muon_dist30 && sel_pion && sel_proton
                    && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Allow 3 showers
        add(kShr3, (nshowers>0&&nshowers<4) && NoBigShowers
                && sel_muon && sel_pion && sel_proton
                && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Looser NTracks cuts
        add(kNTracks5, (nshowers>0&&nshowers<3) && NoBigShowers
                    && sel_muon && sel_pion && sel_proton
                    && (ntracks<6) && sel_2nd_shr && sel_1st_shr && sel_dist);
        add(kNTracks6, (nshowers>0&&nshowers<3) && NoBigShowers
                    && sel_muon && sel_pion && sel_proton
                    && (ntracks<7) && sel_2nd_shr && sel_1st_shr && sel_dist);

        // Wider 1st shower hit window
        bool sel_1st_shr_max400 = shr1_idx>-1 && a.pfpnplanehits_Y->at(shr1_idx)>10
                                               && a.pfpnplanehits_Y->at(shr1_idx)<400;
        bool sel_1st_shr_min5   = shr1_idx>-1 && a.pfpnplanehits_Y->at(shr1_idx)>5
                                               && a.pfpnplanehits_Y->at(shr1_idx)<300;
        add(kShr1HitsMax400, (nshowers>0&&nshowers<3) && NoBigShowers
                          && sel_muon && sel_pion && sel_proton
                          && (ntracks<5) && sel_2nd_shr && sel_1st_shr_max400 && sel_dist);
        add(kShr1HitsMin5,   (nshowers>0&&nshowers<3) && NoBigShowers
                          && sel_muon && sel_pion && sel_proton
                          && (ntracks<5) && sel_2nd_shr && sel_1st_shr_min5 && sel_dist);

        // Looser muon-pion distance
        bool sel_dist40 = (dMuPi<40 && dPrMu<20 && dPrPi<20 &&
                           dS1Mu<100 && dS1Pi<100 && dS1Pr<100);
        bool sel_dist50 = (dMuPi<50 && dPrMu<20 && dPrPi<20 &&
                           dS1Mu<100 && dS1Pi<100 && dS1Pr<100);
        add(kDistMuPi40, (nshowers>0&&nshowers<3) && NoBigShowers
                      && sel_muon && sel_pion && sel_proton
                      && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist40);
        add(kDistMuPi50, (nshowers>0&&nshowers<3) && NoBigShowers
                      && sel_muon && sel_pion && sel_proton
                      && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist50);

        // Looser proton-X distance
        bool sel_dist_pr30 = (dMuPi<30 && dPrMu<30 && dPrPi<30 &&
                              dS1Mu<100 && dS1Pi<100 && dS1Pr<100);
        add(kDistProt30, (nshowers>0&&nshowers<3) && NoBigShowers
                      && sel_muon && sel_pion && sel_proton
                      && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist_pr30);

        // Looser shower-candidate distance
        bool sel_dist_shr150 = (dMuPi<30 && dPrMu<20 && dPrPi<20 &&
                                dS1Mu<150 && dS1Pi<150 && dS1Pr<150);
        add(kDistShr150, (nshowers>0&&nshowers<3) && NoBigShowers
                      && sel_muon && sel_pion && sel_proton
                      && (ntracks<5) && sel_2nd_shr && sel_1st_shr && sel_dist_shr150);

    } // end event loop
    cout << "\r100%\n";

    // -----------------------------------------------------------
    // Print results
    // -----------------------------------------------------------
    double S0 = sc[kBase].S;  // baseline signal for relative efficiency

    printf("\n%-28s  %7s  %7s  %7s  %7s  %7s  %6s  %7s  %5s\n",
           "Scenario","S","Hyp","Bkg","Dirt","B_tot","Eff%","Purity%","S/sqB");
    printf("%s\n", string(100,'-').c_str());
    for (int i=0; i<kNScenarios; i++) {
        printf("%-28s  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %6.2f  %7.2f  %5.3f",
               sc[i].name.c_str(),
               sc[i].S, sc[i].Hyp, sc[i].Bkg, sc[i].Dirt, sc[i].B(),
               sc[i].eff(S0), sc[i].pur(), sc[i].fom());
        if (i==kBase) printf("  <- BASELINE");
        else {
            double dS = sc[i].S - S0;
            if (dS > 0.005) printf("  +%.3f sig", dS);
        }
        printf("\n");
    }
    printf("\nNote: MVA was trained on the baseline preselection.\n"
           "Scenarios that feed the MVA with out-of-distribution events\n"
           "may need retraining for a production analysis.\n");
}
