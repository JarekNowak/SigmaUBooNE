// ============================================================
//  PlotEfficiency.C
//  Selection efficiency as a function of true particle momentum
//  for all six final-state particles: mu+, Sigma0, Lambda,
//  photon, pi-, proton.
//
//  Two selection levels are shown per particle:
//    Pre-MVA : after all pre-selection + distance cuts (Cuts[10])
//    Final   : BDTg > 0.96 (Cuts[13])
//
//  Denominator: true signal events in FV
//    (anti-nu_mu, Sigma0->Lambda(->p pi-) + gamma, true vtx in FV)
//
//  Usage:
//    cd /home/t2k/nowak/MicroBooNE/Sigma
//    root -b -l -q PlotEfficiency.C
// ============================================================

#define Sigma_cxx
#include "Sigma.h"
#include "TMVA/Reader.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

using namespace std;

// SafeGet helpers (identical to Sigma.C)
template<typename T>
static inline float SafeGet(const vector<T>* v, int i, float s=-999.f)
{ return (!v||i<0||i>=(int)v->size()) ? s : (float)v->at(i); }

template<typename T>
static inline int SafeGetI(const vector<T>* v, int i, int s=-999)
{ return (!v||i<0||i>=(int)v->size()) ? s : (int)v->at(i); }

// ============================================================
void PlotEfficiency()
// ============================================================
{
  gStyle->SetOptStat(0);

  // ----------------------------------------------------------
  // Open signal (Hyperon) MC only
  // ----------------------------------------------------------
  TFile *inFile = TFile::Open(
    "/data/nowak/Sigma/"
    "analysisOutputRHC_cttest_Hyperon_ana_hyperon_sigmazero_mSlice_prodgenie_run3b.root",
    "READ");
  if (!inFile || inFile->IsZombie()) {
    cerr << "ERROR: Cannot open Hyperon MC file." << endl; return;
  }
  TTree *tree = (TTree*)inFile->Get("OutputTree");
  if (!tree) { cerr << "ERROR: OutputTree not found." << endl; return; }

  Sigma sigma(tree);

  // ----------------------------------------------------------
  // TMVA Reader — identical variable list to Sigma.C
  // ----------------------------------------------------------
  Float_t r_MuonTrkLength, r_MuonTrkDistVtx, r_MuonTrkScore;
  Float_t r_MuonNHitsY,    r_MuonPhotonPID;
  Float_t r_PionPID,       r_PionTrkLength,    r_PionTrkDistVtx, r_PionTrkScore;
  Float_t r_PionNHitsY,    r_PionDeDxP2,       r_PionPhotonPID;
  Float_t r_ProtonPID,     r_ProtonTrkLength,  r_ProtonTrkDistVtx, r_ProtonTrkScore;
  Float_t r_ProtonNHitsY,  r_ProtonDeDxP2,     r_ProtonPhotonPID;
  Float_t r_Shr1NHitsY,    r_Shr1VtxDist,      r_Shr1OpenAngle;
  Float_t r_Shr2NHitsY,    r_Shr2VtxDist;
  Float_t r_AngleMuonPion, r_AngleProtonPion,  r_AngleProtonMuon;
  Float_t r_AngleShr1Muon, r_AngleShr1Pion,    r_AngleShr1Proton;
  Float_t r_DistMuonPion,  r_DistProtonPion,   r_DistProtonMuon;
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
  mvaReader->AddSpectator("Weight",          &r_Weight);
  mvaReader->BookMVA("BDTg", "SigmaDataset/weights/TMVASigma_BDTg.weights.xml");

  // ----------------------------------------------------------
  // TEfficiency objects: 6 particles × 2 levels
  // Bayesian 68 % CL intervals, flat prior
  // ----------------------------------------------------------
  // TEfficiency title format ";x-title;y-title" — sets axis labels at construction
  auto mkE = [](const char* n, const char* xtitle, int nb, double lo, double hi)
  {
    TString t = TString(";") + xtitle + ";Efficiency";
    TEfficiency *e = new TEfficiency(n, t.Data(), nb, lo, hi);
    e->SetStatisticOption(TEfficiency::kBUniform);
    e->SetConfidenceLevel(0.683);
    return e;
  };

  TEfficiency *ePre[6], *eFin[6];
  ePre[0] = mkE("ePre_mu",  "True p_{#mu^{+}} [GeV/c]",    12, 0.0, 2.4);
  ePre[1] = mkE("ePre_sig", "True p_{#Sigma^{0}} [GeV/c]", 10, 0.0, 2.0);
  ePre[2] = mkE("ePre_lam", "True p_{#Lambda} [GeV/c]",    10, 0.0, 2.0);
  ePre[3] = mkE("ePre_gam", "True p_{#gamma} [GeV/c]",     10, 0.0, 0.5);
  ePre[4] = mkE("ePre_pi",  "True p_{#pi^{-}} [GeV/c]",    10, 0.0, 1.0);
  ePre[5] = mkE("ePre_pr",  "True p_{p} [GeV/c]",          10, 0.0, 1.5);
  eFin[0] = mkE("eFin_mu",  "True p_{#mu^{+}} [GeV/c]",    12, 0.0, 2.4);
  eFin[1] = mkE("eFin_sig", "True p_{#Sigma^{0}} [GeV/c]", 10, 0.0, 2.0);
  eFin[2] = mkE("eFin_lam", "True p_{#Lambda} [GeV/c]",    10, 0.0, 2.0);
  eFin[3] = mkE("eFin_gam", "True p_{#gamma} [GeV/c]",     10, 0.0, 0.5);
  eFin[4] = mkE("eFin_pi",  "True p_{#pi^{-}} [GeV/c]",    10, 0.0, 1.0);
  eFin[5] = mkE("eFin_pr",  "True p_{p} [GeV/c]",          10, 0.0, 1.5);

  for (int k = 0; k < 6; k++) {
    ePre[k]->SetLineColor(kBlue+1);   ePre[k]->SetMarkerColor(kBlue+1);
    ePre[k]->SetMarkerStyle(24);       ePre[k]->SetMarkerSize(1.1);
    eFin[k]->SetLineColor(kRed+1);    eFin[k]->SetMarkerColor(kRed+1);
    eFin[k]->SetMarkerStyle(20);       eFin[k]->SetMarkerSize(1.1);
  }

  // ----------------------------------------------------------
  // Main event loop
  // ----------------------------------------------------------
  Long64_t nentries = tree->GetEntries();
  cout << "Looping over " << nentries << " entries..." << endl;

  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    sigma.LoadTree(jentry);
    sigma.fChain->GetEntry(jentry);
    if (jentry % 50000 == 0)
      cout << "\r  " << (100*(jentry+1)/nentries) << " %" << flush;

    // ---- Truth-level decay products ----
    double proton_mom = 0, pion_mom = 0;
    int mc_decay = 0;
    for (int i = 0; i < (int)sigma.mc_decay_pdg->size(); i++) {
      if (sigma.mc_decay_pdg->at(i) ==  2212) { mc_decay++; proton_mom = SafeGet(sigma.mc_decay_mom, i); }
      if (sigma.mc_decay_pdg->at(i) == -211)  { mc_decay++; pion_mom   = SafeGet(sigma.mc_decay_mom, i); }
    }

    bool inFV = (sigma.mc_nu_pos_x >   3.00 && sigma.mc_nu_pos_x < 253.35 &&
                 sigma.mc_nu_pos_y > -112.53 && sigma.mc_nu_pos_y < 114.47 &&
                 sigma.mc_nu_pos_z >   3.10  && sigma.mc_nu_pos_z < 1033.0);
    bool dead  = (sigma.mc_nu_pos_z > 675.1  && sigma.mc_nu_pos_z <  775.1);

    bool IsSignal = (sigma.mc_nu_pdg == -14) && (sigma.mc_hyperon_pdg == 3212)
                 && (mc_decay == 2) && proton_mom > 0.3 && pion_mom > 0.1;

    if (!IsSignal || !inFV || dead) continue;  // denominator gate

    // True momenta for efficiency axes
    double p[6] = {
      sigma.mc_lepton_mom,             // mu+
      sigma.mc_hyperon_mom,            // Sigma0
      sigma.mc_sigmazero_lambda_mom,   // Lambda
      sigma.mc_sigmazero_photon_mom,   // photon
      pion_mom,                        // pi-
      proton_mom                       // proton
    };

    // ---- Reco-level cuts (identical to Sigma.C) ----
    TVector3 vtx(sigma.reco_primary_vtx_x, sigma.reco_primary_vtx_y, sigma.reco_primary_vtx_z);
    bool sel_FV    = sigma.reco_primary_vtx_inFV;
    int  ntracks   = 0, nshowers = 0;
    int  iMu=-1, iPi=-1, iPr=-1, iSh1=-1, iSh2=-1;
    bool sel_muon=false, sel_pion=false, sel_proton=false;

    for (int i = 0; i < (int)sigma.pfp_trk_shr_score->size(); i++) {
      if (sigma.pfp_trk_shr_score->at(i) < 0.5) nshowers++;
      else                                        ntracks++;
    }

    // Muon candidate
    for (int i = 0; i < (int)sigma.pfp_trk_shr_score->size(); i++) {
      TVector3 st(sigma.trk_start_x->at(i), sigma.trk_start_y->at(i), sigma.trk_start_z->at(i));
      if (sigma.pfp_trk_shr_score->at(i) > 0.8 && (st-vtx).Mag() < 20
          && sigma.trk_length->at(i) > 10 && sigma.trk_llrpid->at(i) > 0.7) {
        sel_muon = true;
        if (iMu==-1 || sigma.trk_length->at(i) > sigma.trk_length->at(iMu)) iMu = i;
      }
    }
    // Pion candidate
    for (int i = 0; i < (int)sigma.pfp_trk_shr_score->size(); i++) {
      if (i==iMu) continue;
      if (sigma.pfp_trk_shr_score->at(i) > 0.5) {
        sel_pion = true;
        if (iPi==-1 || sigma.trk_llrpid->at(i) > sigma.trk_llrpid->at(iPi)) iPi = i;
      }
    }
    // Proton candidate
    for (int i = 0; i < (int)sigma.pfp_trk_shr_score->size(); i++) {
      if (i==iMu || i==iPi) continue;
      if (sigma.pfp_trk_shr_score->at(i) > 0.5) {
        sel_proton = true;
        if (iPr==-1 || sigma.trk_llrpid->at(i) < sigma.trk_llrpid->at(iPr)) iPr = i;
      }
    }
    // Shower candidates
    double sb1=1e9, sb2=1e9;
    for (int i = 0; i < (int)sigma.pfp_trk_shr_score->size(); i++) {
      if (i==iMu||i==iPi||i==iPr) continue;
      if (sigma.pfp_trk_shr_score->at(i) >= 0.5) continue;
      TVector3 st(sigma.trk_start_x->at(i), sigma.trk_start_y->at(i), sigma.trk_start_z->at(i));
      double d = (st-vtx).Mag();
      if (d<sb1) { sb2=sb1; iSh2=iSh1; sb1=d; iSh1=i; }
      else if (d<sb2) { sb2=d; iSh2=i; }
    }

    // Shower quality and track count cuts
    bool sel_1stShr = false;
    if (iSh1>-1) { int h=sigma.pfpnplanehits_Y->at(iSh1); sel_1stShr=(h>10&&h<300); }
    bool sel_nShr  = (nshowers>0 && nshowers<4);
    bool sel_nTrk  = (ntracks<5);

    // Inter-candidate positions and directions
    TVector3 muSt,piSt,prSt,s1St,s2St;
    TVector3 muDir,piDir,prDir,s1Dir;
    if (iMu>-1)  { muSt.SetXYZ(sigma.trk_start_x->at(iMu),  sigma.trk_start_y->at(iMu),  sigma.trk_start_z->at(iMu));
                   muDir.SetXYZ(sigma.trk_dir_x->at(iMu),   sigma.trk_dir_y->at(iMu),   sigma.trk_dir_z->at(iMu)); }
    if (iPi>-1)  { piSt.SetXYZ(sigma.trk_start_x->at(iPi),  sigma.trk_start_y->at(iPi),  sigma.trk_start_z->at(iPi));
                   piDir.SetXYZ(sigma.trk_dir_x->at(iPi),   sigma.trk_dir_y->at(iPi),   sigma.trk_dir_z->at(iPi)); }
    if (iPr>-1)  { prSt.SetXYZ(sigma.trk_start_x->at(iPr),  sigma.trk_start_y->at(iPr),  sigma.trk_start_z->at(iPr));
                   prDir.SetXYZ(sigma.trk_dir_x->at(iPr),   sigma.trk_dir_y->at(iPr),   sigma.trk_dir_z->at(iPr)); }
    if (iSh1>-1) { s1St.SetXYZ(sigma.trk_start_x->at(iSh1), sigma.trk_start_y->at(iSh1), sigma.trk_start_z->at(iSh1));
                   s1Dir.SetXYZ(sigma.trk_dir_x->at(iSh1),  sigma.trk_dir_y->at(iSh1),  sigma.trk_dir_z->at(iSh1)); }
    if (iSh2>-1) { s2St.SetXYZ(sigma.trk_start_x->at(iSh2), sigma.trk_start_y->at(iSh2), sigma.trk_start_z->at(iSh2)); }

    double dMuPi  = (iMu>-1&&iPi>-1)   ? (muSt-piSt).Mag()  : 99999;
    double dPrPi  = (iPr>-1&&iPi>-1)   ? (prSt-piSt).Mag()  : 99999;
    double dPrMu  = (iPr>-1&&iMu>-1)   ? (prSt-muSt).Mag()  : 99999;
    double dS1Mu  = (iSh1>-1&&iMu>-1)  ? (s1St-muSt).Mag()  : 99999;
    double dS1Pi  = (iSh1>-1&&iPi>-1)  ? (s1St-piSt).Mag()  : 99999;
    double dS1Pr  = (iSh1>-1&&iPr>-1)  ? (s1St-prSt).Mag()  : 99999;

    bool sel_dist = (dMuPi<30 && dPrMu<20 && dPrPi<20 && dS1Mu<100 && dS1Pi<100 && dS1Pr<100);

    // Cumulative cuts
    bool Cuts[14] = {false};
    Cuts[0]  = true;
    Cuts[1]  = Cuts[0] && sel_FV;
    Cuts[2]  = Cuts[1] && sel_nShr;
    Cuts[3]  = Cuts[2];
    Cuts[4]  = Cuts[3] && sel_muon;
    Cuts[5]  = Cuts[4] && sel_pion;
    Cuts[6]  = Cuts[5] && sel_proton;
    Cuts[7]  = Cuts[6] && sel_nTrk;
    Cuts[8]  = Cuts[7];
    Cuts[9]  = Cuts[8] && sel_1stShr;
    Cuts[10] = Cuts[9] && sel_dist;

    // BDTg evaluation
    float mva_BDTg = -999.f;
    if (Cuts[10]) {
      r_MuonTrkLength    = (iMu>-1)  ? SafeGet(sigma.trk_length,            iMu)                      : -999.f;
      r_MuonTrkDistVtx   = (iMu>-1)  ? (Float_t)(muSt-vtx).Mag()                                     : -999.f;
      r_MuonTrkScore     = (iMu>-1)  ? SafeGet(sigma.pfp_trk_shr_score,     iMu)                      : -999.f;
      r_MuonNHitsY       = (iMu>-1)  ? (Float_t)SafeGetI(sigma.pfpnplanehits_Y, iMu)                  : -999.f;
      r_MuonPhotonPID    = (iMu>-1)  ? SafeGet(sigma.pfp_photon_bdt_score,  iMu)                      : -999.f;
      r_PionPID          = (iPi>-1)  ? SafeGet(sigma.trk_llrpid,            iPi)                      : -999.f;
      r_PionTrkLength    = (iPi>-1)  ? SafeGet(sigma.trk_length,            iPi)                      : -999.f;
      r_PionTrkDistVtx   = (iPi>-1)  ? (Float_t)(piSt-vtx).Mag()                                     : -999.f;
      r_PionTrkScore     = (iPi>-1)  ? SafeGet(sigma.pfp_trk_shr_score,     iPi)                      : -999.f;
      r_PionNHitsY       = (iPi>-1)  ? (Float_t)SafeGetI(sigma.pfpnplanehits_Y, iPi)                  : -999.f;
      r_PionDeDxP2       = (iPi>-1)  ? SafeGet(sigma.trk_mean_dedx_plane2,  iPi)                      : -999.f;
      r_PionPhotonPID    = (iPi>-1)  ? SafeGet(sigma.pfp_photon_bdt_score,  iPi)                      : -999.f;
      r_ProtonPID        = (iPr>-1)  ? SafeGet(sigma.trk_llrpid,            iPr)                      : -999.f;
      r_ProtonTrkLength  = (iPr>-1)  ? SafeGet(sigma.trk_length,            iPr)                      : -999.f;
      r_ProtonTrkDistVtx = (iPr>-1)  ? (Float_t)(prSt-vtx).Mag()                                     : -999.f;
      r_ProtonTrkScore   = (iPr>-1)  ? SafeGet(sigma.pfp_trk_shr_score,     iPr)                      : -999.f;
      r_ProtonNHitsY     = (iPr>-1)  ? (Float_t)SafeGetI(sigma.pfpnplanehits_Y, iPr)                  : -999.f;
      r_ProtonDeDxP2     = (iPr>-1)  ? SafeGet(sigma.trk_mean_dedx_plane2,  iPr)                      : -999.f;
      r_ProtonPhotonPID  = (iPr>-1)  ? SafeGet(sigma.pfp_photon_bdt_score,  iPr)                      : -999.f;
      r_Shr1NHitsY       = (iSh1>-1) ? (Float_t)SafeGetI(sigma.pfpnplanehits_Y, iSh1)                 : -999.f;
      r_Shr1VtxDist      = (iSh1>-1) ? (Float_t)(s1St-vtx).Mag()                                     : -999.f;
      r_Shr1OpenAngle    = (iSh1>-1) ? SafeGet(sigma.shr_open_angle,       iSh1)                      : -999.f;
      r_Shr2NHitsY       = (iSh2>-1) ? (Float_t)SafeGetI(sigma.pfpnplanehits_Y, iSh2)                 : -999.f;
      r_Shr2VtxDist      = (iSh2>-1) ? (Float_t)(s2St-vtx).Mag()                                     : -999.f;
      r_AngleMuonPion    = (iMu>-1&&iPi>-1)   ? (Float_t)muDir.Angle(piDir)    : -999.f;
      r_AngleProtonPion  = (iPr>-1&&iPi>-1)   ? (Float_t)prDir.Angle(piDir)    : -999.f;
      r_AngleProtonMuon  = (iPr>-1&&iMu>-1)   ? (Float_t)prDir.Angle(muDir)    : -999.f;
      r_AngleShr1Muon    = (iSh1>-1&&iMu>-1)  ? (Float_t)s1Dir.Angle(muDir)    : -999.f;
      r_AngleShr1Pion    = (iSh1>-1&&iPi>-1)  ? (Float_t)s1Dir.Angle(piDir)    : -999.f;
      r_AngleShr1Proton  = (iSh1>-1&&iPr>-1)  ? (Float_t)s1Dir.Angle(prDir)    : -999.f;
      r_DistMuonPion     = (Float_t)dMuPi;
      r_DistProtonPion   = (Float_t)dPrPi;
      r_DistProtonMuon   = (Float_t)dPrMu;
      r_Weight           = 1.0f;
      mva_BDTg = mvaReader->EvaluateMVA("BDTg");
    }

    Cuts[11] = Cuts[10];
    Cuts[12] = Cuts[11];
    Cuts[13] = Cuts[12] && (mva_BDTg > 0.96f);

    bool passP = Cuts[10];   // pre-MVA
    bool passF = Cuts[13];   // final

    for (int k = 0; k < 6; k++) {
      ePre[k]->Fill(passP, p[k]);
      eFin[k]->Fill(passF, p[k]);
    }
  }
  cout << "\nLoop done." << endl;

  // ----------------------------------------------------------
  // Output directory
  // ----------------------------------------------------------
  gSystem->mkdir("figs/EfficiencyPlots", kTRUE);

  const char* labels[6] = {
    "#mu^{+}",
    "#Sigma^{0}",
    "#Lambda",
    "#gamma  (from #Sigma^{0}#rightarrow#Lambda#gamma)",
    "#pi^{-} (from #Lambda#rightarrowp#pi^{-})",
    "p       (from #Lambda#rightarrowp#pi^{-})"
  };
  const char* fnames[6] = {
    "Eff_Muon", "Eff_Sigma0", "Eff_Lambda",
    "Eff_Photon", "Eff_Pion", "Eff_Proton"
  };
  double xmax[6] = {2.4, 2.0, 2.0, 0.5, 1.0, 1.5};

  // ----------------------------------------------------------
  // Individual plots
  // ----------------------------------------------------------
  for (int k = 0; k < 6; k++) {
    TCanvas *c = new TCanvas(fnames[k], fnames[k], 700, 500);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
    gPad->SetGrid();

    ePre[k]->Draw("AP");
    c->Update();
    TGraphAsymmErrors *gPre = ePre[k]->GetPaintedGraph();
    gPre->GetYaxis()->SetRangeUser(0.0, 0.50);
    gPre->GetXaxis()->SetLimits(0.0, xmax[k]);
    gPre->GetXaxis()->SetTitleSize(0.05);
    gPre->GetYaxis()->SetTitleSize(0.05);
    gPre->GetXaxis()->SetTitleOffset(1.1);
    gPre->GetYaxis()->SetTitleOffset(1.1);
    c->Modified(); c->Update();
    eFin[k]->Draw("P SAME");

    TLatex lat;
    lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.033);
    lat.DrawLatex(0.15, 0.91,
      "MicroBooNE BNB run 3b,  4.9#times10^{20} PoT  (signal MC)");
    lat.DrawLatex(0.15, 0.85, labels[k]);

    TLegend *leg = new TLegend(0.50, 0.68, 0.87, 0.82);
    leg->SetBorderSize(0); leg->SetTextSize(0.033);
    leg->AddEntry(ePre[k], "Pre-MVA (after distance cuts)", "lp");
    leg->AddEntry(eFin[k], "Final selection (BDTg > 0.96)",  "lp");
    leg->Draw();

    TString fname = TString("figs/EfficiencyPlots/") + fnames[k] + ".png";
    c->SaveAs(fname);
    delete c;
  }

  // ----------------------------------------------------------
  // Combined 2x3 canvas
  // ----------------------------------------------------------
  TCanvas *cAll = new TCanvas("Eff_All", "Selection Efficiency vs True Momentum", 1400, 900);
  cAll->Divide(3, 2, 0.003, 0.003);
  for (int k = 0; k < 6; k++) {
    TVirtualPad *pad = cAll->cd(k+1);
    pad->SetGrid();
    pad->SetLeftMargin(0.14); pad->SetBottomMargin(0.14);
    ePre[k]->Draw("AP");
    cAll->Update();
    TGraphAsymmErrors *gP = ePre[k]->GetPaintedGraph();
    gP->GetYaxis()->SetRangeUser(0.0, 0.50);
    gP->GetXaxis()->SetLimits(0.0, xmax[k]);
    gP->GetXaxis()->SetTitleSize(0.055);
    gP->GetYaxis()->SetTitleSize(0.055);
    gP->GetXaxis()->SetTitleOffset(1.05);
    gP->GetYaxis()->SetTitleOffset(1.05);
    pad->Modified(); cAll->Update();
    eFin[k]->Draw("P SAME");
    if (k == 0) {
      TLegend *leg = new TLegend(0.33, 0.68, 0.90, 0.82);
      leg->SetBorderSize(0); leg->SetTextSize(0.042);
      leg->AddEntry(ePre[k], "Pre-MVA", "lp");
      leg->AddEntry(eFin[k], "BDTg>0.96", "lp");
      leg->Draw();
    }
  }
  cAll->SaveAs("figs/EfficiencyPlots/Eff_All.png");
  delete cAll;

  inFile->Close();
  delete mvaReader;
  cout << "Saved: figs/EfficiencyPlots/" << endl;
}
