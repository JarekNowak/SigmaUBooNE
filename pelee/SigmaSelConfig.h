#ifndef SigmaSelConfig_h
#define SigmaSelConfig_h
// ============================================================================
//  Shared configuration for the PeLEE-based Sigma0 selection: sample / cut /
//  variable names, per-variable histogram binning, and plot colors.
//
//  Included by BOTH SigmaSelection.C (histogram booking) and PlotSigma.C
//  (plotting), so the two macros can never drift out of sync.
//
//  This replaces the archived SigmaConfig.h.  The variable list has been
//  trimmed to quantities that exist cleanly in standard PeLEE ntuples
//  (tree nuselection/NeutrinoSelectionFilter); the old per-plane track dE/dx
//  and photon-BDT ("GammaPID"/"PhotonPID") slots had no PeLEE equivalent and
//  were dropped.  The four old pass-through cut slots (NoBigShowers, 2ndShower,
//  BDT, BDTAda) are also gone — the BDTg cut will be appended in the BDT pass.
// ============================================================================
#include "TString.h"
#include "Rtypes.h"   // EColor (kRed, kAzure, ...)

namespace SigmaCfg {

  const int Ns = 6;   // samples
  const int Nc = 9;   // cumulative cut stages
  const int Nv = 53;  // variables

  // --- Sample labels (index order matches SigmaSelection.C: 0=signal .. 4=data,
  //     5=EXT beam-off, appended so that no earlier index moves) ---
  static const TString Names[Ns] = {
    "Sigma", "Hyperons", "Background", "Dirt", "Data", "EXT"
  };

  // --- Cut stage labels (real cuts only; cumulative) ---
  //  Final discriminant is a gradient BDT (BDTG) applied at the NTracksCut
  //  preselection.  It subsumes the old hard discriminating cuts (proton-muon
  //  distance, photon shower-size, Λ mass window) — all of which are now BDT
  //  inputs (DistProtonMuon, Shower1NHitsY, LambdaMass, plus Σ⁰ mass, shower
  //  energy, Λ/γ vertex-pointing, ...), so their power is retained and combined
  //  optimally.  Trained by TrainBDT.C on BDTInput.root (dumped here).
  static const TString CutName[Nc] = {
    "AllEvents", "FiducialV", "CosmicVeto", "NShowers", "MuonCandidate",
    "ProtonCandidate", "PionCandidate", "NTracksCut", "BDTG"
  };

  // --- Variable labels (also used as PNG file names) ---
  static const TString VarName[Nv] = {
    /* 0  */ "NTracks", "NShowers", "TopologicalScore",
    /* 3  */ "MuonPID", "MuonTrkLength", "MuonTrkDistVtx", "MuonTrkScore", "MuonNHitsY",
    /* 8  */ "PionPID", "PionTrkLength", "PionTrkDistVtx", "PionTrkScore", "PionNHitsY", "PionCaloEnergyY",
    /* 14 */ "ProtonPID", "ProtonTrkLength", "ProtonTrkDistVtx", "ProtonTrkScore", "ProtonNHitsY", "ProtonCaloEnergyY",
    /* 20 */ "Shower1VtxDist", "Shower1NHitsY", "Shower1DeDxY", "Shower1MoliereAvg",
    /* 24 */ "Shower2VtxDist", "Shower2NHitsY", "Shower2DeDxY",
    /* 27 */ "OpeningAngleMuonPion", "OpeningAngleProtonPion", "OpeningAngleProtonMuon",
    /* 30 */ "OpeningAngleShower1Muon", "OpeningAngleShower1Pion", "OpeningAngleShower1Proton",
    /* 33 */ "DistMuonPion", "DistProtonPion", "DistProtonMuon",
    /* 36 */ "DistShower1Muon", "DistShower1Pion", "DistShower1Proton",
    /* 39 */ "LambdaMass", "ProtonChi2",
    /* 41 */ "SigmaMass",
    /* 42 */ "LambdaPointAngle", "LambdaImpact", "LambdaDecayDist",
    /* 45 */ "PhotonPointAngle", "PhotonImpact",
    /* 47 */ "LambdaGammaAngle", "SigmaMuonAngle",
    /* 49 */ "BDTGScore",
    /* 50 */ "Pi0Mass", "Pi0Gammadot", "Shower2Energy"
  };

  // --- Per-sample plot colors ---
  //   Sigma  : red line (signal, scaled up in stacked plots)
  //   Hyper. : green fill, Background : blue fill, Dirt : gray fill
  //   Data   : black markers, EXT : orange fill
  static const int SampleColor[Ns] = {
    kRed + 1, kGreen + 2, kAzure + 1, kGray + 1, kBlack, kOrange + 1
  };

  // --- Per-variable binning. Returns nbins, [lo,hi), and log-y flag. ---
  inline void Binning(int v, int &nb, double &lo, double &hi, bool &logy)
  {
    nb = 50; lo = 0.; hi = 1.; logy = false;  // default: track/shower scores

    if      (v == 0)                              { nb = 10; lo = 0;  hi = 10;   } // NTracks
    else if (v == 1)                              { nb = 8;  lo = 0;  hi = 8;    } // NShowers
    else if (v == 2)                              { nb = 50; lo = 0;  hi = 1;    } // TopologicalScore
    else if (v == 3 || v == 8 || v == 14)         { nb = 50; lo = -1; hi = 1;    } // *PID (llrpid)
    else if (v == 4 || v == 9 || v == 15)         { nb = 50; lo = 0;  hi = 200;  } // *TrkLength
    else if (v == 5 || v == 10 || v == 16)        { nb = 50; lo = 0;  hi = 50;   } // *TrkDistVtx
    else if (v == 6 || v == 11 || v == 17)        { nb = 50; lo = 0;  hi = 1;    } // *TrkScore
    else if (v == 7 || v == 12 || v == 18)        { nb = 50; lo = 0;  hi = 1000; logy = true; } // *NHitsY (track)
    else if (v == 13 || v == 19)                  { nb = 50; lo = 0;  hi = 1;    } // *CaloEnergyY (GeV)
    else if (v == 20 || v == 24)                  { nb = 50; lo = 0;  hi = 100;  } // shower vtx dist
    else if (v == 21 || v == 25)                  { nb = 50; lo = 0;  hi = 500;  logy = true; } // shower NHitsY
    else if (v == 22 || v == 26)                  { nb = 50; lo = 0;  hi = 10;   } // shower dE/dx
    else if (v == 23)                             { nb = 50; lo = 0;  hi = 0.5;  } // Moliere avg (rad)
    else if (v >= 27 && v <= 32)                  { nb = 50; lo = 0;  hi = 3.1416;} // opening angles (rad)
    else if (v >= 33 && v <= 38)                  { nb = 50; lo = 0;  hi = 200;  } // candidate distances
    else if (v == 39)                             { nb = 60; lo = 1.0; hi = 1.6;  } // Lambda (p pi-) invariant mass (GeV)
    else if (v == 40)                             { nb = 50; lo = 0;  hi = 400;  } // proton candidate chi2_proton
    else if (v == 41)                             { nb = 60; lo = 1.05; hi = 1.65;} // Sigma0 (Λ+γ) invariant mass (GeV)
    else if (v == 42 || v == 45 || v == 47 || v == 48) { nb = 50; lo = 0; hi = 3.1416; } // pointing/decay angles (rad)
    else if (v == 43 || v == 46)                  { nb = 50; lo = 0;  hi = 50;   } // Lambda/photon impact param (cm)
    else if (v == 44)                             { nb = 50; lo = 0;  hi = 50;   } // Lambda decay distance (cm)
    else if (v == 49)                             { nb = 80; lo = -1; hi = 1;    } // BDTG cross-val score
    else if (v == 50)                             { nb = 50; lo = 0;  hi = 0.5;  } // pi0 di-photon mass (GeV)
    else if (v == 51)                             { nb = 50; lo = -1; hi = 1;    } // pi0 gamma-gamma cos(angle)
    else if (v == 52)                             { nb = 50; lo = 0;  hi = 0.5;  } // shower2 energy (GeV)
  }

} // namespace SigmaCfg

#endif // SigmaSelConfig_h
