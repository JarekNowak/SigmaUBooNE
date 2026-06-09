#ifndef SigmaConfig_h
#define SigmaConfig_h
// ============================================================================
//  Shared configuration for the Sigma analysis: sample/cut/variable names,
//  per-variable histogram binning, and plot colors.
//
//  Included by BOTH Sigma.C (histogram booking) and PlotSigma.C (plotting),
//  so the two macros can never drift out of sync.
// ============================================================================
#include "TString.h"
#include "Rtypes.h"   // EColor (kRed, kAzure, ...)

namespace SigmaCfg {

  const int Ns = 5;   // samples
  const int Nc = 14;  // cumulative cut stages
  const int Nv = 83;  // variables

  // --- Sample labels (index order matches Sigma.C: 0=signal .. 4=data) ---
  static const TString Names[Ns] = {
    "Sigma", "Hyperons", "Background", "Dirt", "Data"
  };

  // --- Cut stage labels ---
  static const TString CutName[Nc] = {
    "AllEvents", "FiducialV", "NShowers", "NoBigShowers", "MuonCandidate",
    "PionCandidate", "ProtonCandidate", "NTracksCut", "2ndShower", "1stShower",
    "Distances", "BDT", "BDTAda", "BDTg"
  };

  // --- Variable labels (also used as PNG file names) ---
  static const TString VarName[Nv] = {
    /* 0  */ "NTracks", "NShowers", "NHitsY",
    /* 3  */ "MuonPID", "MuonTrkLength", "MuonTrkDistVtx", "MuonTrkScore",
    /* 7  */ "PionPID", "PionTrkLength", "PionTrkDistVtx", "PionTrkScore",
    /* 11 */ "ProtonPID", "ProtonTrkLength", "ProtonTrkDistVtx", "ProtonTrkScore",
    /* 15 */ "GammaPID", "GammaTrkLength", "GammaTrkDistVtx", "GammaTrkScore",
    /* 19 */ "MuonNHitsY", "PionNHitsY", "ProtonNHitsY", "GammaNHitsY",
    /* 23 */ "MuonNHitsU", "PionNHitsU", "ProtonNHitsU", "GammaNHitsU",
    /* 27 */ "MuonNHitsV", "PionNHitsV", "ProtonNHitsV", "GammaNHitsV",
    /* 31 */ "MuonMeanDeDxPlane0", "MuonMeanDeDxPlane1", "MuonMeanDeDxPlane2",
    /* 34 */ "PionMeanDeDxPlane0", "PionMeanDeDxPlane1", "PionMeanDeDxPlane2",
    /* 37 */ "ProtonMeanDeDxPlane0", "ProtonMeanDeDxPlane1", "ProtonMeanDeDxPlane2",
    /* 40 */ "GammaMeanDeDxPlane0", "GammaMeanDeDxPlane1", "GammaMeanDeDxPlane2",
    /* 43 */ "MuonPhotonPID", "PionPhotonPID", "ProtonPhotonPID", "GammaPhotonPID",
    /* 47 */ "MuonCandShrOpenAngle", "PionCandShrOpenAngle", "ProtonCandShrOpenAngle",
    /* 50 */ "Shower1CandShrOpenAngle", "Shower2CandShrOpenAngle",
    /* 52 */ "Shower1CandDeDxPlane0", "Shower1CandDeDxPlane1", "Shower1CandDeDxPlane2",
    /* 55 */ "Shower2CandDeDxPlane0", "Shower2CandDeDxPlane1", "Shower2CandDeDxPlane2",
    /* 58 */ "OpeningAngleMuonPion", "OpeningAngleProtonPion", "OpeningAngleProtonMuon",
    /* 61 */ "OpeningAngleShower1Muon", "OpeningAngleShower1Pion", "OpeningAngleShower1Proton",
    /* 64 */ "OpeningAngleShower2Muon", "OpeningAngleShower2Pion", "OpeningAngleShower2Proton",
    /* 67 */ "Shower1NHitsY", "Shower2NHitsY",
    /* 69 */ "Shower1VtxDist", "Shower2VtxDist",
    /* 71 */ "DistMuonPion", "DistProtonPion", "DistProtonMuon",
    /* 74 */ "DistShower1Muon", "DistShower1Pion", "DistShower1Proton",
    /* 77 */ "DistShower2Muon", "DistShower2Pion", "DistShower2Proton",
    /* 80 */ "BDTScore", "BDTAdaScore", "BDTgScore"
  };

  // --- Per-sample plot colors ---
  //   Sigma  : red line (signal, scaled up in stacked plots)
  //   Hyper. : green fill, Background : blue fill, Dirt : gray fill
  //   Data   : black markers
  static const int SampleColor[Ns] = {
    kRed + 1, kGreen + 2, kAzure + 1, kGray + 1, kBlack
  };

  // --- Per-variable binning, grouped by variable type. ---
  //   Ranges are sensible starting defaults; tune here and both macros follow.
  //   Returns nbins, [lo,hi), and whether the plot should use log-y.
  inline void Binning(int v, int &nb, double &lo, double &hi, bool &logy)
  {
    nb = 50; lo = 0.; hi = 1.; logy = false;  // default: scores, BDT, cand open angles

    if      (v == 0)                         { nb = 10; lo = 0;  hi = 10;     } // NTracks
    else if (v == 1)                         { nb = 8;  lo = 0;  hi = 8;      } // NShowers
    else if (v == 2)                         { nb = 50; lo = 0;  hi = 1000; logy = true; } // NHitsY
    else if (v >= 19 && v <= 30)             { nb = 50; lo = 0;  hi = 1000; logy = true; } // *NHits[YUV]
    else if (v == 67 || v == 68)             { nb = 50; lo = 0;  hi = 500;  logy = true; } // shower NHitsY
    else if (v == 3 || v == 7 || v == 11 || v == 15) { nb = 50; lo = -1; hi = 1;   } // *PID (llrpid)
    else if (v == 4 || v == 8 || v == 12 || v == 16) { nb = 50; lo = 0;  hi = 200; } // *TrkLength
    else if (v == 5 || v == 9 || v == 13 || v == 17) { nb = 50; lo = 0;  hi = 50;  } // *TrkDistVtx
    else if (v == 69 || v == 70)             { nb = 50; lo = 0;  hi = 100;   } // shower vtx dist
    else if ((v >= 31 && v <= 42) || (v >= 52 && v <= 57)) { nb = 50; lo = 0; hi = 10; } // dE/dx
    else if (v >= 58 && v <= 66)             { nb = 50; lo = 0;  hi = 3.1416;} // opening angles (rad)
    else if (v >= 71 && v <= 79)             { nb = 50; lo = 0;  hi = 200;   } // candidate distances
    else if (v >= 80 && v <= 82)             { nb = 50; lo = -1; hi = 1;     } // BDT scores
    // else: default 50,0,1 -> *TrkScore(6,10,14,18), *PhotonPID(43-46),
    //       CandShrOpenAngle(47-51)
  }

} // namespace SigmaCfg

#endif // SigmaConfig_h
