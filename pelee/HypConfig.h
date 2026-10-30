#ifndef HypConfig_h
#define HypConfig_h
// ============================================================================
//  HypConfig -- operating point for the hypana Sigma0 chain.
//
//  BDTG_CUT sits on the out-of-sample figure-of-merit plateau, located by a
//  0.002-step scan (PeakBDT.C). Rescanned 2026-10-30 after retraining with the GENIE
//  Sigma0 QE sample as signal (47,033 distinct events; the previous hyperon sample
//  was 754 interactions overlaid ~20 times), 95% of each fold, full background:
//
//                          S       B      eff[%]  pur[%]   FoM
//    preselection        44.07  20112.6   100.0    0.219  0.3104
//    FoM peak  +0.608    14.30    504.0    32.4    2.759  0.6280
//    kept      +0.554    17.13    733.2    38.9    2.283  0.6253   <-- this
//
//  The plateau runs from +0.46 to +0.66 within 1.5% in FoM, so the cut was kept,
//  0.4% below the peak. On 2026-10-19 (old signal sample) the same cut gave
//  S 20.19, B 711.6, FoM 0.7463; S is lower now because the QE sample's
//  normalisation is 0.74x the old sample's, and the AUC rose 0.834 -> 0.855.

//  NOTE this buys purity, not signal. At 17 events over 733 background a 2%
//  systematic on the background prediction is +-14 events, comparable to the
//  signal itself, so the systematic budget -- not S/sqrt(S+B) -- decides whether
//  this point or a looser one is actually better (Sensitivity.C does this for the
//  PeLEE chain: the best cut tightens as the systematic grows). Rescan after any
//  change to the tree or the training.
// ============================================================================
namespace HypCfg {
  static const double BDTG_CUT = 0.554;   // FoM plateau, rescanned 2026-10-30
}
#endif
