#ifndef HypConfig_h
#define HypConfig_h
// ============================================================================
//  HypConfig -- operating point for the hypana Sigma0 chain.
//
//  BDTG_CUT sits on the out-of-sample figure-of-merit plateau, located by a
//  0.002-step scan (PeakBDT.C). Rescanned 2026-10-19 after retraining on 95% of
//  each fold (was TMVA's half, i.e. 25% of all events) on the full 14,164-file
//  background (was Run 3 only, 2.54x less POT), with the 2x20 island tuning and
//  the Shower1Energy cap in DumpBDTHyp.C:
//
//                          S       B      eff[%]  pur[%]   FoM
//    preselection        58.68  20112.6   100.0    0.291  0.4132
//    FoM peak  +0.540    21.06    768.0    35.9    2.669  0.7496
//    kept      +0.554    20.19    711.6    34.4    2.759  0.7463   <-- this
//
//  The plateau runs from +0.39 to +0.59 within 1.5% in FoM, and the peak inside
//  it is noise: a tree differing in 8 of 123,000 background values put it at
//  +0.434. So the cut was kept where it was, 0.4% below the peak.
//  Before the retraining (Run 3 background, half-fold training, 2026-09-11) the
//  same cut gave S 19.32, B 670.5, FoM 0.7357.
//
//  NOTE this buys purity, not signal. At 20 events over 712 background a 2%
//  systematic on the background prediction is +-14 events, comparable to the
//  signal itself, so the systematic budget -- not S/sqrt(S+B) -- decides whether
//  this point or a looser one is actually better (Sensitivity.C does this for the
//  PeLEE chain: the best cut tightens as the systematic grows). Rescan after any
//  change to the tree or the training.
// ============================================================================
namespace HypCfg {
  static const double BDTG_CUT = 0.554;   // FoM plateau, rescanned 2026-10-19
}
#endif
