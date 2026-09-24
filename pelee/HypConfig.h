#ifndef HypConfig_h
#define HypConfig_h
// ============================================================================
//  HypConfig -- operating point for the hypana Sigma0 chain.
//
//  BDTG_CUT is the out-of-sample figure-of-merit peak, located by a 0.002-step
//  scan (PeakBDT.C). Rescanned 2026-09-11 after the proton-momentum fix
//  (HyperonProduction trk_energy_proton is momentum, not KE; see
//  hypana-proton-momentum-bug in the memory notes and DumpBDTHyp.C).
//
//                          S       B      eff[%]  pur[%]   FoM
//    preselection        58.68  20198.3   100.0    0.290  0.4123
//    old point +0.490    22.83    969.5    38.9    2.301  0.7248
//    FoM peak  +0.554    19.32    670.5    32.9    2.801  0.7357   <-- this
//
//  The peak is a plateau: the two points above differ by 1.5% in FoM. Moving
//  the cut buys purity (2.30 -> 2.80%) at a cost in signal (22.8 -> 19.3).
//
//  NOTE this buys purity, not signal. At 19 events over 670 background a 2%
//  systematic on the background prediction is +-13 events, comparable to the
//  signal itself, so the systematic budget -- not S/sqrt(S+B) -- decides whether
//  this point or a looser one is actually better. Rescan after any change to the
//  tree or the training.
// ============================================================================
namespace HypCfg {
  static const double BDTG_CUT = 0.554;   // FoM peak, 2026-09-11 (post proton fix)
}
#endif
