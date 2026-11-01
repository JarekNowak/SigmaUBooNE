//____________________________________________________________________________
/*!

\class    genie::HypQELTracePXSec

\brief    Delta S = 1 quasi-elastic hyperon production, nubar N -> l+ Y
          (Y = Lambda, Sigma0, Sigma-), from the full contraction of the
          leptonic and hadronic tensors with the hyperon and nucleon masses
          kept exactly.

          Drop-in replacement for genie::PaisQELLambdaPXSec: same channels,
          same form factors (FormFactorsAlg, LwlynSmithFFDeltaS), same Cabibbo
          angle, same integrator, same dsigma/dQ2 convention. PaisQELLambdaPXSec
          reproduces Llewellyn Smith CCQE when the final baryon has the nucleon
          mass but gets the hyperon-nucleon mass-difference terms wrong: against
          this contraction its Lambda is x0.80 and its Sigma0 and Sigma- x1.85,
          Sigma0/Lambda 0.45 instead of 0.19 (NuWro's Singh model gives 0.21 on
          Ar). See Sigma/README.md "GENIE hyperon cross section" and
          Sigma/pelee/HypXsecWeights.py, of which XSec() is a C++ port.

          dsigma/dQ2 = G_F^2 sin^2(theta_C) L_{mu nu} H^{mu nu} / (256 pi M^2 E^2),
          E the antineutrino energy in the struck-nucleon rest frame,
          L^{mu nu} from vbar(k1) g^mu (1-g5) v(k2) (antineutrino; the neutrino
          channels are not valid processes, as for Pais), H_{mu nu} the
          spin-summed hadronic tensor of
              F1V g^mu + i xiF2V/(2M) sigma^{mu nu} q_nu + FA g^mu g5,  q = k1 - k2.
          Checked against Llewellyn Smith CCQE (nubar 0.41e-38 cm2 at 1 GeV).

\author   J. Nowak (with Claude), 2026-10-31
*/
//____________________________________________________________________________

#ifndef _HYP_QEL_TRACE_PARTIAL_XSEC_H_
#define _HYP_QEL_TRACE_PARTIAL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"

namespace genie {

class QELFormFactorsModelI;
class XSecIntegratorI;

class HypQELTracePXSec : public XSecAlgorithmI {

public:
  HypQELTracePXSec();
  HypQELTracePXSec(string config);
  virtual ~HypQELTracePXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  void Configure (const Registry & config);
  void Configure (string param_set);

  //-- the contraction itself, exposed for validation: L_{mu nu} H^{mu nu} for
  //   nubar(E) + N(M, at rest) -> l+(ml) + Y(MY) at q2 (< 0), given the form factors.
  //   Returns -1 outside the physical region.
  static double Contraction(double E, double q2, double M, double MY, double ml,
                            double F1V, double xiF2V, double FA);

private:
  void   LoadConfig (void);
  double MHyperon   (const Interaction * interaction) const;

  mutable QELFormFactors          fFormFactors;
  const   QELFormFactorsModelI *  fFormFactorsModel;
  const   XSecIntegratorI *       fXSecIntegrator;
  double                          fSin8c2;
};

} // genie namespace
#endif  // _HYP_QEL_TRACE_PARTIAL_XSEC_H_
