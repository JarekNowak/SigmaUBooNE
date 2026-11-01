//____________________________________________________________________________
/*
 genie::HypQELTracePXSec -- see the header for the physics and the validation.
 The process/kinematics/configuration methods are those of
 genie::PaisQELLambdaPXSec (H. Gallagher), unchanged; only XSec() differs.
*/
//____________________________________________________________________________

#include <complex>
#include <cmath>

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"

#include "HypQELTracePXSec.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

namespace {

typedef std::complex<double> cd;
struct M4 { cd a[4][4]; };

M4 zero() { M4 m; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) m.a[i][j] = 0.; return m; }
M4 mul(const M4& x, const M4& y) {
  M4 r = zero();
  for (int i = 0; i < 4; i++) for (int k = 0; k < 4; k++) {
    const cd v = x.a[i][k]; if (v == 0.) continue;
    for (int j = 0; j < 4; j++) r.a[i][j] += v * y.a[k][j];
  }
  return r;
}
M4 add(const M4& x, const M4& y, cd cy = 1.) {
  M4 r; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) r.a[i][j] = x.a[i][j] + cy * y.a[i][j]; return r;
}
M4 scale(const M4& x, cd c) {
  M4 r; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) r.a[i][j] = c * x.a[i][j]; return r;
}
cd trace(const M4& x) { return x.a[0][0] + x.a[1][1] + x.a[2][2] + x.a[3][3]; }
cd traceProd(const M4& x, const M4& y) {   // Tr(x y) without forming the product
  cd t = 0.; for (int i = 0; i < 4; i++) for (int k = 0; k < 4; k++) t += x.a[i][k] * y.a[k][i]; return t;
}

// Bjorken-Drell gamma^mu (upper index), gamma5 = i g0 g1 g2 g3, metric (+,-,-,-)
struct Dirac {
  M4 g[4], g5, one;
  M4 sig[4][4];            // sigma^{mu nu} = i/2 [g^mu, g^nu]
  Dirac() {
    for (int m = 0; m < 4; m++) g[m] = zero();
    one = zero(); for (int i = 0; i < 4; i++) one.a[i][i] = 1.;
    g[0].a[0][0] = 1; g[0].a[1][1] = 1; g[0].a[2][2] = -1; g[0].a[3][3] = -1;
    const cd I(0, 1);
    // g^k = [[0, s_k], [-s_k, 0]]
    cd s[3][2][2] = { { {0, 1}, {1, 0} }, { {0, -I}, {I, 0} }, { {1, 0}, {0, -1} } };
    for (int k = 0; k < 3; k++)
      for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
        g[k + 1].a[i][j + 2] = s[k][i][j];
        g[k + 1].a[i + 2][j] = -s[k][i][j];
      }
    g5 = scale(mul(mul(g[0], g[1]), mul(g[2], g[3])), I);
    for (int m = 0; m < 4; m++) for (int n = 0; n < 4; n++)
      sig[m][n] = scale(add(mul(g[m], g[n]), mul(g[n], g[m]), -1.), 0.5 * I);
  }
  M4 slash(const double p[4]) const {   // g^mu p_mu
    M4 r = scale(g[0], p[0]);
    for (int k = 1; k < 4; k++) r = add(r, g[k], -p[k]);
    return r;
  }
  M4 bar(const M4& x) const {           // g0 x^dagger g0
    M4 d; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) d.a[i][j] = std::conj(x.a[j][i]);
    return mul(mul(g[0], d), g[0]);
  }
};
const Dirac& D() { static const Dirac d; return d; }
const double MET[4] = { 1, -1, -1, -1 };

} // namespace

//____________________________________________________________________________
HypQELTracePXSec::HypQELTracePXSec() :
XSecAlgorithmI("genie::HypQELTracePXSec")
{
}
//____________________________________________________________________________
HypQELTracePXSec::HypQELTracePXSec(string config) :
XSecAlgorithmI("genie::HypQELTracePXSec", config)
{
}
//____________________________________________________________________________
HypQELTracePXSec::~HypQELTracePXSec()
{
}
//____________________________________________________________________________
double HypQELTracePXSec::Contraction(double E, double q2, double M, double MY, double ml,
                                     double F1V, double xiF2V, double FA)
{
  // four-vectors in the struck-nucleon rest frame, nubar along z
  const double q0 = (MY * MY - M * M - q2) / (2 * M);      // energy transfer
  const double El = E - q0;
  const double pl2 = El * El - ml * ml;
  if (El <= ml || pl2 <= 0) return -1;
  const double pl = std::sqrt(pl2);
  const double c = (q2 - ml * ml + 2 * E * El) / (2 * E * pl);
  if (c > 1 + 1e-9 || c < -1 - 1e-9) return -1;
  const double cc = std::max(-1.0, std::min(1.0, c)), sn = std::sqrt(std::max(0.0, 1 - cc * cc));
  const double k1[4] = { E, 0, 0, E };
  const double p1[4] = { M, 0, 0, 0 };
  const double k2[4] = { El, pl * sn, 0, pl * cc };
  double p2[4], ql[4];
  for (int m = 0; m < 4; m++) { p2[m] = k1[m] + p1[m] - k2[m]; ql[m] = MET[m] * (k1[m] - k2[m]); }

  const Dirac& d = D();
  const cd I(0, 1);
  // hadronic vertex Gamma^mu = F1V g^mu + i xiF2V/(2M) sigma^{mu nu} q_nu + FA g^mu g5
  M4 Gam[4], Gbar[4];
  for (int m = 0; m < 4; m++) {
    M4 sq = zero();
    for (int n = 0; n < 4; n++) sq = add(sq, d.sig[m][n], ql[n]);
    Gam[m] = add(add(scale(d.g[m], F1V), sq, I * xiF2V / (2 * M)), mul(d.g[m], d.g5), FA);
    Gbar[m] = d.bar(Gam[m]);
  }
  const M4 P1 = add(d.slash(p1), d.one, M);
  const M4 P2 = add(d.slash(p2), d.one, MY);
  const M4 K1 = d.slash(k1);
  const M4 K2 = add(d.slash(k2), d.one, -ml);
  const M4 PL = add(d.one, d.g5, -1.);                       // 1 - g5

  // L^{mn} = Tr[g^m (1-g5) (k2slash - ml) g^n (1-g5) k1slash]
  // H^{mn} = Tr[(p2slash + MY) Gam^m (p1slash + M) Gbar^n]
  M4 A[4], B[4], C[4], E4[4];
  for (int m = 0; m < 4; m++) {
    A[m] = mul(mul(d.g[m], PL), K2);          // g^m (1-g5) K2
    B[m] = mul(mul(d.g[m], PL), K1);          // g^n (1-g5) K1
    C[m] = mul(mul(P2, Gam[m]), P1);          // P2 Gam^m P1
    E4[m] = Gbar[m];
  }
  cd tot = 0.;
  for (int m = 0; m < 4; m++)
    for (int n = 0; n < 4; n++) {
      const cd L = traceProd(A[m], B[n]);
      const cd H = traceProd(C[m], E4[n]);
      tot += MET[m] * MET[n] * L * H;
    }
  return tot.real();
}
//____________________________________________________________________________
double HypQELTracePXSec::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics  = interaction->Kine();
  const InitialState & init_state  = interaction->InitState();
  const Target &       target      = init_state.Tgt();

  const double E    = init_state.ProbeE(kRfHitNucRest);
  const double q2   = kinematics.q2();
  const double Mnuc = target.HitNucMass();
  const double MY   = this->MHyperon(interaction);
  const double ml   = interaction->FSPrimLepton()->Mass();

  // Antineutrinos only (ValidProcess admits only the nubar channels' final states;
  // a neutrino would need the hermitian-conjugate current).
  if (pdg::IsNeutrino(init_state.ProbePdg())) return 0.;

  fFormFactors.Calculate(interaction);
  const double LH = Contraction(E, q2, Mnuc, MY, ml,
                                fFormFactors.F1V(), fFormFactors.xiF2V(), fFormFactors.FA());
  if (LH <= 0) return 0.;

  double xsec = kGF2 * fSin8c2 * LH / (256. * kPi * Mnuc * Mnuc * E * E);

  //----- The algorithm computes dxsec/dQ^2
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
    xsec *= J;
  }

  //----- free nucleon or nuclear target (simple scaling, as Pais)
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;
  int nuc   = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nuc)) ? target.Z() : target.N();
  xsec *= NNucl;

  return xsec;
}
//____________________________________________________________________________
double HypQELTracePXSec::MHyperon(const Interaction * interaction) const
{
  const XclsTag & xcls = interaction->ExclTag();
  int pdgc  = xcls.StrangeHadronPdg();
  return PDGLibrary::Instance()->Find(pdgc)->Mass();
}
//____________________________________________________________________________
double HypQELTracePXSec::Integral(const Interaction * interaction) const
{
  return fXSecIntegrator->Integrate(this,interaction);
}
//____________________________________________________________________________
bool HypQELTracePXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const XclsTag &      xcls       = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  bool is_exclusive_strange = (xcls.IsStrangeEvent() && !xcls.IsInclusiveStrange());
  if(!is_exclusive_strange) return false;

  if(!proc_info.IsQuasiElastic()) return false;
  if(!proc_info.IsWeak())         return false;

  bool isP = pdg::IsProton ( init_state.Tgt().HitNucPdg() );
  bool isN = pdg::IsNeutron( init_state.Tgt().HitNucPdg() );

  int pdgc = xcls.StrangeHadronPdg();

  return ( (pdgc == kPdgSigmaM && isN) ||
           (pdgc == kPdgLambda && isP) ||
           (pdgc == kPdgSigma0 && isP) );
}
//____________________________________________________________________________
bool HypQELTracePXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state  = interaction->InitState();
  double E = init_state.ProbeE(kRfHitNucRest);

  double MR    = this -> MHyperon  (interaction);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Mnuc  = init_state.Tgt().HitNucP4Ptr()->M();
  double Mnuc2 = TMath::Power(Mnuc,2);

  double ER = ( TMath::Power(MR+ml,2) - Mnuc2 ) / (2*Mnuc);
  return E > ER;
}
//____________________________________________________________________________
void HypQELTracePXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HypQELTracePXSec::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void HypQELTracePXSec::LoadConfig(void)
{
  double thc ;
  GetParam( "CabibboAngle", thc ) ;
  fSin8c2 = TMath::Power(TMath::Sin(thc), 2);

  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
                                             this->SubAlg("FormFactorsAlg"));
  assert(fFormFactorsModel);
  fFormFactors.SetModel(fFormFactorsModel);

  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
