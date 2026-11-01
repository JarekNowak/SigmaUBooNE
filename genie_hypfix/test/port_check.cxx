// Port check: HypQELTracePXSec::Contraction against HypXsecWeights.py (same form factors:
// dipole, G_En = 0, M_A 0.99, g_A 1.267, F 0.463, D 0.804), free nucleon, 1e-40 cm2.
// Expected (python): 1 GeV  Lambda 1.33 Sigma0 0.25 Sigma- 0.50 ; 2 GeV  2.15 0.41 0.81
#include "HypQELTracePXSec.h"
#include <cmath>
#include <cstdio>
void ffs(int pdg, double Q2, double& f1, double& f2, double& fa) {
  const double MN = 0.93827, MV2 = 0.71, MA = 0.99, GA = 1.267, KP = 1.793, KN = -1.913, X = 0.463 / (0.463 + 0.804);
  double tau = Q2 / (4 * MN * MN), GD = 1 / pow(1 + Q2 / MV2, 2);
  double GEp = GD, GMp = (1 + KP) * GD, GEn = 0, GMn = KN * GD;
  double F1p = (GEp + tau * GMp) / (1 + tau), F2p = (GMp - GEp) / (1 + tau), F1n = (GEn + tau * GMn) / (1 + tau), F2n = (GMn - GEn) / (1 + tau);
  double FAd = -GA / pow(1 + Q2 / (MA * MA), 2);      // GENIE convention FA0 = -1.267
  if (pdg == 3122) { f1 = -sqrt(1.5) * F1p; f2 = -sqrt(1.5) * F2p; fa = -(1 + 2 * X) / sqrt(6.) * FAd; }
  else if (pdg == 3212) { f1 = -(F1p + 2 * F1n) / sqrt(2.); f2 = -(F2p + 2 * F2n) / sqrt(2.); fa = (1 - 2 * X) / sqrt(2.) * FAd; }
  else { f1 = -(F1p + 2 * F1n); f2 = -(F2p + 2 * F2n); fa = (1 - 2 * X) * FAd; }
}
int main() {
  const double GF = 1.1663787e-5, s2 = 1 - 0.9742 * 0.9742, hc2 = 0.3894e-27, ml = 0.10566;
  for (double E : {1.0, 2.0}) {
    printf("E = %.1f GeV:", E);
    for (int pdg : {3122, 3212, 3112}) {
      const double M = pdg == 3112 ? 0.93957 : 0.93827, MY = pdg == 3122 ? 1.115683 : pdg == 3212 ? 1.192642 : 1.197449;
      // Q2 range for a free nucleon at rest
      double s = M * M + 2 * M * E, rs = sqrt(s), El = (s + ml * ml - MY * MY) / (2 * rs), pl = sqrt(El * El - ml * ml), Ein = (s - M * M) / (2 * rs);
      double lo = -ml * ml + 2 * Ein * (El - pl), hi = -ml * ml + 2 * Ein * (El + pl);
      const int n = 2000; double sum = 0;
      for (int i = 0; i < n; i++) {
        double Q2 = lo + (i + 0.5) * (hi - lo) / n, f1, f2, fa; ffs(pdg, Q2, f1, f2, fa);
        double LH = genie::HypQELTracePXSec::Contraction(E, -Q2, M, MY, ml, f1, f2, fa);
        if (LH > 0) sum += GF * GF * s2 * LH / (256 * M_PI * M * M * E * E) * (hi - lo) / n;
      }
      printf("  %d %.3f", pdg, sum * hc2 / 1e-40);
    }
    printf("\n");
  }
  return 0;
}
