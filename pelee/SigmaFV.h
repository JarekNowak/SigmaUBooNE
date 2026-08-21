#ifndef SigmaFV_h
#define SigmaFV_h
// ============================================================================
//  Fiducial-volume helper for the PeLEE-based Sigma0 selection.
//
//  Box and dead-region taken verbatim from the archived Sigma.C (truth-FV
//  flags), so the new selection reproduces the old fiducial definition:
//    x in [3.00, 253.35], y in [-112.53, 114.47], z in [3.1, 1033.0]
//    excluding the dead region z in (675.1, 775.1).
//
//  Applied to the reco vertex (reco_nu_vtx_sce_*) for every sample, and to the
//  true neutrino vertex (true_nu_vtx_*) for the MC signal / hyperon samples.
// ============================================================================

namespace SigmaFV {

  const double Xmin =   3.00,  Xmax = 253.35;
  const double Ymin = -112.53, Ymax = 114.47;
  const double Zmin =   3.10,  Zmax = 1033.0;
  const double DeadZmin = 675.1, DeadZmax = 775.1;

  inline bool isFiducial(double x, double y, double z) {
    if (x < Xmin || x > Xmax) return false;
    if (y < Ymin || y > Ymax) return false;
    if (z < Zmin || z > Zmax) return false;
    if (z > DeadZmin && z < DeadZmax) return false;  // dead region
    return true;
  }

} // namespace SigmaFV

#endif // SigmaFV_h
