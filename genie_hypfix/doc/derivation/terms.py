# Coefficients a + b (s-u) + c (s-u)^2 of the dsigma/dQ2 bracket at fixed Q2: exact trace vs GENIE PaisQELLambdaPXSec, per coupling.
import sys, numpy as np
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../../../pelee"))
import HypXsecWeights as H
MMU = H.MMU
def lims(pdg, E):
    Mn, MY = H.NUC_MASS[pdg], H.MASS[pdg]; s = Mn*Mn+2*Mn*E; rs = np.sqrt(s)
    El = (s+MMU**2-MY**2)/(2*rs); pl = np.sqrt(El*El-MMU**2); Ein = (s-Mn*Mn)/(2*rs)
    return -MMU**2+2*Ein*(El-pl), -MMU**2+2*Ein*(El+pl)
def bracket_trace(pdg, E, Q2, coup):
    lo, hi = lims(pdg, E); u = (Q2-lo)/(hi-lo)
    k = H.kinematics(pdg, E, u); s, Ein, pl = k[5], k[6], k[7]
    tr = H.msq_trace(pdg, E, u, coup=coup)*pl/(32*np.pi*s*Ein)/4/(2*Ein*pl)   # dsigma/dQ2 / (GF^2 sin^2)
    M1 = H.NUC_MASS[pdg]
    return tr*2*np.pi*16*M1*M1*E*E          # in Pais bracket units
def pais_coeffs(pdg, Q2, f1, f2, g1):
    M1, M2 = H.NUC_MASS[pdg], H.MASS[pdg]; q2 = -Q2; Mn2 = M1*M1
    F1V, xiF2V, FA = f1, f2, -g1
    v = (M2**2-Mn2-q2)/(2*M1); Mp, Mm = M2+M1, M2-M1
    w1 = (Mm**2-q2)/(4*Mn2)*(F1V+xiF2V)**2+(Mp**2-q2)/(4*Mn2)*FA**2
    w2 = FA**2+(F1V+xiF2V-Mp*xiF2V/(2*M1))**2-q2/Mn2*(xiF2V/2)**2
    w3 = 2*FA*(F1V+xiF2V)
    return -8*Mn2*q2*w1-4*(Mn2*v*v-q2)*w2, -2*q2*w3, w2      # a, b (antineutrino), c
def trace_coeffs(pdg, Q2, coup):
    M1 = H.NUC_MASS[pdg]
    Es = np.array([2.5, 3.5, 5.0, 7.0, 10.0]); x = []; y = []
    for E in Es:
        lo, hi = lims(pdg, E)
        if not (lo < Q2 < hi): continue
        M2 = H.MASS[pdg]; v = (M2**2-M1*M1+Q2)/(2*M1)
        s = M1*M1+2*M1*E; uu = M1*M1+MMU**2+2*v*M1-2*M1*E
        x.append(s-uu); y.append(bracket_trace(pdg, E, Q2, coup))
    c2, c1, c0 = np.polyfit(x, y, 2)
    return c0, c1, c2
cases = {"f1 only": (1, 0, 0), "g1 only": (0, 0, 1), "f2 only": (0, 1, 0), "f1+g1": (1, 0, 1), "f1+f2": (1, 1, 0)}
if __name__ == "__main__":
  for pdg, nm in ((3212, "Sigma0"), (3122, "Lambda")):
      for Q2 in (0.1, 0.5):
          print(f"== {nm}  Q2={Q2}:  coefficient   a (const)        b (s-u)          c (s-u)^2   [trace | GENIE]")
          for lab, (f1, f2, g1) in cases.items():
              t = trace_coeffs(pdg, Q2, (f1, f2, g1)); p = pais_coeffs(pdg, Q2, f1, f2, g1)
              print(f"   {lab:8s}  " + "  ".join(f"{a:8.4f} | {b:8.4f}" for a, b in zip(t, p)))
