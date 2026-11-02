# Closed-form contraction (Gordon form A gamma^mu - B P^mu) checked against the numerical trace for every coupling combination.
import sys, numpy as np
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../../../pelee"))
import HypXsecWeights as H
from terms import trace_coeffs, pais_coeffs
ml2 = H.MMU**2; M = H.MN
def exact_bracket(pdg, x, Q2, f1, f2, g1, lam):
    M1, M2 = H.NUC_MASS[pdg], H.MASS[pdg]
    Sig = M1**2+M2**2+ml2+Q2; s, u = (Sig+x)/2, (Sig-x)/2
    k1p1, k2p2, k1p2, k2p1 = (s-M1**2)/2, (s-ml2-M2**2)/2, (M2**2-u)/2, (M1**2+ml2-u)/2
    k1k2, p1p2 = (ml2+Q2)/2, (M1**2+M2**2+Q2)/2
    p1p1, p2p2 = M1**2, M2**2
    def S(ka, kb, k2a, k2b, ab): return ka*k2b + k2a*kb - k1k2*ab      # S(a,b) with (k1.a,k1.b,k2.a,k2.b,a.b)
    k1P, k2P = k1p1+k1p2, k2p1+k2p2
    S_p2P = S(k1p2, k1P, k2p2, k2P, p1p2+p2p2); S_p1P = S(k1p1, k1P, k2p1, k2P, p1p1+p1p2)
    S_PP = S(k1P, k1P, k2P, k2P, 2*M1**2+2*M2**2+Q2)
    A, B, G = f1+f2*(M1+M2)/(2*M), f2/(2*M), g1
    sym = 2*k1p1*k2p2 + 2*k2p1*k1p2
    val = A*A*(sym-2*M1*M2*k1k2) - 2*A*B*(M1*S_p2P+M2*S_p1P) + B*B*(p1p2+M1*M2)*S_PP + G*G*(sym+2*M1*M2*k1k2) \
          + lam*A*G*(k1p1*k2p2-k1p2*k2p1)
    return 4*val
def exact_coeffs(pdg, Q2, coup, lam):
    xs = np.array([1., 3., 5., 8.]); ys = [exact_bracket(pdg, x, Q2, *coup, lam) for x in xs]
    c2, c1, c0 = np.polyfit(xs, ys, 2); return c0, c1, c2
# fix the interference normalisation lam from f1+g1 at one point
t = trace_coeffs(3212, 0.3, (1, 0, 1)); e0 = exact_coeffs(3212, 0.3, (1, 0, 1), 0.); e1 = exact_coeffs(3212, 0.3, (1, 0, 1), 1.)
lam = (t[1]-e0[1])/(e1[1]-e0[1]); print("lambda =", round(lam, 6))
cases = {"f1": (1, 0, 0), "g1": (0, 0, 1), "f2": (0, 1, 0), "f1+g1": (1, 0, 1), "f1+f2": (1, 1, 0), "f2+g1": (0, 1, 1), "all": (1, 1.5, 1.2)}
worst = 0
for pdg in (3212, 3122):
    for Q2 in (0.1, 0.5, 1.2):
        for lab, cp in cases.items():
            t = np.array(trace_coeffs(pdg, Q2, cp)); e = np.array(exact_coeffs(pdg, Q2, cp, lam)); g = np.array(pais_coeffs(pdg, Q2, *cp))
            d = np.max(np.abs(t-e)); worst = max(worst, d)
            if Q2 == 0.5 and pdg == 3212: print(f"  {lab:6s} trace {np.round(t,4)}  exact {np.round(e,4)}  GENIE {np.round(g,4)}")
print("max |trace - exact| over all cases:", worst)
