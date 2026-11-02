# Symbolic derivation: exact coefficients, GENIE expression with each fix, and the effect of each fix on sigma(E).
import sympy as sp
M1, M2, M, ml2, Q2, x, A, B, G = sp.symbols('M1 M2 M ml2 Q2 x A B G', real=True)
Sig = M1**2+M2**2+ml2+Q2; s, u = (Sig+x)/2, (Sig-x)/2
k1p1, k2p2, k1p2, k2p1 = (s-M1**2)/2, (s-ml2-M2**2)/2, (M2**2-u)/2, (M1**2+ml2-u)/2
k1k2, p1p2 = (ml2+Q2)/2, (M1**2+M2**2+Q2)/2
def S(ka, kb, k2a, k2b, ab): return ka*k2b + k2a*kb - k1k2*ab
k1P, k2P = k1p1+k1p2, k2p1+k2p2
S_p2P = S(k1p2, k1P, k2p2, k2P, p1p2+M2**2); S_p1P = S(k1p1, k1P, k2p1, k2P, M1**2+p1p2)
S_PP = S(k1P, k1P, k2P, k2P, 2*M1**2+2*M2**2+Q2)
sym = 2*k1p1*k2p2 + 2*k2p1*k1p2
val = 4*(A**2*(sym-2*M1*M2*k1k2) - 2*A*B*(M1*S_p2P+M2*S_p1P) + B**2*(p1p2+M1*M2)*S_PP
         + G**2*(sym+2*M1*M2*k1k2) - 4*A*G*(k1p1*k2p2-k1p2*k2p1))
P = sp.Poly(sp.expand(val), x)
c0, c1, c2 = [sp.factor(sp.simplify(P.coeff_monomial(x**k))) for k in (0, 1, 2)]
print("c  (x^2):", c2)
print("b  (x^1):", sp.factor(c1))
for nm, ex in (("A^2", c0.coeff(A, 2).subs({B: 0, G: 0})), ("G^2", sp.expand(c0).coeff(G, 2)),
               ("A*G", sp.expand(c0).coeff(A, 1).coeff(G, 1)), ("A*B", sp.expand(c0).coeff(A, 1).coeff(B, 1)),
               ("B^2", sp.expand(c0).coeff(B, 2))):
    print(f"a[{nm}]:", sp.factor(sp.simplify(ex)))

# ---- GENIE's expression (dsigma/dQ2 bracket) with only A corrected, in the same variables
F1, F2, FA = sp.symbols('F1 F2 FA', real=True)        # GENIE F1V, xiF2V, FA (= -g1)
q2 = -Q2; Mn = M1; Mp, Mm = M2+M1, M2-M1
v = (M2**2 - Mn**2 - q2)/(2*Mn)
Acorr = F1 + Mp*F2/(2*Mn)
def genie(Avec, fix_dim=False):
    w1 = (Mm**2-q2)/(4*Mn**2)*Avec**2 + (Mp**2-q2)/(4*Mn**2)*FA**2
    w2 = FA**2 + (Avec - Mp*F2/(2*Mn))**2 - q2/Mn**2*(F2/2)**2
    w3 = 2*FA*Avec
    second = (Mn**2*v**2 - (Mn**2*q2 if fix_dim else q2))
    return sp.expand(-8*Mn**2*q2*w1 - 4*second*w2 - 2*x*q2*w3 + x**2*w2)     # antineutrino (sign=+1)
exact = sp.expand(val.subs({A: Acorr, B: F2/(2*Mn), G: -FA}))
for lab, gen in (("GENIE as is", genie(F1+F2)), ("A corrected", genie(Acorr)), ("A corrected + q2 -> Mnuc^2 q2", genie(Acorr, True))):
    d = sp.simplify(sp.expand(exact - gen).subs(ml2, 0))
    print(f"{lab:32s} exact - GENIE (ml=0) = {sp.factor(d)}")
d = sp.simplify(sp.expand(exact - genie(Acorr, True)))
print("remaining with ml != 0:", sp.factor(d))

# ---- effect of each fix on sigma(E): integrate the bracket / E^2 over Q2 (same constant everywhere)
import numpy as np, sys
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../../../pelee")); import HypXsecWeights as H
def ffs(pdg, Q2v):
    f1, f2, g1 = H.ff(pdg, Q2v); return f1, f2, -g1
variants = {"GENIE as is": genie(F1+F2), "fix 1 (A)": genie(Acorr), "fix 2 (q2 dim) only": genie(F1+F2, True),
            "fix 1+2": genie(Acorr, True), "exact (incl. m_mu)": exact}
fn = {k: sp.lambdify((M1, M2, ml2, Q2, x, F1, F2, FA), e, "numpy") for k, e in variants.items()}
def sigma(key, pdg, E, MYover=None):
    MY = H.MASS[pdg] if MYover is None else MYover; Mn_ = H.NUC_MASS[pdg]; m2 = H.MMU**2
    s_ = Mn_**2+2*Mn_*E; rs = np.sqrt(s_); El = (s_+m2-MY**2)/(2*rs); pl = np.sqrt(El**2-m2); Ein = (s_-Mn_**2)/(2*rs)
    lo, hi = -m2+2*Ein*(El-pl), -m2+2*Ein*(El+pl)
    xs_, ws_ = np.polynomial.legendre.leggauss(64); tot = 0
    for t, w in zip(xs_, ws_):
        Qv = lo + (hi-lo)*(t+1)/2; f1, f2, fa = ffs(pdg, Qv)
        if MYover is not None: pass
        Sg = Mn_**2+MY**2+m2+Qv; xv = 2*s_-Sg
        tot += w*(hi-lo)/2*max(fn[key](Mn_, MY, m2, Qv, xv, f1, f2, fa), 0)
    return tot/E**2
for pdg, nm in ((3212, "Sigma0"), (3122, "Lambda")):
    for E in (1.0, 2.0):
        ex = sigma("exact (incl. m_mu)", pdg, E)
        print(f"{nm} E={E}: " + "  ".join(f"{k}: {sigma(k, pdg, E)/ex:.3f}" for k in variants))
# equal-mass limit (final baryon at the nucleon mass, Sigma0 couplings)
ex = sigma("exact (incl. m_mu)", 3212, 1.0, MYover=0.93827)
print("M_Y = M_N, E=1: " + "  ".join(f"{k}: {sigma(k, 3212, 1.0, MYover=0.93827)/ex:.3f}" for k in variants))
