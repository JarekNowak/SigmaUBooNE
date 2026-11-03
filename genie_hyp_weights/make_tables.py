#!/usr/bin/env python3
# ============================================================================
#  HypXsecWeights.py -- cross-section weights that correct GENIE's numubar
#  hyperon QE model (genie::PaisQELLambdaPXSec, v3_00_06) to the full Dirac-trace
#  matrix element with the same SU(3) form factors.
#
#  GENIE's Pais formula reproduces Llewellyn Smith CCQE when the final baryon has
#  the nucleon mass (0.5% here), but its hyperon-nucleon mass-difference terms are
#  wrong: at 2 GeV it gives Lambda x0.82 and Sigma0 x1.49 of the trace result,
#  and Sigma0/Lambda 0.45 against 0.19 (NuWro's Singh model gives 0.21 on Ar).
#
#  Weight w(E, u) = [dsigma/dQ2]_trace / [dsigma/dQ2]_Pais for a free nucleon, with
#    E = antineutrino energy in the struck-nucleon rest frame (GENIE's ProbeE(kRfHitNucRest)),
#    u = (Q2 - Q2min) / (Q2max - Q2min) at that E, i.e. (1 - cos theta_CM)/2.
#  Both sides use the same form factors (dipole vector, G_En = 0, M_V^2 0.71, M_A 0.99,
#  g_A 1.267, F 0.463, D 0.804), so the weight carries only the formula difference.
#
#  Also the hyperon polarisation P_L, P_T(E, u) from the same matrix element (spin projector
#  (1 + g5 s)/2 on the hyperon; P_N = 0 without T-odd form factors), in the hyperon rest frame
#  reached by a pure boost from the nucleon rest frame: L along the hyperon, T in the (nubar, Y)
#  plane towards the nubar. GENIE and Geant4 decay the hyperons isotropically; PolWeights.C uses
#  these tables to weight the Sigma0 -> Lambda gamma -> N pi decay angles.
#  Checks: projector against explicit helicity spinors; point-like V-A gives |P| = 1 and P_L = -0.90
#  at u = 0.5, E = 50 GeV (helicity -1 is forbidden only in the exactly-backward corner).
#
#  Output: hyp_xsec_weights.root with TH2D w_<pdg> (weight), PL_<pdg>, PT_<pdg> for 3122, 3212,
#  3112 (x = log10 E/GeV, y = u), read by HypXsecWeight.h and PolWeights.C.
#  Run:    python3 HypXsecWeights.py [out.root]
# ============================================================================
import sys
import numpy as np

MN, MMU = 0.93827, 0.10566
MV2, MA, GA = 0.71, 0.99, 1.267
KP, KN = 1.793, -1.913
XFD = 0.463 / (0.463 + 0.804)
MASS = {3122: 1.115683, 3212: 1.192642, 3112: 1.197449}
NUC_MASS = {3122: 0.93827, 3212: 0.93827, 3112: 0.93957}


def ff(pdg, Q2):
    """SU(3) Cabibbo form factors for N -> Y: (f1, f2, g1) with Gamma = f1 g^mu + i f2/(2M) sig q - g1 g^mu g5."""
    tau = Q2 / (4 * MN**2)
    GD = 1 / (1 + Q2 / MV2)**2
    GEp, GMp, GEn, GMn = GD, (1 + KP) * GD, 0.0, KN * GD
    F1p = (GEp + tau * GMp) / (1 + tau); F2p = (GMp - GEp) / (1 + tau)
    F1n = (GEn + tau * GMn) / (1 + tau); F2n = (GMn - GEn) / (1 + tau)
    FA = GA / (1 + Q2 / MA**2)**2
    if pdg == 3122:
        return -np.sqrt(1.5) * F1p, -np.sqrt(1.5) * F2p, -(1 + 2 * XFD) / np.sqrt(6) * FA
    if pdg == 3212:
        return -(F1p + 2 * F1n) / np.sqrt(2), -(F2p + 2 * F2n) / np.sqrt(2), (1 - 2 * XFD) / np.sqrt(2) * FA
    return -(F1p + 2 * F1n), -(F2p + 2 * F2n), (1 - 2 * XFD) * FA


# --- Dirac matrices (Bjorken-Drell) ---
_I2, _Z = np.eye(2), np.zeros((2, 2))
_S = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
G = [np.block([[_I2, _Z], [_Z, -_I2]]).astype(complex)] + [np.block([[_Z, s], [-s, _Z]]) for s in _S]
G5 = 1j * G[0] @ G[1] @ G[2] @ G[3]
MET = np.diag([1., -1, -1, -1])
ONE = np.eye(4)
SIG = [[0.5j * (G[m] @ G[n] - G[n] @ G[m]) for n in range(4)] for m in range(4)]


def slash(p):
    return sum(MET[m, m] * p[m] * G[m] for m in range(4))


def kinematics(pdg, E, u):
    """CM-frame four-vectors for nubar(E in nucleon rest frame) + N -> mu+ Y at u; None if closed."""
    MY, Mn = MASS[pdg], NUC_MASS[pdg]
    s = Mn**2 + 2 * Mn * E; rs = np.sqrt(s)
    El = (s + MMU**2 - MY**2) / (2 * rs); pl2 = El**2 - MMU**2
    if pl2 <= 0:
        return None
    pl, Ein = np.sqrt(pl2), (s - Mn**2) / (2 * rs)
    c = 1 - 2 * u; st = np.sqrt(max(0., 1 - c * c))
    k1 = np.array([Ein, 0, 0, Ein]); p1 = np.array([rs - Ein, 0, 0, -Ein])
    k2 = np.array([El, pl * st, 0, pl * c]); p2 = k1 + p1 - k2
    q = k1 - k2
    return k1, p1, k2, p2, -(q @ MET @ q), s, Ein, pl


def msq_trace(pdg, E, u, spin=None, coup=None):
    """Spin-summed L_{mu nu} H^{mu nu} for nubar p -> mu+ Y (validated on CCQE: nu 1.06, nubar 0.41 e-38 cm2
    at 1 GeV). L from vbar(k1) g^mu (1-g5) v(k2), H_{mu nu} = sum (ubar2 G_mu u1)(ubar2 G_nu u1)^*.
    spin: hyperon spin four-vector (projector (1 + g5 sslash)/2); k is then a function of the CM kinematics
    and spin(k) must return the four-vector in the same (CM) frame. coup overrides ff() (tests)."""
    k = kinematics(pdg, E, u)
    if k is None:
        return 0.
    k1, p1, k2, p2, Q2, s, Ein, pl = k
    f1, f2, g1 = coup if coup is not None else ff(pdg, Q2)
    ql = MET @ (k1 - k2)
    Gam = [f1 * G[m] + 1j * f2 / (2 * MN) * sum(SIG[m][n] * ql[n] for n in range(4)) - g1 * G[m] @ G5 for m in range(4)]
    Gbar = [G[0] @ g.conj().T @ G[0] for g in Gam]
    P1, P2 = slash(p1) + NUC_MASS[pdg] * ONE, slash(p2) + MASS[pdg] * ONE
    if spin is not None:
        P2 = P2 @ (ONE + G5 @ slash(spin(k))) / 2
    K1, K2 = slash(k1), slash(k2) - MMU * ONE
    tot = 0
    for m in range(4):
        Lm = G[m] @ (ONE - G5) @ K2
        for n in range(4):
            L = np.trace(Lm @ G[n] @ (ONE - G5) @ K1)
            H = np.trace(P2 @ Gam[m] @ P1 @ Gbar[n])
            tot += MET[m, m] * MET[n, n] * L * H
    return tot.real


def dxs_pais(pdg, E, u):
    """GENIE PaisQELLambdaPXSec::XSec (dsigma/dQ2), transcribed; same form factors as ff()."""
    k = kinematics(pdg, E, u)
    if k is None:
        return 0.
    q2 = -k[4]
    f1, f2, g1 = ff(pdg, -q2)
    F1V, xiF2V, FA = f1, f2, -g1          # GENIE's FA0 = -1.267 convention
    M1, M2 = NUC_MASS[pdg], MASS[pdg]; Mn2 = M1**2
    v = (M2**2 - Mn2 - q2) / (2 * M1); s = Mn2 + 2 * M1 * E; uu = Mn2 + MMU**2 + 2 * v * M1 - 2 * M1 * E
    Mp, Mm = M2 + M1, M2 - M1
    w1 = (Mm**2 - q2) / (4 * Mn2) * (F1V + xiF2V)**2 + (Mp**2 - q2) / (4 * Mn2) * FA**2
    w2 = FA**2 + (F1V + xiF2V - Mp * xiF2V / (2 * M1))**2 - q2 / Mn2 * (xiF2V / 2)**2
    w3 = 2 * FA * (F1V + xiF2V)
    sign = 1   # antineutrino
    xs = (-8 * Mn2 * q2 * w1 - 4 * (Mn2 * v * v - q2) * w2 - sign * 2 * (s - uu) * q2 * w3 + (s - uu)**2 * w2) / (16 * Mn2 * E * E)
    return max(xs, 0.)


def boostz(v, beta):
    g = 1 / np.sqrt(1 - beta * beta)
    return np.array([g * (v[0] + beta * v[3]), v[1], v[2], g * (v[3] + beta * v[0])])


def polarization(pdg, E, u, coup=None):
    """Hyperon polarisation (P_L, P_T, P_N) in its rest frame reached by a pure boost from the
    struck-nucleon rest frame. L = hyperon direction there, T = in the (nubar, Y) plane towards
    the nubar, N = L x T. P_a = [H(s_a) - H(-s_a)] / [H(s_a) + H(-s_a)]."""
    k = kinematics(pdg, E, u)
    if k is None:
        return 0., 0., 0.
    MY, Mn = MASS[pdg], NUC_MASS[pdg]
    beta = E / (E + Mn)                              # CM velocity in the nucleon rest frame
    k1l, p2l = boostz(k[0], beta), boostz(k[3], beta)
    pv = p2l[1:]; pm = np.linalg.norm(pv); eL = pv / pm
    kh = k1l[1:] / np.linalg.norm(k1l[1:])
    t = kh - (kh @ eL) * eL
    eT = t / np.linalg.norm(t) if np.linalg.norm(t) > 1e-12 else np.array([1., 0, 0])
    eN = np.cross(eL, eT)
    out = []
    for e in (eL, eT, eN):
        pe = pv @ e
        sl = np.concatenate([[pe / MY], e + pe * pv / (MY * (p2l[0] + MY))])
        sc = boostz(sl, -beta)
        a = msq_trace(pdg, E, u, spin=lambda kk, sc=sc: sc, coup=coup)
        b = msq_trace(pdg, E, u, spin=lambda kk, sc=sc: -sc, coup=coup)
        out.append((a - b) / (a + b) if a + b > 0 else 0.)
    return tuple(out)


# |M|^2 -> dsigma/dQ2 relative normalisation: dsigma/dcos = |M|^2 pl /(32 pi s Ein) * GF^2 sin^2/2 * 1/2 (spin),
# dsigma/dQ2 = dsigma/dcos / (2 Ein pl); Pais: GF^2 sin^2/(2 pi) * xs.  Common GF^2 sin^2 dropped.
def weight(pdg, E, u):
    k = kinematics(pdg, E, u)
    if k is None:
        return 1.
    s, Ein, pl = k[5], k[6], k[7]
    tr = msq_trace(pdg, E, u) * pl / (32 * np.pi * s * Ein) / 2 / 2 / (2 * Ein * pl)
    pa = dxs_pais(pdg, E, u) / (2 * np.pi)
    return tr / pa if pa > 0 else 1.


def total(fn, pdg, E, n=48):
    """sigma(E) up to the common constant, integrating over u (dQ2 = 4 Ein pl du)."""
    xs, ws = np.polynomial.legendre.leggauss(n)
    tot = 0
    for x, w in zip(xs, ws):
        u = 0.5 * (x + 1)
        k = kinematics(pdg, E, u)
        if k is None:
            return 0.
        Ein, pl, s = k[6], k[7], k[5]
        if fn == "pais":
            d = dxs_pais(pdg, E, u) / (2 * np.pi)
        else:
            d = msq_trace(pdg, E, u) * pl / (32 * np.pi * s * Ein) / 4 / (2 * Ein * pl)
        tot += 0.5 * w * d * 4 * Ein * pl
    return tot


if __name__ == "__main__":
    # usage: python3 make_tables.py [out.root] [lepton: mu | e]
    out = sys.argv[1] if len(sys.argv) > 1 else "hyp_xsec_weights.root"
    if len(sys.argv) > 2 and sys.argv[2] == "e":
        MMU = 0.000511          # charged-lepton mass used in the kinematics and the trace (nubar_e)
    import ROOT
    nx, xlo, xhi = 60, np.log10(0.2), np.log10(20.)
    ny = 40
    f = ROOT.TFile(out, "RECREATE")
    for pdg in (3122, 3212, 3112):
        h = ROOT.TH2D(f"w_{pdg}", f"trace/Pais weight, Y={pdg};log_{{10}}(E_{{#nu}}^{{rest}}/GeV);u=(Q^{{2}}-Q^{{2}}_{{min}})/(Q^{{2}}_{{max}}-Q^{{2}}_{{min}})",
                      nx, xlo, xhi, ny, 0, 1)
        hs = ROOT.TH1D(f"sigratio_{pdg}", f"#sigma_{{trace}}/#sigma_{{Pais}}, Y={pdg};log_{{10}}(E/GeV)", nx, xlo, xhi)
        for i in range(1, nx + 1):
            E = 10**h.GetXaxis().GetBinCenter(i)
            for j in range(1, ny + 1):
                h.SetBinContent(i, j, weight(pdg, E, h.GetYaxis().GetBinCenter(j)))
            a, b = total("trace", pdg, E), total("pais", pdg, E)
            hs.SetBinContent(i, a / b if b > 0 else 1.)
        hL = ROOT.TH2D(f"PL_{pdg}", f"P_L, Y={pdg}", nx, xlo, xhi, ny, 0, 1)
        hT = ROOT.TH2D(f"PT_{pdg}", f"P_T (towards the #bar{{#nu}}), Y={pdg}", nx, xlo, xhi, ny, 0, 1)
        for i in range(1, nx + 1):
            E = 10**h.GetXaxis().GetBinCenter(i)
            for j in range(1, ny + 1):
                PLv, PTv, _ = polarization(pdg, E, h.GetYaxis().GetBinCenter(j))
                hL.SetBinContent(i, j, PLv); hT.SetBinContent(i, j, PTv)
        h.Write(); hs.Write(); hL.Write(); hT.Write()
        for E in (0.75, 1.0, 2.0, 5.0):
            a, b = total("trace", pdg, E), total("pais", pdg, E)
            print(f"  Y={pdg}  E={E:4.2f}  sigma trace/Pais = {a / b:.3f}")
    f.Close()
    print("wrote", out)
