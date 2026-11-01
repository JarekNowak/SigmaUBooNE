#!/usr/bin/env python3
# ============================================================================
#  IslandOptimise.py -- cut optimisation of the Lambda-style (island) chain on the
#  expected 90% CL limit, from the island caches (IslandProcess.C), with every weight:
#    GENIE signal: cache weight x hyperon xsec x Sigma0 polarisation (hypw_genie_qe.tsv) x PPFX 0.934
#    NuWro signal: cache weight x polarisation (pol_nuwro.tsv)
#    background:   cache weight x GENIE hyperon-QE xsec weight (bkg_hypqe.tsv)
#  EXT (0.62 at the 2026-10-17 operating point, factorised) and the NuWro non-signal hyperons
#  (0.021) are scaled with the MC background relative to the operating point.
#
#  Expected limit: Bayesian 90% CL on mu (flat prior), lognormal background with width
#  sqrt(rMC^2 + f^2), lognormal signal efficiency with width g = 20%, at the median n of
#  Pois(B) (as Sensitivity.C). B raw >= MINRAW is required so the scan cannot pick a cut set
#  whose background rests on a handful of events; the choice is made on even events and
#  checked on odd ones (and vice versa).
#
#  Run:  python3 IslandOptimise.py [f (0.15)]
# ============================================================================
import sys, itertools
import numpy as np
import uproot
from scipy import stats

D = "/home/t2k/nowak/MicroBooNE/Sigma/pelee/"
F = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
G = 0.20
MINRAW = 5
EXT0, HYPNW0, BMC0 = 0.62, 0.021, 3.92        # current operating point (Sensitivity.C)

def key(r, s, e): return (r.astype(np.int64) * 100000 + s) * 10000000 + e

def table(fn, col=3):
    m = {}
    for l in open(fn):
        if l.startswith('#') or not l.strip(): continue
        x = l.split(); m[(int(x[0]) * 100000 + int(x[1])) * 10000000 + int(x[2])] = float(x[col])
    return m

VARS = ["run", "sub", "evt", "weight", "isSig", "mass", "cosAlpha", "decayLen", "islNPlanes", "llrP",
        "phDist", "sigmaMass", "shrE", "hasPair"]
def load(fn, wtab=None, extra=1.0, sig=None):
    a = uproot.open(D + fn)["isl"].arrays(VARS, library="np")
    if sig is not None:
        m = a["isSig"] == sig; a = {k: v[m] for k, v in a.items()}
    k = key(a["run"], a["sub"], a["evt"])
    w = a["weight"] * extra
    if wtab is not None:
        t = table(D + wtab); w = w * np.array([t.get(int(x), 1.0) for x in k])
    a["w"] = w; a["odd"] = (a["evt"] % 2).astype(bool)
    return a

sg = load("isl_sig_genie_qe.root", "hypw_genie_qe.tsv", 0.934, sig=1)
sn = load("isl_sig_nuwro.root", "pol_nuwro.tsv", 1.0, sig=1)
bk = load("isl_bkg.root", "bkg_hypqe.tsv", 1.0)

# --- expected limit -----------------------------------------------------------
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(24); GH_W = GH_W / GH_W.sum()
SGRID = np.linspace(0, 80, 3201)
def ul_mu(B, rb, S, g=G):
    """median expected 90% CL limit on mu for signal S, background B (frac. width rb)"""
    n = int(stats.poisson.median(B))
    sb = np.sqrt(np.log(1 + rb * rb)); bs = B * np.exp(sb * GH_X - 0.5 * sb * sb)
    se = np.sqrt(np.log(1 + g * g)); es = np.exp(se * GH_X - 0.5 * se * se)
    # likelihood in mu: average over b and e nodes
    mu = SGRID / S
    lam = mu[:, None, None] * S * es[None, None, :] + bs[None, :, None]
    L = (stats.poisson.pmf(n, lam) * GH_W[None, :, None] * GH_W[None, None, :]).sum(axis=(1, 2))
    c = np.cumsum(L); c /= c[-1]
    return mu[np.searchsorted(c, 0.90)]

def yields(sel_s, sel_n, sel_b, half=None):
    def pick(a, m):
        if half is None: return m
        return m & (a["odd"] == half)
    ms, mn, mb = pick(sg, sel_s), pick(sn, sel_n), pick(bk, sel_b)
    f = 2.0 if half is not None else 1.0
    S = sg["w"][ms].sum() * f; SN = sn["w"][mn].sum() * f
    bmc = bk["w"][mb].sum() * f; nraw = int(mb.sum())
    B = bmc * (1 + (EXT0 + HYPNW0) / BMC0)
    rb = np.hypot(1 / np.sqrt(max(nraw, 1)), F)
    return S, SN, B, nraw, rb

def mask(a, c):
    m = (a["hasPair"] == 1) & (np.abs(a["mass"] - 1.115683) < c["mwin"]) & (a["cosAlpha"] > c["cosA"]) \
        & (a["decayLen"] > c["dlen"]) & (a["phDist"] > 2) & (a["islNPlanes"] >= c["npl"])
    if c["llrP"] is not None: m &= a["llrP"] < c["llrP"]
    if c["phmax"] is not None: m &= (a["phDist"] < c["phmax"]) & (a["shrE"] > 0)
    if c["smass"]: m &= (a["sigmaMass"] > 1.10) & (a["sigmaMass"] < 1.25)
    return m

GRID = dict(mwin=[0.015, 0.020, 0.030, 0.040], cosA=[0.90, 0.95, 0.98], dlen=[0.3, 0.5, 1.0],
            npl=[1, 2], llrP=[None, 0.0, -0.2], phmax=[None, 80.0], smass=[False, True])
cuts = [dict(zip(GRID, v)) for v in itertools.product(*GRID.values())]
NOM = dict(mwin=0.020, cosA=0.95, dlen=0.5, npl=2, llrP=None, phmax=None, smass=False)

def evaluate(c, half=None):
    S, SN, B, nraw, rb = yields(mask(sg, c), mask(sn, c), mask(bk, c), half)
    if S <= 0: return None
    return dict(S=S, SN=SN, B=B, nraw=nraw, rb=rb, mu=ul_mu(B, rb, S), muN=ul_mu(B, rb, SN, ) if SN > 0 else np.inf)

def fmt(c):
    return (f"m+-{c['mwin']*1000:.0f} cosA>{c['cosA']:.2f} L>{c['dlen']:.1f} pl>={c['npl']} "
            f"llrP<{'-' if c['llrP'] is None else c['llrP']} phD<{'-' if c['phmax'] is None else int(c['phmax'])} "
            f"mSig{'Y' if c['smass'] else '-'}")

nom = evaluate(NOM)
print(f"f = {F:.2f}, g = {G:.2f}, background raw >= {MINRAW}")
print(f"nominal  {fmt(NOM)}:  S_G {nom['S']:.3f}  S_N {nom['SN']:.4f}  B {nom['B']:.2f} ({nom['nraw']} raw)  "
      f"mu_UL GENIE {nom['mu']:.2f}  NuWro {nom['muN']:.1f}")

res = []
for c in cuts:
    r = evaluate(c)
    if r is None or r["nraw"] < MINRAW: continue
    res.append((r["mu"], c, r))
res.sort(key=lambda x: x[0])
print(f"\n{len(res)} cut sets with >= {MINRAW} raw background; best 15 on the full sample (GENIE mu):")
print(f"  {'cuts':72s}  S_G    S_N     B   raw  mu_G   mu_N")
for mu, c, r in res[:15]:
    print(f"  {fmt(c):72s} {r['S']:.3f} {r['SN']:.4f} {r['B']:5.2f} {r['nraw']:4d} {mu:5.2f} {r['muN']:6.1f}")

# cross-validation: choose on one half, score on the other
print("\ncross-validated choice (choose on one parity, evaluate on the other):")
for h in (False, True):
    best = None
    for c in cuts:
        r = evaluate(c, half=h)
        if r is None or r["nraw"] < MINRAW: continue
        if best is None or r["mu"] < best[0]: best = (r["mu"], c)
    c = best[1]; rt = evaluate(c, half=not h); rn = evaluate(NOM, half=not h)
    print(f"  chosen on {'odd' if h else 'even'}: {fmt(c)}")
    print(f"     on the other half: mu_G {rt['mu']:.2f} (nominal {rn['mu']:.2f}), S_G {rt['S']:.3f} B {rt['B']:.2f} ({rt['nraw']} raw)")
