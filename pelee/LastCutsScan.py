#!/usr/bin/env python3
# ============================================================================
#  LastCutsScan.py -- re-optimise the last three cuts of the Lambda-style chain (photon gap,
#  island planes, proton llrP) and the shower-count ceiling, on the expected 90% CL limit,
#  with the full signal definition (sigdef_*.tsv) and every weight:
#    GENIE signal: cache weight x hyperon xsec x polarisation (hypw_genie_qe.tsv) x PPFX 0.934
#    NuWro signal: cache weight x polarisation (pol_nuwro.tsv)
#    background:   overlay cache weight x GENIE hyperon-QE xsec weight (bkg_hypqe.tsv); since
#                  2026-11-06 its hyperon events are replaced by the hyperon detvar CV (52x data;
#                  flux capped at 3 x xsec, w_dvcv_bkg.tsv; DirtCvCheck.py)
#                  + Sigma0 outside the definition from each model's own sample
#                  + EXT, dirt and NuWro non-signal hyperons scaled with the MC (factorised)
#  Earlier stages fixed at the operating point: pair, |m - m_Lambda| < 20 MeV, cos(alpha) > 0.95,
#  decay length > 0.5 cm. Limit as Sensitivity.C (Bayesian, lognormal b and e, median n).
#  The choice is cross-validated by event parity.  Run: python3 LastCutsScan.py [f (0.15)] [1d|-] [ovl]
#  ("ovl": the overlay-only background model used before 2026-11-06)
# ============================================================================
import sys, itertools
import numpy as np, uproot
from scipy import stats

D = "/home/t2k/nowak/MicroBooNE/Sigma/pelee/"
F = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
ONLY1D = len(sys.argv) > 2 and sys.argv[2] == "1d"     # one cut at a time only, with the mean limit
OVL_ONLY = len(sys.argv) > 3 and sys.argv[3] == "ovl"
MEAN = ONLY1D
G, MINRAW = 0.20, 5
if OVL_ONLY:
    BMC_NOM, EXT_PER_BMC, DIRT_PER_BMC = 2.12, 0.62 * 0.541 / 2.12, 0.073 / 2.12
else:   # factorised EXT / dirt at the operating point, over the MC there (Sensitivity.C)
    BMC_NOM, EXT_PER_BMC, DIRT_PER_BMC = 2.332, 0.62 * 0.598 / 2.332, 0.073 / 2.332
HYP_PER_BMC = 0.0207 * 0.947 / BMC_NOM

def key(r, s, e): return (r.astype(np.int64) * 100000 + s) * 10000000 + e
def table(fn, col=3):
    m = {}
    for l in open(D + fn):
        if l.startswith('#') or not l.strip(): continue
        x = l.split(); m[(int(x[0]) * 100000 + int(x[1])) * 10000000 + int(x[2])] = float(x[col])
    return m

V = ["run", "sub", "evt", "weight", "isSig", "hasPair", "mass", "cosAlpha", "decayLen", "islNPlanes",
     "islPass", "islOK", "llrP", "phDist", "nSh", "nTr", "shrE"]
def load(fn, wtab, extra, defTab=None):
    a = uproot.open(D + fn)["isl"].arrays(V, library="np")
    if defTab:
        m = a["isSig"] == 1; a = {k: v[m] for k, v in a.items()}
    k = key(a["run"], a["sub"], a["evt"])
    t = table(wtab); a["w"] = a["weight"] * extra * np.array([t.get(int(x), 1.0) for x in k])
    if defTab:
        d = table(defTab); a["inDef"] = np.array([d.get(int(x), 0) for x in k]) > 0.5
    a["odd"] = (a["evt"] % 2).astype(bool)
    a["allOK"] = np.array([all(p[i] or not o[i] for i in range(3)) and o.sum() > 0
                           for p, o in zip(a["islPass"], a["islOK"])])
    return a
sg = load("isl_sig_genie_qe.root", "hypw_genie_qe.tsv", 0.934, "sigdef_genie_qe.tsv")
sn = load("isl_sig_nuwro.root", "pol_nuwro.tsv", 1.0, "sigdef_nuwro.tsv")
bk = load("isl_bkg.root", "bkg_hypqe.tsv", 1.0)
if not OVL_ONLY:
    HYPS = [3122, 3212, 3222, 3112, 3322, 3312, 3334]
    ht = table("ovl_hyp_truth.tsv")       # col 3 = hyperon pdg
    keep = ~np.isin(np.abs(np.array([int(ht.get(int(k), 0)) for k in key(bk["run"], bk["sub"], bk["evt"])])), HYPS)
    bk = {k: v[keep] for k, v in bk.items()}
    cv = load("isl_dvcv_bkg.root", "w_dvcv_bkg.tsv", 1.0)        # w = weight x flux (uncapped): redone
    kc = key(cv["run"], cv["sub"], cv["evt"]); wfl, wxs = table("w_dvcv_bkg.tsv", 3), table("w_dvcv_bkg.tsv", 4)
    hp = uproot.open(D + "isl_dvcv_bkg.root")["isl"]["hypPdg"].array(library="np")
    cv["w"] = cv["weight"] * 1.1526e21 / 5.9943e22 * np.array([min(wfl[int(x)], 3.0) * wxs[int(x)] for x in kc])
    cv = {k: v[np.isin(np.abs(hp), HYPS)] for k, v in cv.items()}
    bk = {k: np.concatenate([bk[k], cv[k]]) for k in bk}


GH_X, GH_W = np.polynomial.hermite_e.hermegauss(24); GH_W = GH_W / GH_W.sum()
def ul_mu(B, rb, S, g=G, n=None):
    if n is None: n = int(stats.poisson.median(B))
    sb = np.sqrt(np.log(1 + rb * rb)); bs = B * np.exp(sb * GH_X - 0.5 * sb * sb)
    se = np.sqrt(np.log(1 + g * g)); es = np.exp(se * GH_X - 0.5 * se * se)
    mu = np.linspace(0, 40 / S, 4001)
    lam = mu[:, None, None] * S * es[None, None, :] + bs[None, :, None]
    L = (stats.poisson.pmf(n, lam) * GH_W[None, :, None] * GH_W[None, None, :]).sum(axis=(1, 2))
    c = np.cumsum(L); c /= c[-1]
    return mu[np.searchsorted(c, 0.90)]

def ul_mean(B, rb, S):
    """mean expected limit over n ~ Pois(B): steadier than the median at B ~ 2-4"""
    ns = np.arange(0, int(B + 6 * np.sqrt(B) + 6)); p = stats.poisson.pmf(ns, B)
    return float(sum(pi * ul_mu(B, rb, S, n=int(k)) for k, pi in zip(ns, p) if pi > 1e-4) / p[p > 1e-4].sum())

def mask(a, c):
    m = (a["hasPair"] == 1) & (np.abs(a["mass"] - 1.115683) < 0.020) & (a["cosAlpha"] > 0.95) \
        & (a["decayLen"] > 0.5) & (a["nSh"] <= c["nsh"])
    if c["gap"] is not None: m &= a["phDist"] > c["gap"]
    m &= (a["allOK"] & (a["islNPlanes"] >= 1)) if c["npl"] == "all" else (a["islNPlanes"] >= c["npl"])
    if c["llr"] is not None: m &= a["llrP"] < c["llr"]
    return m

def evaluate(c, half=None):
    f = 1.0 if half is None else 2.0
    def sel(a, extra=None):
        m = mask(a, c)
        if half is not None: m &= a["odd"] == half
        return m if extra is None else m & extra
    ms, mso = sel(sg, sg["inDef"]), sel(sg, ~sg["inDef"])
    mn, mno = sel(sn, sn["inDef"]), sel(sn, ~sn["inDef"])
    mb = sel(bk)
    S, OG = sg["w"][ms].sum() * f, sg["w"][mso].sum() * f
    SN, ON = sn["w"][mn].sum() * f, sn["w"][mno].sum() * f
    bmc, ebmc = bk["w"][mb].sum() * f, np.sqrt((bk["w"][mb] ** 2).sum()) * f
    nraw = int(round((bmc / ebmc) ** 2)) if ebmc > 0 else 0   # effective events (CV events weigh 0.02)
    if S <= 0 or nraw < 1: return None
    B0 = bmc * (1 + EXT_PER_BMC + HYP_PER_BMC + DIRT_PER_BMC)
    rb = np.hypot(np.hypot(ebmc, 0.27 / 0.62 * EXT_PER_BMC * bmc) / B0, F)
    BN, BG = B0 + ON, B0 + OG
    r = dict(S=S, SN=SN, OG=OG, ON=ON, bmc=bmc, nraw=nraw, BN=BN, BG=BG,
             muG=ul_mu(BG, rb, S), muN=ul_mu(BN, rb, SN))
    if MEAN: r["mG"], r["mN"] = ul_mean(BG, rb, S), ul_mean(BN, rb, SN)
    return r

GRID = dict(gap=[None, 0.0, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0], npl=[1, 2, 3, "all"],
            llr=[None, 1.0, 0.5, 0.0, -0.25, -0.5, -1.0], nsh=[3, 2, 1])
NOM = dict(gap=2.0, npl=2, llr=0.0, nsh=3)
def fmt(c):
    return (f"gap>{'-' if c['gap'] is None else c['gap']} pl>={c['npl']} "
            f"llrP<{'-' if c['llr'] is None else c['llr']} nSh<={c['nsh']}")
def line(c, r):
    return (f"  {fmt(c):40s} {r['S']:.3f} {r['SN']:.4f} {r['bmc']:5.2f} {r['nraw']:4d} "
            f"{r['BN']:5.2f} {r['BG']:5.2f}  {r['muG']:5.2f} {r['muN']:6.1f}"
            + (f"   mean {r['mG']:5.2f} {r['mN']:6.1f}" if "mG" in r else ""))
HDR = f"  {'cuts':40s} S_G   S_N    B_MC  raw   B_N   B_G    mu_G   mu_N"

print(f"f = {F:.2f}, g = {G:.2f}; B = {'overlay' if OVL_ONLY else 'overlay non-hyperon + CV hyperon'} MC x (1 + EXT + dirt + hyp, factorised) + out-of-definition Sigma0")
r = evaluate(NOM); print("nominal:"); print(HDR); print(line(NOM, r))

print("\none cut at a time (others nominal):")
print(HDR)
for v in GRID:
    for x in GRID[v]:
        c = dict(NOM); c[v] = x; r = evaluate(c)
        if r: print(line(c, r) + ("   <- nominal" if c == NOM else "") + ("   (raw < 5)" if r["nraw"] < MINRAW else ""))
    print()

if ONLY1D: sys.exit(0)
res = []
for vals in itertools.product(*GRID.values()):
    c = dict(zip(GRID, vals)); r = evaluate(c)
    if r and r["nraw"] >= MINRAW: res.append((r["muG"], c, r))
res.sort(key=lambda x: x[0])
print(f"full grid, best 12 by GENIE limit (B raw >= {MINRAW}):"); print(HDR)
for _, c, r in res[:12]: print(line(c, r))
resN = sorted(res, key=lambda x: x[2]["muN"])
print(f"\nbest 8 by NuWro limit:"); print(HDR)
for _, c, r in resN[:8]: print(line(c, r))

print("\ncross-validated (choose on one parity by GENIE limit, evaluate on the other):")
for h in (False, True):
    best = None
    for vals in itertools.product(*GRID.values()):
        c = dict(zip(GRID, vals)); r = evaluate(c, half=h)
        if r is None or r["nraw"] < MINRAW: continue
        if best is None or r["muG"] < best[0]: best = (r["muG"], c)
    c = best[1]; rt, rn = evaluate(c, half=not h), evaluate(NOM, half=not h)
    print(f"  chosen on {'odd' if h else 'even'}: {fmt(c)}  -> other half mu_G {rt['muG']:.2f} "
          f"mu_N {rt['muN']:.1f} (nominal {rn['muG']:.2f} / {rn['muN']:.1f})")
