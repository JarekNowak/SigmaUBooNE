#!/usr/bin/env python3
# ============================================================================
#  DirtCvCheck.py -- the Lambda-style chain on two new background samples, stage by stage,
#  against the overlay background, with the cuts and weights of LastCutsScan.py (nominal point):
#    dirt  (isl_dirt.root): RHC run3b dirt, old flux -> flux-weighted in IslandProcess.C
#          (fw_dirt_run3b.tsv); here x data POT / dirt POT.
#    CV    (isl_dvcv_bkg.root): hyperon-enriched detvar CV (old GENIE, old flux); weights from
#          CvWeights.C (w_dvcv_bkg.tsv: flux x GENIE hyperon-QE xsec); x data POT / CV POT.
#          Compared with the overlay's OWN hyperon events (ovl_hyp_truth.tsv, HypTruthDump.py)
#          to see whether it can stand in for them at 50x the statistics.
#  Run: python3 DirtCvCheck.py <dirt POT> <CV POT>
# ============================================================================
import sys
import numpy as np, uproot

D = "/home/t2k/nowak/MicroBooNE/Sigma/pelee/"
DATA_POT = 1.1526e21
POT_DIRT, POT_CV = float(sys.argv[1]), float(sys.argv[2])
# Old -> new flux weights are capped: the map (ratios_g4_10_4_cvonly_fornuwrofakedata.root) has
# sparsely filled bins at 14-19 deg off axis, where single bins give w = 12.3-12.4. 0.12% of the CV
# events have w > 3; one of them (w 12.42) carried 17% of the CV yield at the operating point. The
# cap moves the CV normalisation by 0.6% (MuonCandidate CV/overlay-hyperon 1.004 -> 0.994).
FLUX_CAP = 3.0
HYP = {3122: "Lambda", 3212: "Sigma0", 3222: "Sigma+", 3112: "Sigma-"}

def key(r, s, e): return (r.astype(np.int64) * 100000 + s) * 10000000 + e
def table(fn, col=3, cols=None):
    m = {}
    for l in open(fn if fn.startswith("/") else D + fn):
        if l.startswith('#') or not l.strip(): continue
        x = l.split(); k = (int(x[0]) * 100000 + int(x[1])) * 10000000 + int(x[2])
        m[k] = tuple(x[c] for c in cols) if cols else float(x[col])
    return m

V = ["run", "sub", "evt", "weight", "hasPair", "mass", "cosAlpha", "decayLen", "islNPlanes",
     "llrP", "phDist", "nSh", "mode"]
def load(fn, extra=()):
    a = uproot.open(D + fn)["isl"].arrays(V + list(extra), library="np")
    a["k"] = key(a["run"], a["sub"], a["evt"]); return a

# cumulative stages of the nominal chain (LastCutsScan NOM: gap 2, 2 planes, llrP < 0, nSh <= 3)
STAGES = [("MuonCandidate", lambda a: np.ones(len(a["run"]), bool)),
          ("Pair",          lambda a: a["hasPair"] == 1),
          ("LambdaMass",    lambda a: np.abs(a["mass"] - 1.115683) < 0.020),
          ("Alpha+Length",  lambda a: (a["cosAlpha"] > 0.95) & (a["decayLen"] > 0.5)),
          ("PhotonGap",     lambda a: (a["nSh"] <= 3) & (a["phDist"] > 2.0)),
          ("Island2",       lambda a: a["islNPlanes"] >= 2),
          ("PID",           lambda a: a["llrP"] < 0.0)]
def flow(a, w, sel=None):
    m = np.ones(len(w), bool) if sel is None else sel.copy(); out = []
    for _, f in STAGES:
        m &= f(a); out.append((int(m.sum()), float(w[m].sum()), float(np.sqrt((w[m] ** 2).sum()))))
    return out

# ── overlay, as LastCutsScan: cache weight x GENIE hyperon-QE xsec weight ──────────────────
bk = load("isl_bkg.root")
xs = table("bkg_hypqe.tsv"); bk["w"] = bk["weight"] * np.array([xs.get(int(k), 1.0) for k in bk["k"]])
ht = table("ovl_hyp_truth.tsv", cols=(3,)); bk["hyp"] = np.array([int(ht[k][0]) if k in ht else 0 for k in bk["k"]])
isHypO = np.isin(np.abs(bk["hyp"]), list(HYP) + [3322, 3312, 3334])
# ── dirt ─────────────────────────────────────────────────────────────────────────────────────
dt = load("isl_dirt.root"); dt["w"] = np.minimum(dt["weight"], FLUX_CAP) * DATA_POT / POT_DIRT  # weight = flux w
# ── CV ───────────────────────────────────────────────────────────────────────────────────────
cv = load("isl_dvcv_bkg.root", extra=["hypPdg", "nuPdg", "ccnc"])
cw = table("w_dvcv_bkg.tsv", cols=(3, 4, 5))
wf = np.minimum(np.array([float(cw[k][0]) for k in cv["k"]]), FLUX_CAP); wx = np.array([float(cw[k][1]) for k in cv["k"]])
cv["w"] = cv["weight"] * wf * wx * DATA_POT / POT_CV
_, first = np.unique(cv["k"], return_index=True)
cv["uniq"] = np.zeros(len(cv["k"]), bool); cv["uniq"][first] = True
isHypC = np.isin(np.abs(cv["hypPdg"]), list(HYP) + [3322, 3312, 3334])

print(f"\n  data {DATA_POT:.4e} POT; dirt {POT_DIRT:.4e} (x{DATA_POT/POT_DIRT:.4f}); CV {POT_CV:.4e} (x{DATA_POT/POT_CV:.5f})")
print(f"  CV at MuonCandidate: {len(cv['k'])} entries, {len(first)} distinct (run,sub,evt); "
      f"{100*isHypC.mean():.1f}% with a hyperon; flux <w> {wf.mean():.3f}")

rows = [("overlay, all", flow(bk, bk["w"])),
        ("overlay, non-hyperon", flow(bk, bk["w"], ~isHypO)),
        ("overlay, hyperon", flow(bk, bk["w"], isHypO)),
        ("CV, hyperon", flow(cv, cv["w"], isHypC)),
        ("CV, hyperon, distinct", flow(cv, cv["w"], isHypC & cv["uniq"])),
        ("dirt", flow(dt, dt["w"]))]
print("\n  weighted yield per data exposure (raw events in brackets)")
print("  " + f"{'stage':<15}" + "".join(f"{n:>24}" for n, _ in rows))
for i, (s, _) in enumerate(STAGES):
    print("  " + f"{s:<15}" + "".join(f"{r[i][1]:>15.3f} ({r[i][0]:>6})" for _, r in rows))
print("\n  CV / overlay-hyperon and dirt / overlay-all, per stage")
for i, (s, _) in enumerate(STAGES):
    o, c, a, d = rows[2][1][i], rows[3][1][i], rows[0][1][i], rows[5][1][i]
    rc = c[1] / o[1] if o[1] > 0 else np.nan
    erc = rc * np.hypot(o[2] / o[1], c[2] / c[1]) if o[1] > 0 and c[1] > 0 else np.nan
    rd = d[1] / a[1] if a[1] > 0 else np.nan
    erd = rd * np.hypot(a[2] / a[1], d[2] / d[1]) if a[1] > 0 and d[1] > 0 else np.nan
    print(f"  {s:<15} CV/ovl-hyp {rc:6.3f} +- {erc:5.3f}    dirt/ovl {rd:6.4f} +- {erd:6.4f}")

# composition at the Pair+mass stage: by hyperon species and channel
def comp(a, w, sel, hcol):
    m = sel & STAGES[1][1](a) & STAGES[2][1](a); tot = w[m].sum(); out = {}
    for h in HYP:
        for mo, mn in ((0, "QE"), (1, "RES"), (2, "DIS")):
            mm = m & (np.abs(a[hcol]) == h) & (a["mode"] == mo)
            if mm.any(): out[(HYP[h], mn)] = (w[mm].sum() / tot, int(mm.sum()))
    return out
co, cc = comp(bk, bk["w"], isHypO, "hyp"), comp(cv, cv["w"], isHypC, "hypPdg")
print("\n  hyperon composition after the Lambda mass cut (fraction of weighted yield, raw)")
for k in sorted(set(co) | set(cc), key=lambda k: -cc.get(k, (0, 0))[0]):
    print(f"  {k[0]:>7} {k[1]:<4}  overlay {co.get(k,(0,0))[0]:6.3f} ({co.get(k,(0,0))[1]:>4})   CV {cc.get(k,(0,0))[0]:6.3f} ({cc.get(k,(0,0))[1]:>5})")

# for Sensitivity.C: dirt factorised from the Lambda-mass stage; background MC with the CV in
# place of the overlay's own hyperon events
op, dm, om = rows[5][1][-1], rows[5][1][2], rows[0][1]
print(f"\n  dirt at the operating point: {op[1]:.3f} +- {op[2]:.3f} ({op[0]} raw); factorised from the mass"
      f" cut: {dm[1]:.2f} +- {dm[2]:.2f} x {om[-1][1]:.3f}/{om[2][1]:.1f} = {dm[1]*om[-1][1]/om[2][1]:.4f}")
nh, ch = rows[1][1], rows[3][1]
B, eB = nh[-1][1] + ch[-1][1], np.hypot(nh[-1][2], ch[-1][2]); B2 = nh[-2][1] + ch[-2][1]
print(f"  B_MC = overlay non-hyperon {nh[-1][1]:.3f} +- {nh[-1][2]:.3f} ({nh[-1][0]} raw) + CV hyperon"
      f" {ch[-1][1]:.3f} +- {ch[-1][2]:.3f} ({ch[-1][0]} raw) = {B:.3f} +- {eB:.3f} ({100*eB/B:.1f}%)"
      f"   [overlay alone {om[-1][1]:.3f} +- {om[-1][2]:.3f}]")
print(f"  PID pass fraction of this B_MC (for the factorised EXT): {B:.3f} / {B2:.3f} = {B/B2:.3f}")
