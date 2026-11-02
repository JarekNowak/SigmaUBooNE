#!/usr/bin/env python3
# ============================================================================
#  make_figures.py -- figures for GENIE_HYPERON_BUG.md: the GENIE numubar hyperon QE cross section
#  before (genie::PaisQELLambdaPXSec) and after the fix (genie::HypQELTracePXSec), and the
#  phase space the two generations cover.
#    1. sigma(E) from GENIE's own splines (merged hypfix spline file: both algorithms), Ar40 and
#       free nucleon, numubar, for Lambda, Sigma0, Sigma-; and the ratio.
#    2. dsigma/dQ2 shape at fixed E: the Pais expression (unclipped) against the trace, from
#       ../../pelee/HypXsecWeights.py (same form factors on both sides).
#    3. Event rates per POT and Sigma0 kinematics from the generated RHC run3b samples
#       (chunk.summary, chunk.hepmcev): old GENIE, fixed GENIE, and old x the cross-section weight
#       (HypXsecWeight.h table), i.e. what the analysis applies to the old samples.
#    4. Phase space: (E_nu, Q2) of the Sigma0 events with the free-nucleon kinematic limits,
#       and u = (Q2 - Q2min)/(Q2max - Q2min) in the struck-nucleon rest frame.
#  Run: python3 make_figures.py   (writes figs/*.png and prints the numbers quoted in the note)
# ============================================================================
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs"); os.makedirs(FIG, exist_ok=True)
SPL = "/data/uboone/genie_samples/hypfix_splines/gxspl-FNALsmall-hypfix.xml"
OLD = "/data/uboone/genie_samples/genie_uboone_NuMI_RHC_sigmazero_run3b_2026/chunks"
NEW = "/data/uboone/genie_samples/genie_uboone_NuMI_RHC_sigmazero_run3b_2026_hypfix/chunks"
WTAB = os.path.join(HERE, "../../pelee/hyp_xsec_weights.root")
sys.path.insert(0, os.path.join(HERE, "../../pelee"))
import HypXsecWeights as HX

GEV2_TO_1E40CM2 = 0.3894e-27 / 1e-40          # GeV^-2 -> 1e-40 cm^2
NAMES = {3122: r"$\Lambda$", 3212: r"$\Sigma^0$", 3112: r"$\Sigma^-$"}
COL = {"old": "#c0392b", "new": "#1f5fa8", "rw": "#2e8b57"}
MMU = 0.10566
MY = {3122: 1.115683, 3212: 1.192642, 3112: 1.197449}
MNUC = {3122: 0.93827, 3212: 0.93827, 3112: 0.93957}
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

# ---------------------------------------------------------------- 1. splines
def read_splines(path, want):
    """{name: (E, xsec)} for the spline names in want (streamed; the file is 490 MB)."""
    out, cur, E, X = {}, None, [], []
    rx = re.compile(r'<E>\s*([-0-9.eE+]+)\s*</E>\s*<xsec>\s*([-0-9.eE+]+)\s*</xsec>')
    with open(path) as f:
        for line in f:
            if cur is None:
                if "<spline name=" in line:
                    nm = line.split('"')[1]
                    if nm in want: cur, E, X = nm, [], []
            else:
                m = rx.search(line)
                if m: E.append(float(m.group(1))); X.append(float(m.group(2)))
                elif "</spline>" in line:
                    out[cur] = (np.array(E), np.array(X) * GEV2_TO_1E40CM2); cur = None
                    if len(out) == len(want): break
    return out

def sname(alg, tgt, pdg):
    N = 2112 if pdg == 3112 else 2212
    return f"genie::{alg}/Default/nu:-14;tgt:{tgt};N:{N};proc:Weak[CC],QES;strange:{pdg};"

FREE = {3122: 1000010010, 3212: 1000010010, 3112: 1000000010}
want = [sname(a, 1000180400, p) for a in ("PaisQELLambdaPXSec", "HypQELTracePXSec") for p in (3122, 3212, 3112)]
print("reading splines ...")
spl = read_splines(SPL, set(want))

fig, ax = plt.subplots(2, 3, figsize=(12, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
print("sigma ratio fixed/Pais (Ar40, numubar):")
for i, p in enumerate((3122, 3212, 3112)):
    Eo, So = spl[sname("PaisQELLambdaPXSec", 1000180400, p)]
    En, Sn = spl[sname("HypQELTracePXSec", 1000180400, p)]
    nN = 18 if p != 3112 else 22                                    # protons / neutrons in Ar40
    ax[0, i].plot(Eo, So / nN, color=COL["old"], label="GENIE Pais (as used)")
    ax[0, i].plot(En, Sn / nN, color=COL["new"], label="fixed (trace)")
    r = np.interp(Eo, En, Sn) / np.where(So > 0, So, np.nan)
    ax[1, i].plot(Eo, r, color="k"); ax[1, i].axhline(1, color="gray", lw=0.8)
    ax[0, i].set_title(r"$\bar\nu_\mu\,$" + ("p" if p != 3112 else "n") + r"$\to\mu^+$" + NAMES[p])
    ax[1, i].set_xlabel(r"$E_{\bar\nu}$ [GeV]"); ax[1, i].set_xlim(0, 6)
    ax[1, i].set_ylim(0.5, 1.4) if p != 3122 else ax[1, i].set_ylim(0.9, 1.5)
    print("  ", NAMES[p], " ".join(f"E={e}:{np.interp(e, Eo, r):.3f}" for e in (1, 2, 5)))
ax[0, 0].set_ylabel(r"$\sigma$ per target nucleon [$10^{-40}$ cm$^2$]"); ax[1, 0].set_ylabel("fixed / Pais")
ax[0, 0].legend(fontsize=8)
fig.suptitle(r"GENIE $\Delta S=1$ QE cross section on argon, $\bar\nu_\mu$ (G18_10a_02_11a), from the GENIE splines")
fig.tight_layout(); fig.savefig(f"{FIG}/xsec_vs_E.png", dpi=140); plt.close(fig)

# ratio Sigma0 / Lambda
Eo, SoL = spl[sname("PaisQELLambdaPXSec", 1000180400, 3122)]; _, SoS = spl[sname("PaisQELLambdaPXSec", 1000180400, 3212)]
En, SnL = spl[sname("HypQELTracePXSec", 1000180400, 3122)]; _, SnS = spl[sname("HypQELTracePXSec", 1000180400, 3212)]
print(f"  Sigma0/Lambda on Ar at 1 GeV: Pais {np.interp(1, Eo, SoS / np.where(SoL > 0, SoL, np.nan)):.3f}  "
      f"fixed {np.interp(1, En, SnS / np.where(SnL > 0, SnL, np.nan)):.3f}")

# ---------------------------------------------------------------- 2. dsigma/dQ2 shape
fig, ax = plt.subplots(2, 3, figsize=(12, 5.6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
print("Pais clipped to 0 (u range, fraction of the trace sigma there):")
for i, p in enumerate((3122, 3212, 3112)):
    for E, ls in ((1.0, "-"), (3.0, "--")):
        u = np.linspace(0.001, 0.999, 500)
        pa = np.array([HX.dxs_pais(p, E, x) / (2 * np.pi) for x in u])
        tr = []
        for x in u:
            k = HX.kinematics(p, E, x); s_, Ein, pl = k[5], k[6], k[7]
            tr.append(HX.msq_trace(p, E, x) * pl / (32 * np.pi * s_ * Ein) / 4 / (2 * Ein * pl))
        tr = np.array(tr); norm = tr.max()
        ax[0, i].plot(u, pa / norm, color=COL["old"], ls=ls, label=f"GENIE Pais, E={E:g} GeV")
        ax[0, i].plot(u, tr / norm, color=COL["new"], ls=ls, label=f"trace (fixed), E={E:g} GeV")
        ax[1, i].plot(u, np.where(tr > 0, pa / tr, np.nan), color="k", ls=ls, label=f"E={E:g} GeV")
        z = pa <= 0
        if z.any():
            ax[1, i].axvspan(u[z].min(), 1, color=COL["old"], alpha=0.15 if E == 1.0 else 0.08)
            print(f"   {NAMES[p]} E={E:g}: u > {u[z].min():.3f}, {tr[z].sum() / tr.sum():.4f}")
    ax[0, i].set_yscale("log"); ax[0, i].set_ylim(1e-4, 3); ax[0, i].set_title(NAMES[p])
    ax[1, i].axhline(1, color="gray", lw=0.8); ax[1, i].set_ylim(0, 2.2)
    ax[1, i].set_xlabel(r"$u=(Q^2-Q^2_{\min})/(Q^2_{\max}-Q^2_{\min})$")
ax[0, 0].set_ylabel(r"d$\sigma$/d$Q^2$ (trace max = 1)"); ax[1, 0].set_ylabel("Pais / trace")
ax[0, 0].legend(fontsize=7); ax[1, 0].legend(fontsize=7)
fig.suptitle(r"Free-nucleon d$\sigma$/d$Q^2$, same form factors. Shaded: Pais negative, clipped to 0 by GENIE (no events generated)")
fig.tight_layout(); fig.savefig(f"{FIG}/dxs_dq2_shape.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- 3. generated samples
def summaries(base):
    pot, cnt, n = 0.0, {3122: 0, 3212: 0, 3112: 0}, 0
    for d in sorted(glob.glob(f"{base}/c*")):
        if not os.path.exists(f"{d}/DONE"): continue
        for l in open(f"{d}/chunk.summary"):
            x = l.split()
            if x[0] == "pot": pot += float(x[1])
            elif x[0] == "primary_hyperon" and int(x[1]) in cnt: cnt[int(x[1])] += int(x[2])
        n += 1
    return pot, cnt, n

def events(base, maxchunks=None):
    """Sigma0 events, numubar CC: E_nu, Q2, p_mu, cos_mu, p_Sigma, E_rest (struck-nucleon frame)."""
    rows = []
    dirs = [d for d in sorted(glob.glob(f"{base}/c*")) if os.path.exists(f"{d}/DONE")][:maxchunks]
    for d in dirs:
        lines = open(f"{d}/chunk.hepmcev").read().split("\n"); i = 0
        while i < len(lines):
            h = lines[i].split()
            if len(h) != 2: i += 1; continue
            npart = int(h[1]); parts = [lines[i + 1 + k].split() for k in range(npart)]; i += 1 + npart
            P = {}
            for q in parts:
                st, pdg = int(q[0]), int(q[1]); v = np.array([float(q[9]), float(q[6]), float(q[7]), float(q[8])])
                if st == 0: P["nu"] = (pdg, v)
                elif st == 11: P["N"] = v
                elif st == 1 and pdg == -13: P["mu"] = v
                elif st == 1 and pdg == 3212: P["Y"] = v
            if "nu" not in P or P["nu"][0] != -14 or "mu" not in P or "Y" not in P or "N" not in P: continue
            nu, N, mu, Y = P["nu"][1], P["N"], P["mu"], P["Y"]
            q = nu - mu; Q2 = -(q[0]**2 - q[1:] @ q[1:])
            pmu = np.linalg.norm(mu[1:]); cosmu = (mu[1:] @ nu[1:]) / (pmu * np.linalg.norm(nu[1:]))
            mN = np.sqrt(max(N[0]**2 - N[1:] @ N[1:], 1e-6))
            Erest = (nu[0] * N[0] - nu[1:] @ N[1:]) / mN
            rows.append((nu[0], Q2, pmu, cosmu, np.linalg.norm(Y[1:]), Erest))
    return np.array(rows), len(dirs)

def q2_limits(p, E):
    Mn, MYp = MNUC[p], MY[p]; s = Mn * Mn + 2 * Mn * E; rs = np.sqrt(s)
    El = (s + MMU**2 - MYp**2) / (2 * rs); pl2 = El**2 - MMU**2
    pl = np.sqrt(np.clip(pl2, 0, None)); Ein = (s - Mn * Mn) / (2 * rs)
    lo, hi = -MMU**2 + 2 * Ein * (El - pl), -MMU**2 + 2 * Ein * (El + pl)
    return np.where(pl2 > 0, lo, np.nan), np.where(pl2 > 0, hi, np.nan)

# weight table (as HypXsecWeight.h: bilinear in (log10 E_rest, u), clamped, capped at 5)
wt = uproot.open(WTAB)["w_3212"]; wv, xe, ye = wt.values(), wt.axes[0].edges(), wt.axes[1].edges()
from scipy.interpolate import RegularGridInterpolator
xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
interp = RegularGridInterpolator((xc, yc), wv)
def xsec_weight(Erest, Q2):
    lo, hi = q2_limits(3212, Erest); u = np.clip((Q2 - lo) / (hi - lo), 0, 1)
    x = np.clip(np.log10(Erest), xc[0], xc[-1]); uu = np.clip(u, yc[0], yc[-1])
    return np.minimum(interp(np.stack([x, uu], 1)), 5.0), u

print("reading generated samples ...")
po, co, no = summaries(OLD); pn, cn, nn = summaries(NEW)
print(f"  old: {no} chunks, {po:.4e} POT;  fixed: {nn} chunks, {pn:.4e} POT")
print("  rate per POT fixed/old:", {NAMES[k]: round((cn[k] / pn) / (co[k] / po), 3) for k in co})
ev_o, nco = events(OLD, maxchunks=1500); ev_n, ncn = events(NEW)
# POT of the chunks actually read
pot_o = sum(float(open(f"{d}/chunk.summary").readline().split()[1]) for d in
            [d for d in sorted(glob.glob(f"{OLD}/c*")) if os.path.exists(f"{d}/DONE")][:nco])
pot_n = pn
w_o, u_o = xsec_weight(ev_o[:, 5], ev_o[:, 1]); _, u_n = xsec_weight(ev_n[:, 5], ev_n[:, 1])
SCALE = 1e21
print(f"  Sigma0 numubar events read: old {len(ev_o)} ({nco} chunks), fixed {len(ev_n)}")
print(f"  Sigma0 per 1e21 POT: old {len(ev_o) / pot_o * SCALE:.1f}, old x weight {w_o.sum() / pot_o * SCALE:.1f}, "
      f"fixed {len(ev_n) / pot_n * SCALE:.1f}")
print(f"  fraction at u > 0.8: old {np.mean(u_o > 0.8):.4f}, fixed {np.mean(u_n > 0.8):.4f}; "
      f"u > 0.9: old {np.mean(u_o > 0.9):.5f}, fixed {np.mean(u_n > 0.9):.4f}")

# rates bar chart
fig, ax = plt.subplots(figsize=(5, 3.4))
xs = np.arange(3); ks = (3122, 3212, 3112)
ax.bar(xs - 0.2, [co[k] / po * SCALE for k in ks], 0.4, color=COL["old"], label="GENIE Pais (as generated)")
ax.bar(xs + 0.2, [cn[k] / pn * SCALE for k in ks], 0.4, color=COL["new"], label="fixed")
for x, k in zip(xs, ks):
    ax.text(x, max(co[k] / po, cn[k] / pn) * SCALE * 1.03, f"x{(cn[k] / pn) / (co[k] / po):.2f}", ha="center", fontsize=9)
ax.set_xticks(xs, [NAMES[k] for k in ks]); ax.set_ylabel(r"primary hyperons per $10^{21}$ POT (generated)")
ax.set_title("NuMI RHC run3b, argon, all flavours", fontsize=9); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/rates.png", dpi=140); plt.close(fig)

# 1D kinematics
specs = [(0, r"$E_{\bar\nu}$ [GeV]", np.linspace(0, 6, 49)), (1, r"$Q^2$ [GeV$^2$]", np.linspace(0, 2.5, 51)),
         (2, r"$p_\mu$ [GeV/c]", np.linspace(0, 5, 51)), (3, r"$\cos\theta_\mu$", np.linspace(-1, 1, 41)),
         (4, r"$p_{\Sigma^0}$ [GeV/c]", np.linspace(0, 2, 41))]
fig, ax = plt.subplots(2, 3, figsize=(12, 6.4)); ax = ax.ravel()
for j, (c, lab, b) in enumerate(specs):
    for arr, w, pot, key, nm in ((ev_o[:, c], None, pot_o, "old", "GENIE Pais"), (ev_n[:, c], None, pot_n, "new", "fixed"),
                                 (ev_o[:, c], w_o, pot_o, "rw", "Pais x weight")):
        h, _ = np.histogram(arr, b, weights=w); h = h / pot * SCALE
        ax[j].step(b[:-1], h, where="post", color=COL[key], label=nm, lw=1.4 if key != "rw" else 1.1,
                   ls="-" if key != "rw" else "--")
    ax[j].set_xlabel(lab); ax[j].set_ylabel(r"$\Sigma^0$ per bin per $10^{21}$ POT")
ax[0].legend(fontsize=8)
for arr, pot, key, nm in ((u_o, pot_o, "old", "GENIE Pais"), (u_n, pot_n, "new", "fixed")):
    h, b = np.histogram(arr, np.linspace(0, 1, 41)); ax[5].step(b[:-1], h / pot * SCALE, where="post", color=COL[key], label=nm)
h, b = np.histogram(u_o, np.linspace(0, 1, 41), weights=w_o); ax[5].step(b[:-1], h / pot_o * SCALE, where="post", color=COL["rw"], ls="--")
ax[5].set_xlabel(r"$u=(Q^2-Q^2_{\min})/(Q^2_{\max}-Q^2_{\min})$, nucleon rest frame"); ax[5].set_yscale("log")
fig.suptitle(r"Generated $\bar\nu_\mu$ CC $\Sigma^0$, NuMI RHC run3b, absolute per POT")
fig.tight_layout(); fig.savefig(f"{FIG}/kinematics.png", dpi=140); plt.close(fig)

# 2D phase space: lab (E_nu, Q2) and struck-nucleon rest frame (E_rest, u)
fig, ax = plt.subplots(2, 3, figsize=(13.5, 7.6))
Eg = np.linspace(0.45, 6, 300); lo, hi = q2_limits(3212, Eg)
planes = [((ev_o[:, 0], ev_o[:, 1]), (ev_n[:, 0], ev_n[:, 1]), np.linspace(0.4, 6, 57), np.linspace(0, 3, 61),
           np.linspace(0.4, 6, 15), np.linspace(0, 3, 13), r"$E_{\bar\nu}$ [GeV]", r"$Q^2$ [GeV$^2$]", True),
          ((ev_o[:, 5], u_o), (ev_n[:, 5], u_n), np.linspace(0.4, 6, 57), np.linspace(0, 1, 51),
           np.linspace(0.4, 6, 15), np.linspace(0, 1, 11), r"$E_{\bar\nu}$ in the struck-nucleon rest frame [GeV]", r"$u$", False)]
for r, (o, n, bx, by, cx, cy, lx, ly, lab) in enumerate(planes):
    Ho, _, _ = np.histogram2d(*o, [bx, by]); Ho = Ho / pot_o * SCALE
    Hn, _, _ = np.histogram2d(*n, [bx, by]); Hn = Hn / pot_n * SCALE
    vmax = max(Ho.max(), Hn.max())
    for a, H, t in ((ax[r, 0], Ho, "GENIE Pais (before)"), (ax[r, 1], Hn, "fixed (after)")):
        m = a.pcolormesh(bx, by, np.where(H > 0, H, np.nan).T, norm=LogNorm(vmax * 1e-4, vmax), cmap="viridis")
        if lab: a.plot(Eg, lo, "w-", lw=1); a.plot(Eg, hi, "w--", lw=1); a.set_ylim(0, 3)
        a.set_title(t); a.set_xlabel(lx); a.set_ylabel(ly)
    fig.colorbar(m, ax=ax[r, 1], label=r"$\Sigma^0$ per bin per $10^{21}$ POT")
    Co, _, _ = np.histogram2d(*o, [cx, cy]); Cn, _, _ = np.histogram2d(*n, [cx, cy])
    R = np.where((Co > 0) & (Cn > 0), (Cn / pot_n) / np.where(Co > 0, Co / pot_o, 1), np.nan)
    m2 = ax[r, 2].pcolormesh(cx, cy, R.T, vmin=0.3, vmax=1.7, cmap="RdBu_r")
    if lab: ax[r, 2].plot(Eg, hi, "k--", lw=1); ax[r, 2].set_ylim(0, 3)
    ax[r, 2].set_title("after / before, per POT (coarse bins)"); ax[r, 2].set_xlabel(lx)
    fig.colorbar(m2, ax=ax[r, 2], label="ratio")
fig.suptitle(r"Phase space of the generated $\Sigma^0$ (NuMI RHC run3b). Top: lab, free-nucleon $Q^2_{\min}$ (solid), $Q^2_{\max}$ (dashed). Bottom: nucleon rest frame")
fig.tight_layout(); fig.savefig(f"{FIG}/phase_space.png", dpi=140); plt.close(fig)
print("figures in", FIG)
