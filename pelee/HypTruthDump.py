#!/usr/bin/env python3
# ============================================================================
#  HypTruthDump.py -- per-event truth of every event with a hyperon (mc_hyperon_pdg a hyperon code; -999 = none) in a
#  hypana file list: "run sub evt hyperon_pdg nu_pdg ccnc mode q2". Reads only these branches, so
#  it is quick even on the full overlay list. Used to split the overlay background into its
#  hyperon and non-hyperon parts when the island cache (isl_bkg.root) predates the hypPdg branch.
#  Run: python3 HypTruthDump.py lists/bkg_overlay.list ovl_hyp_truth.tsv
# ============================================================================
import sys, uproot, numpy as np
lst, out = sys.argv[1], sys.argv[2]
B = ["run", "subrun", "event", "mc_hyperon_pdg", "mc_nu_pdg", "mc_ccnc", "mc_mode", "mc_nu_q2"]
files = [l.strip() for l in open(lst) if l.strip()]
nEv = nHyp = bad = 0
with open(out, "w") as fo:
    fo.write(f"# run sub evt hyperon_pdg nu_pdg ccnc mode q2  -- HypTruthDump.py {lst}\n")
    for i, f in enumerate(files):
        try:
            a = uproot.open(f)["ana/OutputTree"].arrays(B, library="np")
        except Exception:
            bad += 1; continue
        nEv += len(a["run"]); m = np.isin(np.abs(a["mc_hyperon_pdg"]), [3122, 3212, 3222, 3112, 3322, 3312, 3334]); nHyp += int(m.sum())
        for j in np.nonzero(m)[0]:
            fo.write(f"{a['run'][j]} {a['subrun'][j]} {a['event'][j]} {a['mc_hyperon_pdg'][j]} {a['mc_nu_pdg'][j]} "
                     f"{a['mc_ccnc'][j]} {a['mc_mode'][j]} {a['mc_nu_q2'][j]:.5f}\n")
        if i % 1000 == 0: print(f"  {i}/{len(files)}  events {nEv}  hyperon {nHyp}", flush=True)
print(f"  {len(files)} files ({bad} unreadable), {nEv} events, {nHyp} with a hyperon ({100.0*nHyp/max(nEv,1):.2f}%)")
