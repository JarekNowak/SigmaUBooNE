# GENIE hyperon QE cross-section weights

Per-event weights that correct GENIE's antineutrino ΔS = 1 quasi-elastic hyperon cross section,
for productions already generated with GENIE as distributed.

The channels are ν̄ p → l⁺ Λ, ν̄ p → l⁺ Σ⁰ and ν̄ n → l⁺ Σ⁻ (`QEL-CC-LAMBDA`,
`genie::PaisQELLambdaPXSec`).

*J. Nowak (MicroBooNE Σ⁰ analysis), 2026-11. Self-contained: ROOT only.*

## Why

`PaisQELLambdaPXSec::XSec` has two errors when the final baryon is heavier than the nucleon:

1. **The γ^μ coefficient** is `(F1V + xiF2V)`, valid only for equal masses. It should be
   `F1V + xiF2V·(M_N+M_Y)/2M_N`.
2. **The `w2` term** of `xsec` is missing a factor M².

Against the full contraction L_μν H^μν, computed with the same form factors, GENIE gives:
- Λ: ×0.79 (1 GeV) to ×0.84 (5 GeV);
- Σ⁰ and Σ⁻: ×1.52 to ×1.43.

The Q² shape is also wrong, and the expression is clipped to zero near Q²_max.

The bug is present in every release checked: v3_00_04, v3_00_06, v3_06_02 and `master` (2026-11).
It affects the tunes in use, including G18_10a_02_11a (MicroBooNE) and AR23_20i (SBND), which both
use `PaisQELLambdaPXSec` for `QEL-CC-LAMBDA`. The full write-up, figures and the derivation of the
faulty lines are in the SigmaUBooNE repository, `genie_hypfix/doc/`. The GENIE issue text is
`genie_hypfix/doc/GENIE_ISSUE.md`.

## What the weight is

w = [dσ/dQ²]_correct / [dσ/dQ²]_GENIE for a free nucleon, tabulated against two variables:
- E_rest: the antineutrino energy in the struck-nucleon rest frame (0.2–20 GeV, logarithmic);
- u = (Q² − Q²_min)/(Q²_max − Q²_min), at that energy.

Both sides use the same form factors:
- dipole vector form factors (M_V² = 0.71 GeV²), G_E^n = 0;
- dipole axial form factor with M_A = 0.99 GeV and g_A = 1.267;
- SU(3) couplings F = 0.463, D = 0.804;
- the same Cabibbo angle.

The correct side is the full contraction, evaluated numerically with Dirac matrices
(`make_tables.py`). A ratio of differential cross sections at fixed (E, Q²) corrects the
**kinematics as well as the rate**: in a two-body reaction, E and Q² fix the final state up to a
flat azimuth.

## Files

| file | content |
|---|---|
| `hyp_xsec_weights_mu.root` | tables `w_3122`, `w_3212`, `w_3112` (TH2D, log₁₀E_rest × u) for ν̄μ. Also: σ ratios `sigratio_*` (TH1D in log₁₀E) and the hyperon polarisation in production, `PL_*` (longitudinal) and `PT_*` (transverse, towards the ν̄), from the same matrix element |
| `hyp_xsec_weights_e.root` | the same, with the electron mass, for ν̄e |
| `GenieHypXsecWeight.h` | header-only lookup: loads both tables and picks one by neutrino flavour |
| `example_gst.C` | weights a GENIE `gst` file and writes a friend tree |
| `validate_closure.C` | the closure test below (MicroBooNE generator files) |
| `make_tables.py` | regenerates the tables (`python3 make_tables.py out.root [mu\|e]`); change the form factors or M_A here |

## Which events

Weight **only** GENIE `QEL-CC-LAMBDA` events:
- **probe:** an antineutrino (GENIE has this channel for ν̄ only);
- **process:** charged current, quasi-elastic;
- **final state:** a Λ (3122), Σ⁰ (3212) or Σ⁻ (3112) among the primary particles.

| format | selection |
|---|---|
| gst | `neu < 0 && cc && qel`, hyperon in `pdgf[]` (GENIE does not re-interact hyperons) |
| art `simb::MCTruth` | `GetNeutrino().CCNC() == 0`, `Mode() == 0` (QE), `Nu().PdgCode() < 0`, a primary (status 1, mother the vertex) hyperon. The struck nucleon is in `simb::GTruth::fHitNucP4` |
| GHEP | `ProcInfo().IsQuasiElastic() && IsWeakCC() && ExclTag().IsStrangeEvent()` |

Do **not** weight associated production (RES or DIS Λ K⁺, Σ K, ...). That comes from other
GENIE models and is not affected. Everything else gets w = 1.

## Usage

```cpp
#include "GenieHypXsecWeight.h"
GenieHypXsecWeight::Weights W;
W.load("hyp_xsec_weights_mu.root", "hyp_xsec_weights_e.root");

// per selected event: four-vectors in GeV
double w = W.weight(nuPdg, hyperonPdg, nubar, hitNucleon, lepton);
// no struck-nucleon record: E_rest = E_nu (nucleon at rest), Q2 = -(nubar - lepton)^2
double w0 = W.weight(nuPdg, hyperonPdg, Enu, Q2);
```

```bash
root -l -b -q 'example_gst.C("gntp.0.gst.root","hypw_friend.root")'   # gst->AddFriend("hypw","hypw_friend.root")
```

Use the struck nucleon as generated. Its Fermi motion shifts E_rest, and the tables are per free
nucleon. With the nucleon at rest the weight is approximate where the ratio varies quickly: near
threshold and near Q²_max.

## Validation

1. **Closure on generated samples** (`validate_closure.C`). We compared two MicroBooNE NuMI RHC
   productions on argon:
   - GENIE as distributed: 300 chunks, 1.24×10²³ POT;
   - GENIE with the corrected cross section: 443 chunks, 1.98×10²³ POT.

   ν̄μ CC Σ⁰ per 10²¹ POT:

   | sample | Σ⁰ per 10²¹ POT |
   |---|---|
   | GENIE as distributed | 193.0 |
   | **GENIE × this weight** | **128.6 ± 0.8** |
   | GENIE corrected | **127.4 ± 0.8** |

   Weighted over corrected, in bins of E_ν:

   | E_ν (GeV) | 0–1 | 1–2 | 2–3 | 3–4 | 4–5 | > 5 |
   |---|---|---|---|---|---|---|
   | weighted / corrected | 1.034 ± 0.017 | 0.998 ± 0.012 | 0.998 ± 0.025 | 1.08 ± 0.05 | 0.93 ± 0.08 | 0.97 ± 0.11 |

   An earlier version of the same comparison found agreement to 0.4% in the E_ν, Q², p_μ, cos θ_μ
   and p_Σ distributions (SigmaUBooNE `genie_hypfix/doc/`).
2. **GENIE `gst` end to end** (`example_gst.C`). `gevgen` v3_00_04, tune G18_10a_02_11a,
   2000 ν̄μ events on Ar40, generator list `LambdaCCQE`, flux ∝ E·exp(−E/1.2 GeV) over 0.5–5 GeV.
   Then `gntpc -f gst`. Mean weight per hyperon:

   | hyperon | events | mean weight |
   |---|---|---|
   | Λ | 832 | 1.216 |
   | Σ⁰ | 345 | 0.680 |
   | Σ⁻ | 823 | 0.676 |

   These match the GENIE spline ratios, corrected/Pais: Λ 1.27 / 1.22 and Σ 0.67 / 0.69 at
   1 / 2 GeV.
3. **Corrected cross section.** The trace calculation reproduces Llewellyn Smith CCQE (ν 1.06,
   ν̄ 0.41 ×10⁻³⁸ cm² at 1 GeV). It also reproduces a closed form, from the Gordon decomposition,
   to 3×10⁻¹² for every coupling combination. Splines made by GENIE with the corrected algorithm
   agree with these tables' σ ratios to 1–2%. The residual is GENIE's own G_E^n ≠ 0, which is not in
   the tables.

## Caveats

- **Capped weights:** weights are capped at 5. Near Q²_max GENIE's expression is negative and
  clipped to zero, so GENIE generates no events there. That region holds ≤ 0.16% (Λ) and ≤ 0.02%
  (Σ) of the correct cross section, and no weight can recover it.
- **Tune dependence:** the tables use M_A = 0.99 GeV with dipole vector form factors. The weight is
  a ratio, so the form-factor choice largely cancels: M_A from 0.961 to 0.99 GeV moves it by less
  than 0.3%. AR23_20i has M_A = 0.995 GeV, F, D and FA(0) as here. The vector form factors differ
  (GENIE G_E^n ≠ 0), which is about a 2% effect on Σ in the splines and is not in the table. If you
  need better than ~2%, regenerate the tables with your form factors.
- **Lepton-mass terms:** the corrected cross section includes them; GENIE omits them (0.1–0.7%).
  Electron and muon have their own tables.
- **Polarisation:** the weight does not polarise the hyperons. GENIE produces them unpolarised and
  Geant4 decays them isotropically. The predicted polarisation (`PL_*`, `PT_*`) is in the ROOT file
  if your analysis uses decay angles. A decay-angle weight for Σ⁰ → Λγ → pπ⁻ is `pelee/PolWeights.C`
  in the SigmaUBooNE repository.
- **No re-interaction:** GENIE does not re-interact hyperons in the nucleus. This weight does not
  change that.
- **Statistics:** weighting costs statistics where w is far from 1. For a large new production,
  generating with the corrected cross section is better: `genie_hypfix/` (GENIE v3_00_04), or the
  5-line patch in `GENIE_ISSUE.md`.
