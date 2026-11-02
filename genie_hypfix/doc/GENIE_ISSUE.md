**Title:** PaisQELLambdaPXSec: wrong hyperon–nucleon mass terms in the ν̄ ΔS=1 QE cross section (Λ ×0.8, Σ⁰/Σ⁻ ×1.5)

---

### Summary

`genie::PaisQELLambdaPXSec::XSec` computes the ν̄ ΔS = 1 quasi-elastic cross sections
(ν̄p → l⁺Λ, ν̄p → l⁺Σ⁰, ν̄n → l⁺Σ⁻). It is wrong whenever the final baryon is heavier than the
nucleon. Against the full contraction L_μν H^μν, computed with the same SU(3) form factors, GENIE
gives:

| | 1 GeV | 2 GeV |
|---|---|---|
| Λ | 0.79× | 0.82× |
| Σ⁰ | 1.52× | 1.49× |
| Σ⁻ | 1.52× | 1.49× |

Σ⁰/Λ comes out 0.45 instead of 0.24. NuWro's hyperon model (`Singh_Model`) agrees with the full
contraction.

**Affected:** v3_00_04 (incl. the MicroBooNE build v3_00_04_ub4), v3_00_06, v3_06_02, and `master`
(checked 2026-11-02). The lines are `src/Physics/Strange/XSection/PaisQELLambdaPXSec.cxx`
111–115 in v3_06_02.

### Issues

**1. The γ^μ coefficient assumes equal baryon masses (the main error, factor ~1.5).**
For on-shell spinors with M_Y ≠ M_N, the Gordon identity gives
F1 γ^μ + i(ξF2/2M) σ^{μν} q_ν = [F1 + ξF2 (M_N+M_Y)/(2M)] γ^μ − (ξF2/2M)(p_N + p_Y)^μ.
The code uses `(F1V + xiF2V)` for the γ^μ coefficient, which holds only for M_Y = M_N = M. It
enters `w1`, `w2` (through `(F1V + xiF2V - Mp*xiF2V/(2*Mnuc))`, which should reduce to F1V) and
`w3`.

**2. A missing factor M² in the `w2` term of `xsec` (1–2%).**
`-4*(Mnuc2*v2 - q2)*w2` subtracts q² (GeV²) from M²ν² (GeV⁴). It should be
`-4*Mnuc2*(v2 - q2)*w2`. This term is also wrong in the equal-mass limit, by 1% at 1 GeV, which
is why the formula appears to reproduce Llewellyn Smith CCQE when M_Y = M_N.

**3. A negative dσ/dQ², hidden by clipping.**
Because of 1 and 2, the expression turns negative near Q²_max, and `xsec = TMath::Max(xsec, 0.)`
sets it to zero:

| | 1 GeV | 3 GeV |
|---|---|---|
| Λ | u > 0.90 | positive |
| Σ⁰ | u > 0.965 | u > 0.98 |
| Σ⁻ | u > 0.965 | u > 0.98 |

with u = (Q² − Q²_min)/(Q²_max − Q²_min). No events are generated there. The region holds 0.16%
(Λ) and 0.02% (Σ) of the correct cross section. Once 1 and 2 are fixed, the expression is positive
everywhere.

**4. Lepton-mass terms are dropped (0.1–0.7%; minor).**
The `w`-terms omit the m_l² contributions of the contraction, and `Fp` is commented out. After
fixes 1 and 2 these are the only difference from the exact result: 0.7% for Σ⁰ at 1 GeV, ≤ 0.2%
above 2 GeV.

### Proposed patch (fixes 1 and 2)

```cpp
double FM = F1V + Mp * xiF2V / (2 * Mnuc);        // gamma^mu coefficient for M_Y != M_N
double w1 = (Mm2 - q2)/(4*Mnuc2)*TMath::Power(FM, 2) + (Mp2 - q2)/(4*Mnuc2) * FA2;
double w2 = FA2 + TMath::Power(F1V, 2) - q2 / Mnuc2 * TMath::Power((xiF2V / 2), 2);
double w3 = 2 * FA * FM;
double xsec = Gf*fSin8c2 / (16*Mnuc2*E2) * (-8*Mnuc2*q2*w1 - 4*Mnuc2*(v2 - q2)*w2 - sign*2*(s - u)*q2*w3 + (s-u)*(s-u)*w2);
```

(`FM` must be defined after `Mp`.) With a massless lepton this equals the exact contraction
identically (symbolic check). We also have a drop-in algorithm, `HypQELTracePXSec`, which
evaluates the full contraction including the lepton-mass terms. It has the same form factors,
Cabibbo angle, integrator and nuclear scaling, and we can contribute it if preferred.

Effect on σ, relative to the exact contraction:

| | as is | fix 1 | fixes 1+2 |
|---|---|---|---|
| Σ⁰, 1 / 2 GeV | 1.523 / 1.489 | 0.990 / 0.998 | 1.007 / 1.002 |
| Λ, 1 / 2 GeV | 0.788 / 0.823 | 0.978 / 0.995 | 1.001 / 1.001 |

### Validation

- **Independent contraction.** L_μν H^μν is evaluated from explicit Dirac matrices. It reproduces
  Llewellyn Smith CCQE (ν 1.06, ν̄ 0.41 ×10⁻³⁸ cm² at 1 GeV). It also reproduces a closed form
  derived from the Gordon decomposition to 3×10⁻¹² for every coupling combination (f₁, f₂, g₁ and
  their interferences).
- **Term by term.** At fixed Q², the coefficients of (s−u)⁰, (s−u)¹ and (s−u)² isolate the
  errors:
  - the constant term is wrong for every coupling;
  - the (s−u)² term is wrong only with f₂;
  - the interference term is wrong only with f₂·g₁.
- **Splines.** We made splines with GENIE (tune G18_10a_02_11a) using the corrected
  cross section. New/Pais for ν̄μ on argon, at 1 / 2 / 5 GeV:
  - Λ: 1.271 / 1.219 / 1.189;
  - Σ⁰: 0.672 / 0.686 / 0.698.
- **Generation.** NuMI RHC on argon, 1.98×10²³ POT corrected against 1.25×10²⁴ POT Pais. Primary
  hyperons per POT, corrected/Pais:
  - Λ 1.257, Σ⁰ 0.663, Σ⁻ 0.670.

  Reweighting the old events by the ratio of the differential cross sections reproduces the new
  generation to 0.4% in rate and in every kinematic distribution.

### Impact

Every Λ, Σ⁰ and Σ⁻ prediction from this channel is affected: Λ is low by ~20%, Σ⁰ and Σ⁻ are high
by ~50%, and Σ⁰/Λ is ×1.9 too large. Published comparisons of hyperon production with GENIE, for
example MicroBooNE's Λ measurement (arXiv:2212.07888), used the uncorrected model.

### Not in scope (noted for completeness)

- The produced hyperons are unpolarised. The polarisation predicted by the same matrix element is
  not propagated to the decays.
- INTRANUKE does not re-interact hyperons.

Both are separate from this bug.

Documentation, figures and derivation scripts are available on request (or attached).
