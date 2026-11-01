# hypfix — correct ΔS = 1 quasi-elastic hyperon cross section for GENIE

> **Versioned copy.** The working copy runs from
> `~/MicroBooNE/production/nuwro_prod/genie_gen/hypfix/` (not under git); its scripts use paths
> relative to that location (`../../setup_env.sh`, `$SCRIPTDIR/hypfix/env_fix.sh`). Here the
> `generation/` directory holds the files that live one level up there
> (`gen_genie_chunk_hypfix.sh`, `submit_genie_hypfix.sh`, `genie_lambdaccqe_numi_{rhc,fhc}_hypfix.fcl`).
> Not copied: the build products (`lib/`, `test/port_check`), the test splines and logs, and the
> merged spline file `/data/uboone/genie_samples/hypfix_splines/gxspl-FNALsmall-hypfix.xml`
> (490 MB; `make_splines.sh` + `merge_splines.sh` regenerate it).

GENIE's `genie::PaisQELLambdaPXSec` (v3_00_04_ub4 as used by uboonecode v08_00_00_85, and
v3_00_06) gets the hyperon–nucleon mass-difference terms wrong. Against the full contraction
of the leptonic and hadronic tensors with the same SU(3) form factors, its ν̄ Λ cross section
is ×0.80 and its Σ⁰ and Σ⁻ ×1.5–1.85; Σ⁰/Λ comes out 0.45 instead of 0.19. NuWro's
`Singh_Model` agrees with the contraction (0.21 on Ar). With the final baryon at the nucleon
mass the Pais formula does reproduce Llewellyn Smith CCQE, which is why it looked fine.
The analysis-side evidence is in `Sigma/README.md`, "GENIE hyperon cross section".

`genie::HypQELTracePXSec` is a drop-in replacement: same channels (ν̄p→l⁺Λ, ν̄p→l⁺Σ⁰,
ν̄n→l⁺Σ⁻), same form factors (`LwlynSmithFFDeltaS`), same Cabibbo angle, same integrator,
same dσ/dQ² convention and nuclear scaling. Only `XSec()` differs: it evaluates
L_{μν}H^{μν} numerically from Dirac matrices (a C++ port of `Sigma/pelee/HypXsecWeights.py`).

GENIE is not rebuilt. The class lives in its own library, which ROOT autoloads through its
rootmap, and an XML override on `GXMLPATH` points the `QEL-CC-LAMBDA` generator at it.

## Files

| file | what |
|---|---|
| `HypQELTracePXSec.h/.cxx`, `LinkDef.h` | the algorithm |
| `build.sh` | builds `lib/libGHypFix.so` + dictionary + rootmap (run in the SL7 container) |
| `xml/ModelConfiguration.xml` | G18_10a copy with `XSecModel@…/QEL-CC-LAMBDA` → `genie::HypQELTracePXSec/Default` |
| `xml/master_config.xml` | GENIE's copy plus the entry for `HypQELTracePXSec.xml` |
| `xml/HypQELTracePXSec.xml` | algorithm config (as `PaisQELLambdaPXSec.xml`) |
| `env_fix.sh` | source after `setup_env.sh`: prepends `xml/` to `GXMLPATH`, `lib/` to `LD_LIBRARY_PATH` |
| `make_splines.sh` | `gmkspl` for all 44 targets of the production spline file × ν̄e, ν̄μ, ν̄τ, 250 knots to 1000 GeV |
| `merge_splines.sh` | production `gxspl-FNALsmall.xml` + the new splines → `gxspl-FNALsmall-hypfix.xml` |
| `test/port_check.cxx`, `test/run_port.sh` | C++ contraction against the Python trace |
| `test/spl_test.sh` | quick `gmkspl` on p, n, Ar40 |

Generation: `../genie_lambdaccqe_numi_{fhc,rhc}_hypfix.fcl` (the standalone fcls with
`XSecTable` pointing at the merged file) run by `../gen_genie_chunk_hypfix.sh` (sources
`env_fix.sh`). `gen_genie_chunk.sh` itself is untouched: running jobs read it from disk.

## Validation

1. **Port.** `test/run_port.sh` gives free-nucleon σ (10⁻⁴⁰ cm²) Λ 1.330 / Σ⁰ 0.252 / Σ⁻ 0.499
   at 1 GeV and 2.145 / 0.406 / 0.809 at 2 GeV, identical to the Python trace, which is itself
   checked against Llewellyn Smith CCQE (ν 1.06, ν̄ 0.41 ×10⁻³⁸ cm² at 1 GeV).
2. **Inside GENIE.** Splines made by GENIE with the new algorithm, divided by the production
   Pais splines, ν̄μ: Λ 1.268 / 1.218 / 1.189 at 1 / 2 / 5 GeV (expected from the weight
   tables 1.269 / 1.216 / 1.185); Σ⁰ 0.670 / 0.685 / 0.698 and Σ⁻ 0.665 / 0.682 / 0.695
   (expected 0.656–0.687; GENIE's own G_E^n ≠ 0 raises Σ by ~2%). Argon scales exactly as the
   free nucleon. With GENIE's form factors Σ⁰/Λ on a free proton is 0.24 at 1 GeV.

## What it does not change

Geant4 decays Σ⁰ → Λγ and Λ → Nπ isotropically, so the Σ⁰ polarisation still needs the
decay-angle weight (`Sigma/pelee/PolWeights.C`). GENIE still does not re-interact hyperons.
