# MicroBooNE Sigma0 analysis

Three selections live in this repo. They are fully independent — different ntuple
formats, different samples, separate inputs and outputs — and can be run side by side.

| chain | format | signal MC | S | B | purity | role |
|---|---|---|---|---|---|---|
| `SigmaSelection.C` + BDT | PeLEE | NuWro run3b | 2.85 | 841 | 0.34% | **the analysis**: the only chain that can see data, EXT and dirt |
| `DumpBDTHyp.C` + BDT | hypana | GENIE | 19.3 | 670 | 2.8% | truth studies; candidate purity, backtracking (Run 3 background only) |
| `SigmaLambdaSel.C` | hypana | GENIE | 1.33 | 3.76 | **26%** | **best performing**; needs wire activity, so MC-only for now (NuWro: S 0.142) |

NuWro yields include hyperon FSI and GENIE yields do not, a model difference of roughly
4x in the Σ⁰ rate (see "The three signal samples disagree on the absolute Σ⁰ rate"),
so compare the chains by purity and efficiency, not by S across rows.

All yields are POT-scaled to the full RHC data exposure 1.1526e21 (see
[POT normalisation](#pot-normalisation)). The three are not alternatives to choose
between yet: the Lambda-style chain wins decisively on purity and on robustness to
background systematics, but until data/EXT/dirt exist in hypana format it cannot
produce a result on its own.

## Environment

Anything that loads TMVA from an **ACLiC-compiled** macro (i.e. `SigmaSelection.C+`
once BDT weights exist) needs a shim on this machine: `/etc/flexiblasrc.d/netlib.conf`
registers a BLAS backend whose library is not installed, and cling aborts with exit 1
and no message beyond `libflexiblas_netlib.so: cannot open shared object file`.
Interpreted TMVA (`TrainBDT.C`) is unaffected, which makes it look like a code bug.

```bash
mkdir -p ~/lib/flexiblas-shim
ln -sf /usr/lib64/flexiblas/libflexiblas_openblas-openmp.so \
       ~/lib/flexiblas-shim/libflexiblas_netlib.so
export LD_LIBRARY_PATH=$HOME/lib/flexiblas-shim:$LD_LIBRARY_PATH
```

Setting `FLEXIBLAS=` alone does not fix it. Two other traps: a stale
`SigmaSelection_C.so` causes a separate silent exit-1 (delete it and use `.C++` to
force a clean rebuild), and macro output is block-buffered, so a crashed run leaves an
**empty** log — absence of output says nothing about where it died.

## Main analysis — PeLEE ntuples

Sigma0 selection on standard PeLEE ntuples (tree
`nuselection/NeutrinoSelectionFilter`), NuWro hyperon-enhanced signal + all six RHC
overlay runs (Run 1+2+3+4a+4b+4c) + RHC dirt + run3b NuMI EXT + NuWro fake-data, with a gradient-boosted BDT
as the final discriminant under two-fold cross-validation (folds trained on even/odd
event numbers, each applied to the events it did not train on). All work is in
`pelee/`; run from there:

```bash
root -l -b -q 'SigmaSelection.C++'  # 1st pass -> SigmaHists.root + BDTInput.root (BDTG passes through)
root -l -b -q TrainBDT.C            # train 2 folds -> SigmaBDT_A/, SigmaBDT_B/
root -l -b -q 'SigmaSelection.C+'   # 2nd pass -> applies the cross-validated BDTG
root -l -b -q 'PlotSigma.C+'        # plots -> figs/stack/, figs/shape/
```

The BDT bootstraps itself: with no weights present the BDTG stage passes through and
dumps the training tree; after training, re-running applies the cut.

Source: `SigmaSelection.C`, `SigmaSelConfig.h`, `SigmaFV.h`, `PidAssign.h`,
`TrainBDT.C`, `PlotSigma.C`. Samples/POT are an editable block at the top of
`SigmaSelection.C` (default exposure = full RHC data POT 1.1526e21; background overlay
POT 1.5463e22). `ScanPeleeBDT.C` relocates the operating point from `BDTInput.root`.
Detector systematics: `./run_detvar.sh` (see "Detector variations" below).

### Key elements

- **Cut chain** (9 stages): FiducialV, CosmicVeto (`topological_score`>0.20),
  NShowers, MuonCandidate, ProtonCandidate, PionCandidate, NTracksCut, BDTG.
- **The proton and pion candidates are OPTIONAL** (`requireProtonPion=false`). They
  were never PID cuts — they demanded ≥1 and ≥2 non-muon tracks — and because roughly
  half the signal has no reconstructed pion they removed more signal than background
  (−31.6%/−18.3% and −52.9%/−43.5%), more than halving the figure of merit. Events
  without a pair are kept, every pair-dependent BDT input carries a sentinel, and the
  input `HasPair` flags them. Those two cut-flow rows therefore repeat MuonCandidate.
  Set `requireProtonPion=true` to restore the old chain.
- **Joint proton/pion assignment** (`PidAssign.h`): each non-muon track is scored
  separately for the two roles with a quadratic logistic model and the best ordered
  pair is taken. Against HyperonProduction truth this raises pion-candidate purity
  from 33.5% to 42.5% at no cost in efficiency. The old rule ranked tracks by a single
  monotonic discriminant and took the extremes, which cannot isolate a pion whose
  ionisation sits in a band between minimum-ionising and proton-like.
- **Sigma0 mass** reconstructed from the per-PFP shower energy `shr_energy_y_v` (the
  scalar shower-energy branch is mostly sentinel).
- **57 BDT inputs**: track PID/kinematics, inter-candidate distances/angles, Lambda &
  Sigma0 masses, shower energy/size, Lambda/gamma vertex-pointing, pi0 (di-photon)
  variables, the second-shower block, and `HasPair`.
- **Operating point** `BDT_MIN = 0.616`, the out-of-sample FoM peak of the 2026-10-15
  training (before EXT existed; with EXT in B the peak moves to 0.690 for +1.2% FoM, flat).

### Results

Scaled to RHC data POT, with the unbiased two-fold cross-validation (2026-10-17, all
samples below; B = hyperons + overlay + dirt + EXT):

| Stage | Raw signal | Eff [%] | S | B | Purity [%] | S/sqrt(S+B) |
|---|---|---|---|---|---|---|
| MuonCandidate | 11107 | 51.61 | 6.469 | 28978.3 | 0.022 | 0.0380 |
| NTracksCut | 10931 | 50.79 | 6.370 | 26548.8 | 0.024 | 0.0391 |
| **BDTG (final)** | **4904** | **22.79** | **2.850** | **840.7** | **0.338** | **0.0981** |

| B at the BDTG point | before 2026-10-17 | now | raw events now |
|---|---|---|---|
| other hyperons (NuWro) | 1.74 | 1.74 | 2982 |
| overlay | 730.70 | 715.71 | 9966 |
| dirt | 17.48 (FHC!) | 13.83 | 37 |
| EXT | — | **109.46** | 146 |
| total | 749.9 | 840.7 | |

EXT is 13% of the final background; at MuonCandidate it is 6,844 against 20,348 from
the overlay. Signal is unchanged (only background samples changed).

**Caveat on the operating point.** At the peak a 2% systematic on the background is
±17 events against 2.9 signal, and purity never exceeds ~0.5% anywhere on this
curve. `S/sqrt(S+B)` cannot see that, so the systematic budget — not the counting
metric — should ultimately set `BDT_MIN`. Rescan with `ScanPeleeBDT.C` after any
change to the training or the preselection.

### Sample changes of 2026-10-17

- **The dirt was FHC.** `/data/uboone/dirt/neutrinoselection_filt_run1_dirt_overlay.root`
  (4.21e20, called "Run 1 dirt" here) has the true-E_ν spectra of the FHC dirt samples
  (χ²/ndf ≈ 1) and not of the RHC ones (χ²/ndf 90–660): ⟨E⟩ ν_μ 1.60 GeV against 2.15
  for RHC, ν̄_μ 2.28 against 1.72. It is replaced by the RHC dirt of runs 1, 2b, 3b, 4a and
  4b from `/data/uboone/temp/custom_pelee_ntuples/dirt_overlay/` (1.7409e21 POT). Equal
  run/subrun/event numbers between dirt samples do NOT mean equal events: they share the
  overlaid EXT (5 of 139,525 matched events have the same neutrino).
- **The RHC dirt is on the OLD NuMI flux**, so it takes the same flux map as NuWro. Its
  `ppfx_cv` profile in true E_ν matches the old-flux overlay in every bin (ν̄_μ 0.886–0.910
  from 0 to 5 GeV) and not the current-flux overlay (0.904–1.012).
- **Run 3 overlay was missing**: `Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa..ae`, 5.5185e21
  POT, +55% background statistics.
- **EXT added** (sample 5): `/data/uboone/EXT/neutrinoselection_filt_run3b_beamoff.root`,
  32,649,128.65 gates (analyser's table), scaled by 2.4477e7 beam-on triggers (the table's
  23,534,937 per 1.10823e21, taken at 1.1526e21): 0.7497. Data files have **empty**
  `backtracked_*` vectors whenever there are PFPs, so the ragged-vector guard used to drop
  every neutrino-slice EXT event; it no longer requires them for data. The 1,655 repeated
  run/subrun/event in the file (1,619 with identical reco) are dropped.
- **Overlay signal removed from the background**: truth Σ⁰ CC events in the true FV inside the
  inclusive overlay (1,038 raw) were counted in B while the NuWro sample models the signal;
  they are now dropped, as the hypana chain does.

### Detector variations

The run 3 RHC detector-variation set (`custom_pelee_ntuples/detvar/run3`: CV, 2.48e21 POT, and
LYAttenuation, LYDown, LYRayleigh, Recombination, SCE, WireMod X / YZ / θXZ / θYZ) re-simulates
the same generated events, but some variations lost files (θYZ keeps 81%). `run_detvar.sh`
runs `SigmaSelection.C` on each (its `detvarFile` argument: one file, truth signal in the true
FV vs everything else, a per-event tree `dv`) and `DetVar.C` compares every variation with the
CV on the events they share, with the CV weight on both sides. These samples are on the old
flux, so only the ratios are used.

| variation | S var/CV at BDTG | B var/CV at BDTG |
|---|---|---|
| LYAttenuation | 1.000 ± 0.000 | 1.001 ± 0.004 |
| LYDown | 0.979 ± 0.021 | 0.993 ± 0.004 |
| LYRayleigh | 0.997 ± 0.027 | 0.985 ± 0.006 |
| Recombination | 1.001 ± 0.076 | 0.975 ± 0.021 |
| SCE | 0.851 ± 0.099 | 0.999 ± 0.025 |
| WireMod θXZ | 0.881 ± 0.048 | 1.032 ± 0.018 |
| WireMod θYZ | 0.865 ± 0.064 | 1.001 ± 0.017 |
| WireMod X | 0.959 ± 0.057 | 1.004 ± 0.017 |
| WireMod YZ | 0.958 ± 0.070 | 0.995 ± 0.017 |
| **quadrature** | **24.2%** (MC noise alone 17.7%) | **4.5%** (MC noise alone 4.8%) |

Background: no variation is resolved; the quadrature sum equals its MC noise, so the
background detector systematic at the BDT point is ≲5%. That alone is ±38 events against
S = 2.85. Signal: 64 raw CV signal events at the BDT point, so most of the 24% is MC noise;
the angular wire modifications (θXZ 2.5σ, θYZ 2.1σ) are the only resolved effects. Earlier
stages (MuonCandidate) give 7.8% signal / 2.3% background.

## Truth-study chain — HyperonProduction ("hypana") ntuples

The PeLEE truth record stops at the primary vertex and never contains the
Sigma0 → Lambda gamma → p pi- gamma chain, so it cannot say what a candidate really
is. HyperonProduction backtracks every PFP, which is what makes candidate purity,
the pion candidate and the photon measurable at all. `HypInput.h` reads those ntuples
through the same variable names the PeLEE selection uses.

```bash
root -l -b -q HypCutFlow.C                   # POT-weighted cut flow + candidate purity
root -l -b -q 'DumpBDTHyp.C("BDTHyp.root")'  # training tree from lists/
root -l -b -q AblateBDT.C                    # 2-fold training + variable-block ablation
root -l -b -q CurveBDT.C                     # efficiency curve -> bdt_efficiency_curve.png
root -l -b -q PeakBDT.C                      # locate the FoM peak
root -l -b -q BkgDisplay.C                   # event displays of surviving background
```

File lists are repo-tracked in `pelee/lists/`; regenerate by globbing `*_hypana.root`
under the GENIE run3b signal and RHC overlay hypana directories. Operating point in
`HypConfig.h` (`BDTG_CUT = 0.554`).

**Neutrino-slice requirement.** PeLEE applies one implicitly — 58.3% of its overlay
events have `nslice==0`, no PFPs and no usable vertex, so they die at the fiducial
cut. hypana attaches a Pandora vertex to those events regardless, so they survived and
inflated the background 2.4x. The matching quantity is `flash_match_nu_slice_ID >= 0`
(missing in 57.0% of hypana events), added as an explicit `SliceID` stage. With it the
two chains agree to 6-8% throughout.

### What the truth says

Candidate purity on signal at preselection: muon 93.5%, photon 72.1%, proton 66.4%,
pion 46.7%. Only **51.4%** of selected signal events contain a reconstructed true pi-
at all — a reconstruction limit no selection or classifier work can reach, and the
binding constraint on any pion-based approach.

The Lambda mass separates far less than its ranking suggests, and not for the reason
first assumed. HyperonProduction's `trk_energy_proton` is the proton **momentum**, not
the kinetic energy the PeLEE branch of the same name carries; `HypInput.h` and
`DumpBDTHyp.C` now convert it. With the proton right, the signal peak is centred
(median 1.111 GeV) and 54.6% of paired signal sits within 20 MeV of the true mass --
but so does 34.9% of paired background, up from 18.8%. The paired background is
largely genuine p + pi from Delta production at the primary vertex, whose invariant
mass with correct momenta sits just above the p+pi threshold, exactly where the
Lambda is. Fixing the reconstruction sharpened both and left the classifier
unchanged (AUC 0.828). The **main PeLEE analysis was never affected**: its branch is
genuinely kinetic energy. `BkgDisplay.C` shows this one event at a time — of ten background
events that fake a Lambda, none contained a true pi-.

### Measured dead ends

- **More signal MC does not help.** The out-of-sample AUC did not move when the signal
  sample grew 25%.
- **Lambda displacement variables do not work** at this resolution: 1.85 cm flight
  against 5.54 cm resolution, truth correlation r = −0.014.
- **Better candidate assignment helps background almost equally**, because the
  background's Lambda mass is combinatoric and benefits from the same improvement.
- **Correcting the proton momentum helps background almost equally** too, for the
  same reason (Delta -> p pi sits at the Lambda mass). Third instance of the pattern.
  It also barely matters for the mass itself: the Lambda -> p pi Q-value is 38 MeV,
  so the proton is nearly at rest in the Lambda frame and the invariant mass is
  dominated by M_p plus the pion energy. That is why a 2.3x proton error shifted the
  mass only ~20 MeV.
- **Correcting the pion momentum makes the mass WORSE.** The pion candidate uses the
  muon-hypothesis range momentum, which is 20% low for a pion. Removing that bias
  (power-law fit, or a 1.25 scale) roughly doubles the Lambda mass width (7.8 -> 12.5
  or 16.0 MeV half-IQR) because the pion carries most of the decay kinetic energy and
  scaling it up amplifies its spread. A sharp peak 10 MeV low beats a centred broad
  one for a classifier. Keep the muon-hypothesis estimator.
- **The seeded island test does not help the BDT**, though the test itself works (see
  the Lambda-style chain below). It reaches only ~8% of the BDT's sample, because 70%
  of events have no proton/pion pair to seed it with and the pair is optional there.
  It needs wire activity, which only the hypana ntuples carry. *Superseded:* an
  earlier version of this note said the test was "capped at 27.2%" by short Lambda
  flight. That was a tuning artefact — see `IslandTune.C`.

## Lambda-style selection — `SigmaLambdaSel.C` (hypana)

The highest-purity chain in the repo, and the one whose sensitivity survives
background systematics. It follows the Lambda CCQE analysis (Thorpe & Nowak,
*Hyperon\_Technote\_v2*) rather than the topology+PID chain of `SigmaSelection.C`,
because at equal efficiency that analysis reaches 91% purity against our 2.8%.

```bash
root -l -b -q 'SigmaLambdaSel.C("lists/sig_genie.list",0,true)'    # signal, raw counts x 1.9284e-2
root -l -b -q 'SigmaLambdaSel.C("lists/bkg_overlay.list",0,false)' # background, raw counts x 1.6352e-1
```

Source: `SigmaLambdaSel.C`, `LambdaStyle.h` (range momenta, V-fit secondary vertex,
pair features), `Island.h` (connectedness test), `SigmaFV.h`, `HypInput.h`.

### What it does differently

- **Every ordered pair is scored**, not each track separately. `PidAssign.h` scores a
  track for a role using log10(distance to the *primary* vertex) with coefficient
  −0.2758 on the pion — it penalises displaced tracks, i.e. exactly the signature of a
  Lambda that flew. The pair selector's only topological input is the separation of the
  two track *starts*, which is displacement-agnostic.
- **The precuts are a preference, not a veto.** Measured on the full samples they are
  dead weight as a veto: the same FoM with ~20% less signal. The best precut-passing
  pair is used when one exists, otherwise the best pair of all (`usePrecuts=true`
  restores the veto).
- **The decay vertex is fitted** (a V fit, validated to recover 3.78 cm for a true
  4.00 cm flight), giving a decay length and the pointing angle α. The old
  closest-approach displacement variables measured uncorrelated with truth (r = −0.014).
- **The connectedness island test is a cut**, applied last.

### Operating point

`LAM_MASS_WIN = 0.020`, `COS_ALPHA_MIN = 0.95`, decay length > 0.5 cm, **photon start
> 2 cm from the vertex** (`photonGap`, default on since 2026-10-17), island on
**≥ 2 planes**, no precut veto. Chosen from a full-statistics scan (`IslandAna.C`) and
checked out-of-sample by event parity (`IslandXVal.C`), by the RHC extension
(`IslandScan2()`, `IslandHandles.C`) and, for the gap, on the extension alone. The gap
is explained in [Two optional handles](#two-optional-handles-applymichel-photongap).

Full-statistics background (`lists/bkg_overlay.list`, 14,164 files, 7.0488e21 POT,
weight 0.1635), GENIE signal. These are `SigmaLambdaSel.C` on the full lists (logs
`sls_gap_sig.log`, `sls_gap_bkg.log`), and every row agrees with the `IslandOpPoint.C`
cache to the event:

| stage | raw sig | raw bkg | S | B |
|---|---|---|---|---|
| AllEvents | 15043 | 1840942 | 290.1 | 301031 |
| SliceID | 7801 | 776242 | 150.4 | 126931 |
| MuonCandidate | 3043 | 123000 | 58.7 | 20113 |
| PairCandidate | 921 | 51238 | 17.8 | 8378 |
| LambdaMass ±20 MeV | 504 | 16391 | 9.72 | 2680 |
| PointingAlpha | 144 | 706 | 2.78 | 115.4 |
| PhotonGap > 2 cm | 139 | 549 | 2.68 | 89.8 |
| **Island (≥2 planes)** | **69** | **23** | **1.33 ± 0.16** | **3.76 ± 0.78** |

Purity 26%, absolute efficiency 0.46%, and **98.6% of surviving signal has true Lambda
flight > 2 cm** (NuWro 93.9%) against 43% in the BDT chain. The displaced-vertex bias
is cured, which was the point of the rewrite. Errors are MC statistics, √Σw². Without
the gap the same point is S = 1.39 ± 0.16, B = 6.21 ± 1.01 (72/38 raw), purity 18%.

The **one-plane** point (`islandPlanes=1`) gives S = 1.89 ± 0.19, B = 14.9 ± 1.6 on
98/91 raw events (purity 11.3%); without the gap 1.99 and 21.4 (103/131).

Odd/even cross-validation (`IslandXVal.C`, f = 15%, gap on) has both folds choosing
≥ 2 planes with no precut veto. Out of sample: S 1.16, B 4.7 (60/29 raw), purity 20%,
FoM 0.457, against 0.380 without the gap.

**The background extension raised the two-plane B by 1.5x (numbers without the gap),
and the old value was the low one.** On the Run 3 list alone (the first 3,491 files, 2.7795e21 POT) the same code
gives B = 4.15 ± 1.31 on 10 raw events, the number this table used to quote. Every stage
up to PointingAlpha agrees per POT between the two samples to 1–4%; only the last step
differs. Given the 38 events, the Run 3 list holds 10 against 15.0 expected from its
POT share (1.6σ low); at one plane it holds 47 against 51.6. The connectivity tuning and
the choice of two planes were made on those 10 events, so the low value is partly the
tuning sample's own fluctuation. The extension alone, which played no part in the tuning,
gives B = 7.6 ± 1.4 (28 raw). The island-pass rate after PointingAlpha does not depend on
the era within the statistics (two-plane: Run 3 list 3.7 ± 1.2%, ext Run 1 3.3 ± 1.9%,
Run 2 8.3 ± 2.1%, Run 3 6.4 ± 1.9%; χ² ≈ 6.6 for 4 dof). Across the 123,000
MuonCandidate rows no (run, subrun, event) repeats, so the extension is statistically
independent of the original sample. Before the extension, one cross-validation fold
chose the precut veto; with full statistics neither does.

**With the NuWro signal** (`isl_sig_nuwro.root`, flux-weighted, see
[The NuWro signal sample](#the-nuwro-signal-sample-run3b-2026-10-15)) the same point gives
S = 0.142 ± 0.009 on 245 raw events (0.217 ± 0.015, 358 raw, at one plane), 9.4x below
GENIE. The selection is the same for both generators. NuWro's efficiency from true Σ⁰ to
the two-plane point is 0.42% against 0.46%, the gap keeps 97% of it against 96%, and
half of the PointingAlpha survivors pass the island test in both. The factor is the normalisation: hyperon FSI, which NuWro has and
GENIE does not, plus the flux weight.

### Why S/√(S+B) is the wrong metric here

With a fractional systematic *f* on the background, using S/√(S+B+(fB)²):

| f | BDT chain (S 19.3, B 670) | **Lambda-style (S 1.33, B 3.76)** | no gap (S 1.39, B 6.21) | Run 3 only, no gap (B 4.15) |
|---|---|---|---|---|
| 0% | 0.736 | **0.590** | 0.504 | 0.593 |
| 7% | 0.359 | **0.586** | 0.498 | 0.590 |
| 15% | 0.186 | **0.572** | 0.477 | 0.571 |
| 30% | 0.095 | **0.527** | 0.417 | 0.508 |

The crossover is at ~7% (the BDT column predates the background extension and the
NuWro retrain), and the technote quotes 30–60% on RES/DIS backgrounds. At any
realistic systematic the high-purity chain wins, and the full-statistics background
does not change that. Every operating point in this repo was
historically chosen with f = 0, which systematically favoured the BDT chains.

### Island finding — `Island.h`

Build a wire × time grid of filtered activity, seed clusters at the muon / proton /
pion track starts converted to (wire, tick) per plane, grow them, and require the
**proton and pion clusters to merge while the muon stays separate**. Per-plane, so a
bad view can be fallen back on. Coordinate transforms were derived and validated in
`ctcalib.C` (none existed): `wire = (A·y + B·z + C)/0.3` reproduces 98.7% of stored
vertex wires within 2 wires; `tick = 18.1509·x + 818.1`.

**Connectivity is tuned** (`CDW = 2` wires, `CDT = 20` ticks; was 3 × 40). The original
40 ticks is 2.2 cm of drift, so vertex activity bridged the decay V to the muon
whenever the Lambda flew less than ~3 cm. Seeds were never the problem — they land on
activity 98% of the time. Retuning lifted the signal pass rate at 2–3 cm flight from
45% to 80%, i.e. into the "highly efficient above 2 cm" regime the technote describes:

| true flight | 3×40 (old) | 2×20 (now) |
|---|---|---|
| 1–2 cm | 19% | 43% |
| 2–3 cm | 45% | 80% |
| 3–5 cm | 70% | 85% |
| 5–10 cm | 86% | 91% |

On **one** plane the FoM is a flat plateau across every connectivity — looser settings
admit background as fast as signal. The gain is on **two** planes (6.8× rejection
against 3.7×). Minimum cluster size, seed radius and an ADC occupancy threshold all
only remove signal. Do not trust a grid optimiser here: its own pick (cdw 6, minclust
20, seed radius 5) fails the parity check (0.393 vs 0.264) and is overfit to 11 events.

**2×20 survives an out-of-sample test (2026-10-17).** The tuning used the Run 3 list
only, so the RHC extension was cached separately (`islcache_bkg_ext.root`, 2,261
pre-island events against 2,298 expected from the Run 3 rate) and the connectivity grid
was re-scanned at the two-plane point, choosing on the Run 3 list and scoring on the
extension (`IslandScan2()` in `IslandTune.C`, log `islscan2.log`):

| connectivity | S (raw) | B Run 3 (raw) | B extension (raw) | FoM₁₅ Run 3 / ext |
|---|---|---|---|---|
| 1 × 20 | 1.41 (73) | 4.56 (11) | 8.10 (30) | 0.555 / 0.425 |
| **2 × 20** | **1.39 (72)** | **4.15 (10)** | **7.56 (28)** | **0.571 / 0.434** |
| 2 × 30 | 1.29 (67) | 4.56 (11) | 6.48 (24) | 0.514 / 0.438 |
| 3 × 30 | 1.23 (64) | 4.98 (12) | 5.94 (22) | 0.474 / 0.437 |
| 3 × 40 (old) | 1.04 (54) | 3.32 (8) | 5.67 (21) | 0.485 / 0.382 |

The extension's FoM is a plateau (0.425–0.438 over 1–3 wires × 20–30 ticks) and 2×20
is 0.004 from its top, so the choice was not overfit. The drop from 0.571 to 0.434 is
the same at every grid point: the Run 3 list is low across the whole grid because its
events are shared between neighbouring settings. It is a property of that sample, not
of the parameter choice. The one-plane scan (`islscan1.log`) is flat on the extension
as well (0.33–0.36). The two-plane signal pass rate is the same for both generators:
50.5% of NuWro (`islcache_sig_nuwro.root`, flux-weighted) against 50% of GENIE.

### Tooling

The island test reads raw wire windows, which is what makes hypana loops slow. Two
caches make everything downstream cheap:

```bash
root -l -b -q 'IslandProcess.C("lists/sig_genie.list",0,true,"isl_sig_genie.root",1.9284e-2)'
root -l -b -q 'IslandProcess.C("lists/bkg_overlay.list",0,false,"isl_bkg.root",1.6352e-1)'
root -l -b -q 'IslandProcess.C("lists/sig_nuwro.list",0,true,"isl_sig_nuwro.root",8.6055e-4,"<flux_rw>/nuwro_rhc_run3b_reco2_fluxweights.tsv")'
root -l -b -q 'IslandCache.C("lists/sig_genie.list",0,true,"islcache_sig.root",1.9284e-2)'
root -l -b -q 'IslandOpPoint.C("isl_sig_genie.root","isl_bkg.root")'   # S, B +- MC stat at both points
```

The pass costs ~1 ms/event. The full background (1.84M events) takes ~15 min when the
list is split 8 ways (`split -n l/8`) and the outputs are `hadd`ed. For NuWro the last
argument is the per-event flux-weight table (`FluxWeights.h`), and `weight` becomes POT
scale × flux weight. Current caches (2026-10-17, tuned 2×20 island): `isl_sig_genie.root`,
`isl_sig_nuwro.root`, `isl_bkg.root` (full list), `isl_bkg_run3.root` (Run 3 list only,
weight 0.4147). The September caches were made with the old 3×40 connectivity and are
kept as `isl_*_3x40.root`. `PhotonMichel.C` numbers below were measured on those.
`islcache_*.root` holds the raw wire windows used by `IslandTune.C`: `islcache_sig.root`
(GENIE), `islcache_sig_nuwro.root` (flux-weighted), `islcache_bkg.root` (Run 3 list,
weight 0.4147) and `islcache_bkg_ext.root` (extension, files 3492–14164, weight 0.26998).
Each is weighted to data on its own POT. `IslandScan2()` combines them, and
`IslandHandles.C` evaluates the Michel and gap handles on the `isl_*` caches.

- `IslandProcess.C` → one row per event at MuonCandidate: per-plane island result, V-fit
  observables, Σ mass, photon and Michel variables, true candidate PDGs, true flight.
- `IslandCache.C` → the raw wire windows plus seeds for the pre-island superset, so
  island *parameters* can be re-scanned without touching the files again.

Consumers: `IslandAna.C` (cut flow, placement, flight dependence, background
composition, threshold scans), `IslandTune.C` (parameter grid; `IslandPoint()` for one
set with parity split), `IslandXVal.C` (odd/even cross-validation of the operating
point), `IslandSideband.C` (sideband definitions), `PhotonMichel.C` (the two handles
below).

### Sidebands

| region | S | B | raw B | S/B |
|---|---|---|---|---|
| signal region (blind) | 1.89 | 14.9 | 91 | 12.7% |
| mass sideband, +40 to +200 MeV | 0.02 | 10.5 | 64 | 0.2% |
| α sideband, 0 < cos α < 0.8 | 0.50 | 58.2 | 356 | 0.9% |
| island-fail region | 0.79 | 74.9 | 458 | 1.1% |

(one-plane point for statistics; full background, tuned island and photon gap in every
region, 2026-10-17. `IslandSideband.C(...,false)` drops the gap.) All three are RES/DIS-dominated like the signal
region. The α sideband is the workhorse. Two things to know: the mass sideband exists
only *above* the peak — the p+π invariant mass has a kinematic floor at 1.078 GeV — and
the signal peak sits in the −20…0 MeV bin because the pion runs on the muon range
hypothesis.

### Two optional handles (`applyMichel`, `photonGap`)

`photonGap` is **on by default since 2026-10-17**, part of the operating point above;
`applyMichel` stays off.


- **Michel tag.** The signal muon is a μ⁺ 100% of the time; the surviving background
  muon is a true **μ⁻ 78%** (wrong-sign ν_μ CC in RHC). A μ⁺ always decays, a μ⁻
  captures 75% of the time in argon. Requiring a contained muon to have another PFP
  starting within 5 cm of its end keeps 78% of signal and 48% of background, and what
  it removes is **87% true μ⁻**. Validated in background: contained true μ⁺ tagged 25%
  vs true μ⁻ 16%.
- **Photon conversion gap.** 18.7% of background "photons" start within 2 cm of the
  vertex against 3.5% of signal; requiring > 2 cm keeps 95% of signal and removes 21%
  of background.

Both were measured on the Run 3 list in September (one-plane point, old 3×40 island)
and looked unmeasurable on the two-plane point, where that list has 8–10 MC events.
**With the full background the photon gap is a real gain at two planes**
(`IslandHandles.C`, log `islhandles.log`; the extension is out of sample because the
2 cm threshold was fixed in September):

| two planes | S GENIE (raw) | S NuWro (raw) | B Run 3 (raw) | B ext (raw) | B full (raw) | FoM₁₅ ext / full |
|---|---|---|---|---|---|---|
| no handle | 1.39 (72) | 0.147 (253) | 4.15 (10) | 7.56 (28) | 6.21 (38) | 0.434 / 0.477 |
| Michel tag | 1.08 (56) | 0.107 (184) | 2.07 (5) | 4.32 (16) | 3.43 (21) | 0.448 / 0.494 |
| **photon gap > 2 cm** | **1.33 (69)** | **0.142 (245)** | 3.32 (8) | **4.05 (15)** | **3.76 (23)** | **0.555 / 0.572** |
| Michel + gap | 1.02 (53) | 0.104 (179) | 2.07 (5) | 2.97 (11) | 2.62 (16) | 0.499 / 0.525 |

The gap keeps 96% of GENIE and 97% of NuWro signal and halves the extension's background
(+28% FoM out of sample, purity 26%). On the Run 3 list it removed 2 of 10 events where
the pooled rate predicts 3.9, a 1.3σ fluctuation that hid it. The events it removes are
mostly not photons: at two planes 11 of 15 are hadrons (5 p, 4 π±, 1 μ, 1 e)
reconstructed as a shower at the vertex. The threshold sits on a plateau
(FoM 0.54–0.60 from 1 to 8 cm). 5 cm scores highest on 14 events but costs 18% of
signal, so the pre-set 2 cm is kept. At one plane the gap gives 0.330 → 0.406 on the
extension.
The **Michel tag** is marginal on its own (+3% out of sample, costing 22% of signal)
and lowers the FoM when added to the gap, so it stays off.
Σ-mass windows were tested and cost more signal than they are worth.

**The photon is otherwise unused in this chain** — `NShowers` already requires a
shower, so `requirePhoton` is a no-op. Σ mass and photon energy remain untouched.

## POT normalisation

Every $S$ and $B$ in this repo is the expected count in the full RHC dataset,
**1.1526e21 POT** (`/data/uboone/pot.txt`, RunAll_rhc).

| sample | format | generated POT | scale |
|---|---|---|---|
| GENIE hyperon signal, `sig_genie.list` (4851 files) | hypana | 5.9770e22 | 1.9284e-2 |
| RHC overlay background, `bkg_overlay.list` (14,164 files) | hypana | 7.0488e21 | 1.6352e-1 |
| &nbsp;&nbsp;— Run 3 slice only, `bkg_overlay_run3.list` (3491 files) | hypana | 2.7795e21 | 4.1468e-1 |
| NuWro Σ⁰ signal, `sig_nuwro.list` (2500 files) | hypana | 1.3394e24 | 8.6053e-4 |
| NuWro Σ⁰ signal, `Hyperons_Nuwro_run3b_WithCR_pelee_2500.root` | PeLEE | 1.3394e24 | 8.6053e-4 |
| &nbsp;&nbsp;— superseded Run 1-EXT sample, `sig_nuwro_run1ext.list` (2981 files) | hypana | 3.4141e23 | 3.3760e-3 |
| RHC overlay (Run 1+2+3+4a+4b+4c) | PeLEE | 1.5463292e22 | 7.4538e-2 |
| RHC dirt (runs 1, 2b, 3b, 4a, 4b; old flux, × flux map) | PeLEE | 1.7409194e21 | 0.6621 |
| run3b NuMI EXT (beam-off) | PeLEE | 32,649,128.65 gates | 0.7497 |
| run3b NuMI EXT, ext_unbiased skim (hypana, in production) | hypana | 162,912 gates | 150.2 |
| &nbsp;&nbsp;— superseded: Run 1 "dirt" (it is FHC) | PeLEE | 4.21274e20 | 2.736 |

**The NuWro scale also carries a flux weight of 0.677** on top of the POT scale: see
"NuWro is on the old NuMI flux" below. The table lists the POT scale alone.

**Never extrapolate POT from a file subset.** An early NuWro scale was taken from the
first 300 files and was 52% high; POT/file is not uniform across a list. For the run3b
sample the hypana `ana/MetaTree` sum (1.3394e24 over 2500 files) and the PeLEE `SubRun`
sum over the merged ntuple (1.33940e24, 249832 entries) agree to 4 digits.

The hypana overlay used to be only the **Run 3** slice (2.78e21 of the 9.94e21 the
PeLEE chain uses), which is why its per-event weight was 0.41 and why nothing below
~2 background events could be resolved. Since 2026-09-24 `bkg_overlay.list` also holds
the RHC overlay extension (Runs 1/2/3/4a). On 2026-10-16 the 140 files recovered by
the residue pass were added, completing it at 10,673 files, 4.2693e21 POT and
1,116,279 events, so the list is **14,164 files at 7.0488e21 POT** -- 71% of the
full RHC, a per-event weight of 0.164, and a background MC statistical error 0.62x
its former size. That is +1.13% of POT against the 14,024-file state, so background
yields quoted before that date are 1.11% high. POT/event agrees between the two
slices (3.825e15 extension vs 3.829e15 existing), so they combine without a
relative normalisation. Paths:
`/data/uboone/SigmaFiles/prod_numi_extunbiased_nu_overlay_g4.10.4_flux_run3_rhc_pandora_reco2_reco2/hypana/`.

### The NuWro signal sample (run3b, 2026-10-15)

The Run 1-EXT `prod_hyperon_sigmazero_reco1_prodnuwro_mcc9` sample is **retired**. Its
replacement is generated and overlaid locally by `~/MicroBooNE/production/nuwro_prod/`:
2500 jobs × 100 NuWro Σ⁰ events on swizzled NuMI **run3b ext_unbiased** data, full
gen→reco2 chain, then hypana + PeLEE ntuples in
`/data/uboone/prod_nuwro_sigmazero_overlay_numi_rhc_run3b_extunbiased_mcc9/ana/`.
All 2500 jobs completed (249,832 events, 168 dropped by the reco2 watchdog), and
`check_ext_unique.py --outputs` confirms 250,000 distinct EXT events with none reused.

| | old (Run 1 EXT) | new (run3b EXT) |
|---|---|---|
| files / events | 2981 / 64,716 | 2500 / 249,832 |
| generated POT | 3.4141e23 | 1.3394e24 |
| true Σ⁰ CC ν̄μ | 17,714 (27.4%) | 57,841–58,545 (23.5%) |
| hypana efficiency at `NTracksCut` | 20.01% | 20.43% |
| raw signal surviving | 3,545 | 11,958 |

The cosmics now come from the same run period as the overlay background, and the
selection behaves the same on them (efficiency moves 0.4 points, candidate purity
improves slightly), so the gain is 3.4× the surviving MC statistics.

**NuWro is on the old NuMI flux.** NuWro's flux files were converted from the old
(bugged-geometry) gsimple flux, while the GENIE overlay, the dirt and the data use the
current g4 10.4 flux. Every NuWro event therefore carries the official MicroBooNE
"New/Old" ratio weight in (E_ν, angle to the beam axis). Measured over the finished
production the mean weight is **0.6774**, and — the question worth asking, since hyperon
events sit higher in energy — it is **0.6770 for the Σ⁰ CC signal subset** (⟨E⟩ 1.67 GeV
against 1.54 GeV for the rest), so the correction is flat to 0.1%. The spread is wider
for signal (rms 0.297 vs 0.226), costing a little effective statistics: 49.9k of 58.5k.
`SigmaSelection.C` computes the weight on the fly from `nu_e`/`true_nu_p*`;
`HypCutFlow.C` takes the per-event table as its `fluxWeights` argument
(`nuwro_prod/flux_rw/nuwro_rhc_run3b_reco2_fluxweights.tsv`). Details and the FHC
numbers: `nuwro_prod/flux_rw/README.md`.

### What the new sample does to the selection

Both chains were run twice, with the same code and the same background, changing only
the signal sample (and then the BDT). The selection itself is unaffected by the move to
run3b cosmics; what changes is statistics and normalisation.

PeLEE chain, scaled to 1.1526e21 POT:

| | old sample, old BDT | new sample, old BDT | new sample, **retrained** |
|---|---|---|---|
| raw signal at `AllEvents` | 6,446 | 21,521 | 21,521 |
| efficiency at `NTracksCut` | 50.05% | 50.79% | 50.79% |
| BDT training tree `TreeS` | 3,226 | 10,931 | 10,931 |
| BDTG cut | +0.562 | +0.502 | **+0.616** |
| S | 4.200 | 2.639 | 2.850 |
| B | 777.9 | 967.1 | **749.9** |
| purity | 0.537% | 0.272% | 0.379% |
| signal eff of preselection | 38.6% | 41.4% | **44.7%** |
| S/√(S+B) | 0.1502 | 0.0847 | 0.1039 |

S is 1.7x smaller than it was purely through normalisation — 0.677 (flux weight) x 0.864
(this production's lower Σ⁰ fraction) — so comparing the first column with the third
measures the sample's normalisation, not its selection. The middle column is the right
baseline for the retraining, and against it the retrained BDT keeps MORE signal with LESS
background: +23% in S/√(S+B), purity +39%.

**The retraining confirms that this chain was signal-statistics-limited.** Scores are
cross-validated throughout (even→foldB, odd→foldA), so this is out-of-sample:

| | AUC | B/B₀ at 50% signal eff | at 40% | at 20% |
|---|---|---|---|---|
| old sample, old weights | 0.8429 | 0.0681 | 0.0414 | 0.0125 |
| new sample, old weights | 0.8333 | 0.0688 | 0.0432 | 0.0132 |
| new sample, retrained | **0.8656** | **0.0467** | **0.0293** | **0.0084** |

At matched signal efficiency the retrained BDT rejects ~32% more background. The gain
folds together 3.4x the training signal and training on cosmics from the matching run
period; separating the two needs a third training on a 3,226-event subsample, which has
not been run. Weights before the retraining are kept in `pelee/pelee_bdt_prerun3b_backup/`.

hypana chain (`HypCutFlow.C`, signal pass only), with the flux weight applied to all
58,545 signal events at ⟨w⟩ = 0.6770: S = 34.105 at `AllEvents` → **6.970** at
`NTracksCut`, efficiency 20.44%, against 11.968 → 20.01% for the retired sample, on
3.4x the raw survivors (11,958 against 3,545). Candidate purity improves slightly
(muon 95.2→96.0%, proton 66.6→67.2%, pion 43.1→43.5%, photon 73.1→75.1%).

### The three signal samples disagree on the absolute Σ⁰ rate

Measured with one truth definition (`mc_hyperon_pdg==3212 && mc_nu_pdg==-14 && CC`) and
one fiducial volume, Σ⁰ events per 1.1526e21 POT with the true vertex in the FV:

| sample | in-FV Σ⁰ per full-RHC POT | vs the plain overlay |
|---|---|---|
| plain RHC overlay, `bkg_overlay.list` (unfiltered, trustworthy POT) | 83.3 | 1.00 |
| GENIE hyperon sample, `sig_genie.list` | 110.1 | 1.32 |
| NuWro run3b, POT scale only | 18.9 | 0.23 |
| NuWro run3b, with the 0.677 flux weight | 12.8 | **0.15** |

The fiducial fractions are the same in all three (37–38%), so this is not a generation-
volume artefact, and it is not new — the retired sample gives 22 before the flux weight.

**Most of the NuWro gap is hyperon FSI, not normalisation.** The generator configuration
(`nuwro_old/configsets/ultrapure-sigma0-numi-rhc.txt`) sets `dyn_hyp_cc=1` with every
other channel off and, within it, `hyp_lambda=0`, `hyp_sigma_zero=1`, `hyp_sigma_minus=0`:
NuWro generated the **Σ⁰ channel and nothing else**, so its quoted POT/event is the
exposure per Σ⁰ *produced*, which is exactly the normalisation a Σ⁰ signal sample wants.
Everything else in the sample is FSI: of the ν̄μ CC events only 25.3% still carry a Σ⁰,
27.9% carry a Λ (Σ⁰N → ΛN in the cascade), 10% a Σ±, and 37.3% no strange particle at all.
The MicroBooNE Λ paper (arXiv:2212.07888) states the other half of this outright — "the
NuWro predictions include the effect of final state interactions, while GENIE does not
take them into account". GENIE's hyperons leave the nucleus untouched, so for GENIE
produced = surviving, and the comparison that matters is Σ⁰ **produced**:

| | in-FV Σ⁰ produced per full-RHC POT |
|---|---|
| GENIE overlay (no hyperon FSI, so = surviving) | 83.3 |
| NuWro, POT scale only | 79.1 |
| NuWro, with the 0.677 flux weight | 53.6 |

That is agreement to 0.95x before the flux weight and 0.64x after — the ordinary
generator spread, not a factor of six. The same FSI explains the sign flip between
channels: in Fig. 2 of the Λ paper NuWro sits **above** GENIE (≈2.1 against ≈1.5
×10⁻⁴⁰ cm²/Ar), because FSI *adds* Λ, while for Σ⁰ it is far below, because FSI *removes*
Σ⁰. Both follow from the one modelling difference.

**A NuWro control run closes the normalisation question.** The same flux, target and
config with `dyn_hyp_cc=0 dyn_qel_cc=1` — ordinary CC quasi-elastic instead of the hyperon
channel — gives 7.539e15 POT per event over the box. The overlay's own CC QEL rate, 20,277
events with the true vertex in the FV in 7.1479e20 POT, converts to 8.835e15 on the same
footing (37.0% box-mass FV fraction, 0.677 old-flux factor). **NuWro's generic event rate
is therefore within 17% of GENIE's**, which is the ordinary CCQE model difference, so its
flux × target normalisation is sound and the hyperon numbers can be taken at face value.
Referenced to each generator's own CCQE rate, NuWro's Σ⁰ *production* is 0.81× GENIE's.

The 6.5× gap in *surviving* Σ⁰ then decomposes exactly:

| factor | size |
|---|---|
| hyperon FSI — NuWro keeps 23.9% of the Σ⁰ it produces, GENIE keeps all | **4.19×** |
| old → current NuMI flux weight | 1.48× |
| residual production-rate difference | 1.05× |
| product | 6.51× (measured 6.51×) |

**What is still open** is the residual ~1.3–1.6x between the GENIE hyperon sample, the
GENIE overlay and the Λ paper's own GENIE prediction — the 1.61x of the Λ cross-check.
Practical consequence: a Σ⁰ yield quoted from the NuWro sample is a yield *with* hyperon
FSI and one from the GENIE samples is a yield *without* it, and the difference between
them is a model choice of roughly a factor of four, not a bug to be fixed. The measured
Λ cross section, 2.0 +2.2/−1.7 ×10⁻⁴⁰ cm²/Ar, brackets both generators and cannot
currently distinguish them.

### Sample bookkeeping notes

- The **NuWro sample is Σ⁰-filtered** (Σ⁰ 24.9%, Λ 27.9% at the hyperon level, and 37.3%
  of events lose the hyperon to FSI) where the GENIE sample is Λ-dominated (Λ 50.6%,
  Σ⁰ 7.4%). It can validate Σ yields but never a Λ rate.
- hypana and PeLEE disagree on raw signal count (17714 vs 6446) only because PeLEE
  quotes true-Σ⁰-in-FV as its efficiency denominator. At the slice stage they agree
  event-for-event (9426 vs 9427).
- The Lambda cross-check (`LambdaSel.C`) reproduces the published Lambda CCQE chain on
  our files: efficiency 0.79× theirs (expected — our pair selector is provisional and
  their decay BDT is replaced by cuts), but **signal normalisation 1.61×**, which is
  GENIE-against-GENIE and still unexplained. Candidates in order: GENIE version/tune
  (hyperon production is axial-mass sensitive), NuMI flux version, the equivalent-
  unfiltered exposure quoted for the hyperon-filtered production, FV definition.

## Production status — extending the hypana background

The hypana chain is MC-statistics limited: 10 background events at the two-plane point.
The fix is more RHC overlay in hypana format. Runs 1, 2 and 3 RHC overlay **reco1** are
being transferred to `/data/uboone/SigmaFiles/`:

| sample | directory |
|---|---|
| Run 1 | `NuMI_Nuw_Flux_Run1_RHC_Reco1_jh/` |
| Run 2 | `Run2_RHC_Nu_Overlay_reco1/` |
| Run 3 | `Run3_RHC_Nu_Overlay_Reco1_fix_gen_stash/` |

These are **reco1c**, not reco2: they carry `recob::Wire`, `gaushit`, optical, CRT and
full truth, but **no Pandora products** (104 branches for Runs 1–2, 124 for Run 3 which
adds the CRT; the working reco2 input has 272). So they need reco2 before
HyperonProduction can run:

```
reco1c ──[reco2_local, 5-FCL chain]──▶ reco2 ──[run_HyperonProduction.fcl]──▶ *_hypana.root
```

Both stages already exist locally (`~/MicroBooNE/production/reco2_local/` and
`~/MicroBooNE/HypProductionNuMIFlux/slurm/`), running `uboonecode v08_00_00_85 -q
e17:prof` inside the SL7 apptainer on SLURM.

**The chain is verified end to end** (2026-09-13, Run 1 file, run 5435 subrun 18):
reco1c → reco2 in 4 min wall (5 steps, clean) → HyperonProduction in 14 s → an
ntuple with **142 branches, byte-identical in branch set to the existing Run 3 hypana
files**, 11 `ct_*` island branches with the wire windows populated (100 wires,
40–78 hits, 3 vertex wires per event). Events and subruns are preserved at every
stage, and `sumdata::POTSummary` survives into reco2. No run-period problem: the
same fcl chain works for a Run 1 file despite having been written for the Run 1-era
prodnuwro sample, because MicroBooNE's period dependence comes from the conditions
database keyed on run number.

### How much background this actually buys

**Measured, not assumed.** reco1 holds ~8.6 MB/event (five files sampled: 8.5, 8.8,
8.3, 7.8, 8.6) and ~40 events/file, and the POT is ~2.6e15/event for Run 1 against
3.83e15 for the existing Run 3 overlay.

| | events | POT | vs existing |
|---|---|---|---|
| existing hypana overlay (Run 3) | 725,995 | 2.7795e21 | — |
| extension processed 2026-09-24 (10,533 files) | 1,095,671 | 4.1906e21 | +151% |
| residue pass recovered 2026-10-16 (140 files) | 20,608 | 0.0787e21 | +1.1% |
| **extension complete** (10,673 files) | **1,116,279** | **4.2693e21** | **+154%** |
| combined `bkg_overlay.list` (14,164 files) | 1,842,274 | 7.0488e21 | **2.54x** |
| to reach the full RHC 9.94e21 | +760,000 | +2.9e21 | needs ~7 TB more reco1 |

So the ×3.6 ceiling implied by the PeLEE/hypana POT ratio requires ~22 TB of reco1
and a few thousand core-days of reco2. What has landed so far adds 5%, which does
**not** change the two-plane background estimate (10 MC events → ~10.5). Decide the
target exposure before committing the CPU: doubling the background is worth ~6 TB and
~500 core-days; the full RHC exposure is a different scale of undertaking. Disk is not
the constraint (81 TB free on `/data`), CPU is.

Notes for the full production:

- **`jobtag` collides.** `slurm_reco2.sbatch` derives the working directory from
  `basename(dirname(input))`, which for this path layout is the last two digits of the
  run — `54` for run 4954 *and* 5054. Concurrent tasks sharing a tag share a CWD, and
  art's `mem.db`/`time.db` then collide ("disk I/O error"). Append
  `$SLURM_ARRAY_TASK_ID` to `WORKDIR` before submitting at scale.
- **`submit_slurm.sh` rewrites one shared manifest in place**, and a pending array reads
  its manifest at launch. Use a per-sample manifest and call `sbatch` directly when
  another array is queued.
- Its `find` pattern is `*reco1c.root`; these files end `_reco1c_<timestamp>_merged.root`,
  so pass an explicit list with `-l`.
- **reco2 is self-sufficient**: it keeps `sumdata::POTSummary` and the
  `recob::Wire`/`gaushit` collections the island test needs (confirmed on the test
  output), so reco1 can be deleted once reco2 is verified. Verify by **POT
  conservation**, not file count — a dropped subrun gives a full-looking directory with
  a quietly wrong normalisation.
- **Cost per file**: 4 min wall for a 4-event file, most of it fixed cost (container
  start, service construction, Pandora geometry). The marginal per-event cost is much
  lower and is being measured on the 36-event test files; do not extrapolate from the
  small-file number.

## Old analysis — OutputTree / TMVA (archived)

The original selection on `HyperonProduction`/`OutputTree` ntuples with the TMVA BDT.
Self-contained in `archive/`; all inputs/outputs resolve relative to that directory.

```bash
bash archive/runSigma.sh           # selection -> archive/SigmaHists.root (+ cut-flow)
bash archive/runSigma.sh --plot    # also plots -> archive/figs/
```

BDT bootstrap (the trained weights were removed; regenerate if needed):

```bash
bash archive/runSigma.sh                 # 1st pass -> archive/TMVA_input.root
root -b -l -q archive/RunTMVA.C          # train    -> archive/SigmaDataset/weights/
bash archive/runSigma.sh                 # 2nd pass -> activates the BDTg cut
```
