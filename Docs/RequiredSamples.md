# Required reco2 samples for the final Σ⁰ analysis

*Compiled 2026-11-01. The SAM definitions were found with `samweb` (read-only, in the SL7 container
on lapx). File counts, event counts and sizes come from `samweb list-files --summary`
(TB = 10¹² bytes).*

## Scope

The Λ-style (island) chain gives the analysis result. It runs on HyperonProduction ntuples
(hypana), which are made from **reco2 artroot**: the island test needs the raw wire activity, which
PeLEE ntuples don't carry. So every sample that enters the Λ-style result must be available here
at reco2 (or at reco1, if we run reco2 ourselves).

The PeLEE chain is the cross-check. It needs only PeLEE ntuples, which already exist for data, EXT,
dirt, overlay and detector variations in both beams (`/data/uboone/temp/custom_pelee_ntuples`).

The data, EXT and dirt definitions below were found by following the SAM parentage of those PeLEE
ntuples up to their reco2 parent files. They are therefore the same events as the PeLEE chain.

**Blindness.** The analysis is blind. Beam-on data is listed last (Section 6). It is not
transferred or processed through hypana until the selection, the background model and the
systematics are frozen and unblinding is approved. Everything before Section 6 is MC, beam-off or
dirt.

**Transfer.** Use the tools and instructions in [`../transfer/README.md`](../transfer/README.md).
They prestage on a gpvm (`gpvm_prestage.sh`), then copy on lapx (`fetch_sample.sh`,
`fetch_loop.sh`). Every sample below has a key in `transfer/samples.tsv`.

---

## 1. Already here

| sample | SAM definition / source | held |
|---|---|---|
| NuWro Σ⁰ signal, RHC run3b | own production (`prod_nuwro_sigmazero_overlay_numi_rhc_run3b_extunbiased_mcc9`) | all, hypana done |
| GENIE Σ⁰ QE signal, RHC run3b (Pais cross section, weighted) | own production | all, hypana done |
| GENIE hyperon-enriched, RHC run3b | `cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2` | 4,851 / 4,851 files |
| Overlay, new flux, Run3 RHC reco2 | `prod_numi_extunbiased_nu_overlay_g4.10.4_flux_run3_rhc_pandora_reco2_reco2` | 3,491 / 3,491 |
| Overlay, new flux, Run1 RHC reco1 (our reco2) | `NuMI_Nuw_Flux_Run1_RHC_Reco1_jh_reco1_reco1` | 2,242 / 2,265 |
| Overlay, Run2 RHC reco1 (our reco2) | `Run2_RHC_Nu_Overlay_reco1_reco1_reco1` | 3,915 / 13,843 |
| Overlay, Run3 RHC reco1 (our reco2) | `Run3_RHC_Nu_Overlay_Reco1_fix_gen_stash_reco1_reco1` | 4,456 / 7,159 |
| Overlay, Run4a RHC reco1 (our reco2) | `run_4a_numi_rhc_g4_10_4_nu_overlay_reco1_reco1` | 60 / 1,130 |
| EXT run3b reco2 | `prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2` | 16,570 / 3.33M events |

The hypana background overlay is 14,164 files, 7.05×10²¹ POT: 71% of the RHC overlay exposure.

## 2. RHC: required for a Λ-style result

| # | sample | SAM definition | files | events | TB | why |
|---|---|---|---|---|---|---|
| R1 | **Dirt, run3b** (old flux, FluxRW) | `prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b_v28_all` | 8,303 | 392k | 2.39 | dirt is only bounded today (≤ 8.7% of B) |
| R2 | **Detector variations, hyperon-enriched**, 10 sets | `cthorpe_make_hyperon_events_numi_rhc_run3b_detvar_reco2_<V>_reco2` | ~4,910 each | ~287k each | 26.4 total | detector systematics on signal and on the real-Λ background (46% of B) |
| R3 | **Overlay remainder** (reco1, our reco2) | the three reco1 definitions in Section 1, missing files only | ~13,700 | ~1.5M | ~13.5 | background MC statistics: B rests on 12 raw events (±29%) |
| R4 | NuWro hyperon overlay, run3b | `cthorpe_prod_numi_uboone_overlay_rhc_mcc9_run3b_v51_nuwro_hyperon_nuwro_reco2_reco2` | 496 | 28.6k | 0.19 | NuWro-model hyperon background; content to check (may cover part of the planned NuWro S1 sample) |

`<V>` = CV, LYDown, LYRayleigh, LYAttenuation, Recomb2, SCE, WireModX, WireModYZ, WireModThetaXZ,
WireModThetaYZ.

**EXT (R5): decision needed.** The current EXT estimate is factorised (0.34 ± 0.15, 14% of B).

- **Option A, direct estimate:** a large EXT reco2 set (each row below is a separate SAM definition).

  | | SAM definition | files | events | TB |
  |---|---|---|---|---|
  | Run3b remainder | `prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2` | 77,400 | 3.33M | 18.05 |
  | Run2b | `prod_extnumi_swizzle_crt_inclusive_v4b_offbeam_run2_reco2_new2_run2_reco2_all` | 64,869 | 1.05M | 5.58 |
  | Run4b | `Run_4b_NuMI_beam_off_Pandora_Reco2_run4b_reco2_all` | 4,278 | 1.54M | 8.23 |
  | Run1 (RHC runs only) | `prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2` | 36,052 | 1.32M | 6.71 (whole set) |

- **Option B, data-driven:** a cosmic estimate from sidebands. This needs beam-on data, so it
  waits for unblinding or for a sideband definition agreed under the blinding policy.

## 3. FHC

| # | sample | SAM definition | files | events | TB |
|---|---|---|---|---|---|
| F1 | **Overlay Run1+2, new flux, reco1** (our reco2, as RHC) | `prodgenie_numi_fhc_nu_overlay_newflux_run12_reco1` | 14,900 | 1.47M | 12.13 |
| F2 | **Dirt Run1** (old flux) | `prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all` | 11,634 | 572k | 3.15 |
| F3 | **EXT Run1** (FHC runs; also needed as overlay input for our FHC productions) | `prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2` | 36,052 | 1.32M | 6.71 |
| F4 | **Detector variations, hyperon-enriched run1**, 10 sets | `cthorpe_make_hyperon_events_run1_detvar_<V>_hyperon_reco2_reco2` | ~1,900 each | ~95k each | 7.6 total |
| F5 | GENIE hyperon-enriched run1 | `cthorpe_make_hyperon_events_reco2_numi_fhc_run1_pt1_hyperon_reco2_reco2` | 1,889 | 94k | 0.69 |
| F6 | GENIE hyperon run1 v51 | `cthorpe_prod_numi_uboone_overlay_fhc_mcc9_run1_v51_GENIE_hyperon_real_GENIE_reco2_reco2` | 484 | 24k | 0.15 |
| F7 | NuWro hyperon run1 v51 | `prod_numi_uboone_overlay_fhc_mcc9_run1_v51_nuwro_hyperon_nuwro_reco2_reco2` | 495 | 25k | 0.15 |

F3 overlaps with the run1 EXT prestage needed by `samples/fhc_run1.conf`, so one transfer serves both.

## 4. Optional: other run periods, inclusive variations

These are needed only if the MC must cover each run period rather than one period per beam
(run3b RHC, run1 FHC) scaled to all of them.

| sample | SAM definition | files | TB |
|---|---|---|---|
| Overlay Run4b RHC reco2 | `run4b_numi_rhc_nu_overlay_Pandora_reco2_v08_00_00_63_run4b_reco2` | 4,506 | 4.22 |
| Dirt Run4a RHC | `run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_reco2` | 249 | 0.18 |
| Dirt Run4b RHC | `numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_reco2` | 1,142 | 0.62 |
| Overlay Run4c FHC | `numi_run4c_fhc_nu_overlay_pandora_unified_reco2_run4c_reco2` | 2,276 | 2.31 |
| Overlay Run4d FHC | `run4d_numi_fhc_nu_overlay_pandora_unified_reco2_run4d_fhc_reco2` | 3,431 | 3.55 |
| Overlay Run5 FHC | `NuMI_Run5_FHC_nu_overlay_pandora_Reco2_run5_reco2` | 4,757 | 5.44 |
| EXT Run2a FHC | `prod_extnumi_swizzle_inclusive_v4_run2a_reco2_run2a_all_reco2` | 45,345 | 7.84 |
| EXT Run4c FHC | `run4c_numi_beam_off_pandora_reco2_run4c_reco2_all` | 3,792 | 5.99 |
| EXT Run4d FHC | `run4d_numi_beam_off_pandora_reco2_run4d_reco2_all` | 10,409 | 19.75 |
| EXT Run5 FHC | `run5_numi_beam_off_pandora_reco2_all_run5_reco2_all` | 10,796 | 13.76 |
| Detector variations, inclusive run3 RHC (8 found) | `numi_run3_standard_nu_overlay_{cv,LYDown,Recomb2,SCE,WireModThetaXZ}_reco2_v08_00_00_54_run3_reco2`, `numi_run3_standard_nu_overlay_{LYrayleigh,LYattenuation,WireModThetaYZwithsplines}_reco2_run3_reco2` | ~19,800 each | ~55 |
| Detector variations, inclusive run1 FHC | `prodgenie_numi_nu_overlay_v08_00_00_53_{CV,SCE,WireModThetaXZ}_300k_reco2_run1_reco2` | 9,047 (CV) | 1.86 (CV) |

The inclusive detector variations matter only if the non-hyperon background needs its own
detector systematic. The hyperon-enriched sets (R2, F4) cover the signal and the real-Λ background.
WireModX and WireModYZ (run3), and the remaining run1 variations, were not found under these names.
They can be traced by parentage from their PeLEE ntuples.

## 5. Own productions (generated here, not transferred)

| sample | status |
|---|---|
| Fixed-GENIE Σ⁰ QE (`genie::HypQELTracePXSec`), FHC run1 then RHC run3b | generating; replaces the weighted GENIE QE signal |
| NuWro Σ⁰ signal, FHC run1 | hepmc done; overlay needs F3 |
| NuWro hyperon background S1 (QE Λ, Σ⁻) and S2 (associated strangeness), both beams | planned (`production/nuwro_prod/PLAN_NUWRO_BKG.md`); check R4/F7 first |

## 6. Beam-on data (blind: transfer only after unblinding is approved)

RHC exposure 1.15×10²¹ POT; FHC exposure 1.10×10²¹ POT.

| beam | run | SAM definition | files | events | TB |
|---|---|---|---|---|---|
| RHC + FHC | Run1 (split by run number) | `prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2` | 37,123 | 895k | 4.71 |
| RHC | Run2b | `prod_numi_swizzle_inclusive_v4b_run2_beamon_run2_reco2_run2_beam_good_reco2` | 62,611 | 283k | 1.57 |
| RHC | Run3b | `prod_numi_mcc9_v08_00_00_45_run3b_run3b_reco2_beam_good_reco2` | 44,597 | 1.14M | 6.25 |
| RHC | Run4a | `run4a_NuMI_RHC_beam_on_pandora_reco2_reco2_reco2_beam_good` | 166 | 43k | 0.23 |
| RHC | Run4b | `Run_4b_NuMI_beam_on_Pandora_Reco2_run4b_reco2_beam_good` | 3,738 | 843k | 4.55 |
| FHC | Run2a | `prod_numi_swizzle_inclusive_v4_run2_reco2_run2a_beam_good_reco2` | 32,556 | 764k | 3.73 |
| FHC | Run4c | `Run_4c_NuMI_beam_on_Pandora_Reco2_run4c_reco2_beam_good` | 4,228 | 393k | 2.13 |
| FHC | Run4d | `run4d_NuMI_beam_on_pandora_unified_reco2_run4d_reco2_beam_good` | 2,626 | 229k | 1.23 |
| FHC | Run5 | `run5_numi_beam_on_pandora_reco2_all_run5_reco2_beam_good` | 4,025 | 749k | 4.04 |

Total: about 28 TB.

- **Before unblinding:** an open data subset may be used for data/MC checks of the island test and
  the sidebands, if the blinding policy allows. PeLEE has open-trigger Run1 ntuples
  (`custom_pelee_ntuples/open_trigger`); their reco2 parents would be the matching subset.
- **Run4a looks incomplete:** the definition holds 43k events against 740k Run4a triggers. Check
  it, or find the complete Run4a set, before unblinding.
- **Run4c (RHC part):** `/data/uboone/pot.txt` lists RHC POT for Run4c. Its RHC runs are in the
  Run4c definition above (labelled FHC in the PeLEE production); split by horn current per run.

## 7. Totals

| block | TB |
|---|---|
| RHC required (R1–R4) | ~42 |
| RHC EXT, option A | ~39 |
| FHC required (F1–F7) | ~31 (incl. EXT run1) |
| Optional (Section 4) | ~120 |
| Beam-on (Section 6, after unblinding) | ~28 |

## Caveats

- **GENIE cross section.** The hyperon-enriched samples (R2, F4–F6) are GENIE with the Pais
  hyperon cross section, so they need the cross-section and polarisation weights (`HypXsecWeight.h`).
- **Duplicated events.** The old RHC hyperon-enriched sample was 754 distinct events, each
  duplicated about 20 times. Check the number of distinct events in R2/F4 before relying on their
  statistics.
- **Dirt flux.** The dirt samples (R1, F2) are on the old flux and need the flux reweighting
  already used in the PeLEE chain.
- **Run1 mixes horn modes.** Run1 definitions contain both FHC and RHC runs; select by run number.
