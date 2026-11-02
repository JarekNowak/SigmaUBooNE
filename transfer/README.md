# Transferring the reco2 samples from FNAL

These tools copy the samples listed in `../Docs/RequiredSamples.md` from FNAL dCache to `/data` on
lapx, for the Λ-style (hypana) chain. The files are on tape at FNAL, and a lapx user token cannot
recall from tape. The transfer therefore has two halves:

1. **On a uboone gpvm:** create a SAM sub-definition and prestage it from tape to dCache disk.
2. **On lapx:** copy the files that are online over xrootd (resumable), then make hypana ntuples.

Beam-on data (block `D`) is blind: both scripts refuse it without `--unblind`, which is only to be
used once unblinding is approved.

## Files

| file | where | what |
|---|---|---|
| `samples.tsv` | — | the samples: key, block (R RHC required, F FHC required, E EXT decision, O optional, D beam-on), blind flag, tier (reco2, or reco1 that needs our reco2), SAM dimensions |
| `gpvm_prestage.sh` | gpvm | create `<SAMUSER>_sigma_<key>[_N]` and run `samweb prestage-dataset` in the background |
| `fetch_sample.sh` | lapx | list a sample in SAM, find the dCache copies, copy those online; resumable, `.part` until complete |
| `fetch_worker.sh` | lapx (container) | the SAM/xrootd part of `fetch_sample.sh`; not run directly |
| `fetch_loop.sh` | lapx | renew the token and re-run `fetch_sample.sh` every 30 min until a sample is complete |
| `status.sh` | lapx | local against listed files per sample (`-s`: also the SAM size) |
| `sample_lookup.sh` | — | shared helpers; `SAMUSER` (default `janow`), and the destination `$DATA_ROOT/<category>/<tier>/<key>` (default root `/data/uboone`) |

## Step 1: prestage (FNAL gpvm)

```bash
ssh <user>@uboonegpvm03.fnal.gov
git clone -b pelee-analysis <Sigma repo> && cd Sigma/transfer   # or copy this directory
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup sam_web_client
htgettoken -a htvaultprod.fnal.gov -i uboone       # if samweb asks for credentials

./gpvm_prestage.sh -l                       # keys, blocks, blind flags
./gpvm_prestage.sh -s dirt_rhc_run3b        # create the sub-definition only, print its size
./gpvm_prestage.sh dirt_rhc_run3b           # create + prestage (background, prestage_logs/)
./gpvm_prestage.sh ext_rhc_run3b 20000      # first 20000 files only: janow_sigma_ext_rhc_run3b_20000
./gpvm_prestage.sh -b R                     # every RHC-required key
```

Set `SAMUSER` if your FNAL user name is not `janow`. Recalls of several TB take hours to days.
dCache keeps prestaged files on disk only for a while, so prestage in priority order and copy
each block soon after:

1. `R`: dirt, hyperon-enriched detector variations, overlay remainder;
2. `F`;
3. `E` (once the EXT strategy is decided);
4. `O`.

Progress: `samweb project-summary` (the project name is in the log), or `./fetch_sample.sh <key> -M 0`
on lapx, which reports how many files are online.

## Step 2: copy (lapx)

```bash
cd /home/t2k/nowak/MicroBooNE/Sigma/transfer
../../production/nuwro_prod/get_token.sh    # first time: open the CILogon URL; then renews silently ~1 week
./fetch_sample.sh dirt_rhc_run3b -M 0       # list only: local / online / tape-only / stat errors
./fetch_sample.sh dirt_rhc_run3b            # copy everything online (4 parallel; -P to change)
./fetch_sample.sh ext_rhc_run3b -N 20000    # the limited sub-definition from step 1
nohup ./fetch_loop.sh dv_hyp_rhc_CV dv_hyp_rhc_LYDown > loop_dv.log 2>&1 &   # long transfers
./status.sh                                 # progress of every sample
```

Each sample lands in `/data/uboone/<category>/<tier>/<key>/`, where the category is dirt, detvars, EXT, overlay, hyperon or beam_on, and the tier is reco2 or reco1. For example, `/data/uboone/dirt/reco2/dirt_rhc_run3b/`. Each directory holds:
- `pnfs.list`: the SAM file list;
- `files.list`: the complete local files, which is the input list for the next step;
- `.xrdcp_errors.log`.

Re-running is safe: files already here are skipped. `--relist` rebuilds `pnfs.list` from SAM, after
the sub-definition changes.

**"stat error"** means the token expired or is missing. Run `get_token.sh` again.

**"tape-only"** means the file is not prestaged yet. Run step 1, or wait.

**Disk.** `fetch_sample.sh` stops if `/data` has less than 5 TB free. The required RHC and FHC blocks
are about 73 TB against about 77 TB free on 2026-11-01. Check with `df -h /data`, and delete or
skim the reco2 after the hypana ntuples are made (keep the EXT reco2 and skims: they may be
reprocessed).

## Step 3: hypana ntuples

- **reco2 samples:** run HyperonProduction on `files.list`, the same way as the existing hypana
  productions (`production/nuwro_prod/ana_file.sh`, `submit_ana.sh`).
- **reco1 samples** (`ovl_*_reco1`): run our reco2 first, as for the RHC background extension (see
  memory `hypana-background-extension`: STEP_TIMEOUT of 10 h, residue passes).

Then add the ntuples to `pelee/lists/` and the POT and weights:
- **POT:** measure it file by file from the reco2 `POTSummary`.
- **Weights:** dirt and the hyperon-enriched samples need the flux and hyperon cross-section weights
  (`Docs/RequiredSamples.md`, Caveats).

## Checks after each sample

- `./status.sh`: local equals listed.
- File count and events against the SAM summary (`./status.sh -s`).
- For the hyperon-enriched detector variations: the number of **distinct** events. The old
  hyperon-enriched sample turned out to be 754 events duplicated about 20 times.
