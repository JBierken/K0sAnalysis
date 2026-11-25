# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository. It is a living document; extend it as the work goes on.

## What this is
The analysis framework for a CMS displaced-vertexing study with V0s (K0S, Λ). It compares data and simulation for the number of reconstructed K0S as a function of the radial displacement Δ2D (`_RPV`) in Z→μμ events, after normalizing simulation to data at small Δ2D. Results go into AN-20-111 and DP-2024-010 (TRK POG). The study was first done on Run 2 (pre-UL/UL); active work is on Run 3 (2022–2026).

## Related repository: the ntuplizer
- **Current (Run 3):** `/user/jbierken/CMSSW_15_0_17/src/V0Analysis/V0Analyzer`, a separate git repo with its own `CLAUDE.md`. It is a CMSSW `EDAnalyzer`, forked from heavyNeutrino/multilep, that rebuilds V0s from `packedPFCandidates`+`lostTracks` (`src/K0Analyzer.cc`). It writes `blackJackAndHookers/blackJackAndHookersTree` plus `hCounter`, `nVertices`, `nTrueInteractions`. CRAB output goes to `/pnfs/iihe/cms/store/user/jbierken/K0sAnalysis/NTuples/MINIAOD/{data,sim}/v<N>/`.
- **Legacy (Run 2):** `ntuplizer/heavyNeutrino` in this repo. It is kept for reference only and is not built or used.
- Jet tagging in the ntuplizer is being reworked (DeepJet → ParticleNet/UParT for Run 3). The ntuplizer output must keep matching what `skimming/skim_ztomumu.cc` reads. See [Interface contract](#ntuplizer--skimmer-interface-contract).

## Pipeline
All steps run with CMSSW sourced (`cmsenv` in `CMSSW_15_0_17`, python3, ROOT, uproot, awkward). Most scripts use relative paths (`sys.path.append('../')`, `python3 ../plotting/...`), so **run them from their own directory**.

1. **Skim** (`skimming/`): `make -f make_skim_ztomumu`, then `./skim_ztomumu <in> <out> <nentries|-1> [veto_tagger]`, or batch via `skim_submit.py -i <samplelist> -o <outdir>`. The selection keeps events with exactly 2 POG-medium muons (pT > 30/25 GeV, |η| < 2.4, opposite charge), `_passTrigger_mm`, |mμμ − mZ| < 10 GeV and no b-tagged jet (loose WP). The veto tagger and all loose WPs come from one per-period table, `bTagSettingsTable` (look it up with `bTagSettings(campaign, year)`). Defaults: **RobustParT for 2022–2023BPix, UParT for 2024–2026, and the DeepJet column for Run 2**. The optional 4th argument (`DeepJet|PNet|RobustParT|UParT|default`) overrides the default. `_nimloth_nBJets_<tagger>` is written for all four taggers (−1 if the tagger isn't evaluated in that ntuple or has no WP). WPs without an official BTV value are marked `PLACEHOLDER`. It adds `_nimloth_*` (event) and `_celeborn_*` (muon) branches. Campaign and year are parsed from the **input file path** (`getYearFromFileName`).
2. **V0 building** (`v0building/`): `v0builder.py -i <in> -o <out> [-s legacy]`, batch via `v0builder_submit.py`, merge via `merge.py`/`mergeeras.py`. Selection is in `v0selections.py::selection_legacy`. It writes trees `laurelin` (K0S) and `telperion` (Λ) with `_RPV`, `_RPVSig`, `_mass`, `_pt`, `_eta`, `_phi`, `_weight`, `_nTrueInt`, …, and copies `hCounter`.
3. **Analysis** (`analysis/`), driven by `mcvsdata_submit.py`:
   ```bash
   cd analysis
   python3 mcvsdata_submit.py -i <filedir> -v run3 -c config_ksvars.py -o <outdir> \
       [-e default|perera|<era> ...] [--runmode local|condor] [--fitstrategy legacy|adaptive]
   ```
   - `-v/--version` ∈ `run2preul`, `run2ul`, `run3`. `-e default` uses the hard-coded era list in the submit script, which gets edited often.
   - `mcvsdata_getfiles.py` maps eras to input files (MC dir names in `mcdirdict`, data in `datadirdict`), luminosities (`tools/lumitools.py`) and the DY cross-section (6688 pb for Run 3).
   - For each era × variable × bkgmode × normalization × binning, the submit script writes `config.json`, `variable.json`, optionally `sidevariable.json`/`normvariable.json`/`yvariable.json`, into `<outdir>/<era>/<var>_<bkg>_<norm>_<bins>/`. It then runs `mcvsdata_fill.py` (→ `histograms.root`) and `../plotting/mcvsdataplotter[2d].py` (→ png).
   - Variable configs: `config_ksvars.py`, `config_ks2dvars.py`, `config_lavars.py`, `config_controlvars.py`, `config_ksvars_detector.py` (each exports a `config` dict). `configs_old/` is unused.
   - `mcvsdata_submit_all.py` is an older Run 2 convenience looper with hard-coded `/pnfs/.../llambrec/` paths.

## `analysis/mcvsdata_fill.py` (core logic)
- `get_histogram(...)` processes one file and returns `(counts, errors, confidence, confidence_error)`.
- MC weights: `_weight`, multiplied by the pileup reweight (`reweighting/pileup/pileupreweighter.py`, `PileupReweighter(campaign, year_pu)`). They are normalized to their own sum, then scaled by xsec×lumi. Per-era names (e.g. `2024C`, `2022postEEF`) are mapped back to the pileup year by stripping era letters.
- **Shared-MC split (2024/2025/2026):** all three years use the 2024 MC sample. `get_mc_split_info` in `mcvsdata_getfiles.py` gives each year in `MC2024_SHARED_YEARS` a `splitrange` = [lo, hi) proportional to its luminosity (from `tools/lumitools.py`); per-era entries use their year's range. `get_split_mask` in `mcvsdata_fill.py` keeps events whose splitmix64 hash of event/run number falls in that range, so the subsets are disjoint and statistically independent. Changing a year's luminosity or the year list reshuffles which events each year gets. Every branch read must be masked with `splitmask` consistently.
- Pileup profiles are read from this checkout (`reweighting/pileup/<year>/Pileup_ratio.root`, made with `reweighting/make_pileup.py -e <year> [-m]`). Until 2026-11-23 the reweighter read them from `CMSSW_14_0_15`, whose `2025/Pileup_ratio.root` was a copy of the 2024 ratio.
- **Background modes:** `None` (plain histogram) or `sideband`. In sideband mode, every bin of the main (and y) variable gets an invariant-mass fit via `fitting/count_peak.py::count_peak_unbinned(mode='hybrid', fit_strategy=...)`, and fit plots go to `--sideplotdir`.
- **Normalization modes:** `lumi`, `yield`, `eventyield` (uses the `nimloth` tree), and `range`. `range` scales total simulation to data in the `normvariable` range, typically `_RPV` ∈ [0, 0.5] cm, and is the main physics result.
- If `--normmode` is not given, simulation xsec and lumi are set to 1.

## Fitting (`fitting/count_peak.py`, `tools/fittools.py`)
The signal model is a Gaussian or double Gaussian on a pol1/pol2 background in the K0S mass window 0.44–0.56 GeV.
- `--fitstrategy legacy` (default) chooses the model from effective-entry thresholds (1000 and 20000).
- `--fitstrategy adaptive` fits all four combinations and keeps the best χ²/ndf. See `ADAPTIVE_STRATEGY_SUMMARY.md`.
- Each TF1 needs a unique name (`funcname`) when several are fitted on one histogram.

## Ntuplizer ↔ skimmer interface contract
`skimming/skim_ztomumu.cc` breaks quietly when this contract drifts. Check it whenever the V0Analyzer changes:
- The tree path is `blackJackAndHookers/blackJackAndHookersTree`. Histograms in that directory are copied over.
- Fixed array sizes must match the ntuplizer headers: `nL_max = 20` (`LeptonAnalyzer.h`), `nJets_max = 100` (`JetAnalyzer.h`).
- Branches read: `_nL`, `_nLight`, `_nMu`, `_lPOGMedium`, `_lPt`, `_lEta`, `_lPhi`, `_lE`, `_lCharge`, `_passTrigger_mm`, `_nJets`, `_jetPt`, `_jetEta`, `_jetPhi`, and the b-vs-all scores `_jetDeepFlavor`, `_jetParticleNet_BvsAll`, `_jetRobustParT_BvsAll`, `_jetUParT_BvsAll` (each optional; −1 means not evaluated). The skim assumes muons are stored first in the lepton arrays.
- To add a tagger, extend `Tagger`/`taggerNames`/`taggerBranches` and add a column to `bTagSettingsTable`. To add a period, add a row there and a match in `getYearFromFileName`. Tagger availability in the v1 ntuples: RobustParT is evaluated in 2022/2023 only, UParT in 2024/2025 only, DeepJet and PNet everywhere.
- Downstream, `v0building/v0builder.py` reads `_V0*` branches, `hCounter` and `nTrueInteractions` from the same directory.

## Skimmer pitfalls (fixed 2026-11-12, keep in mind)
- Any branch read in the skim event loop must be re-enabled after the `SetBranchStatus("_jet*", 0)`-style disables, because `GetEntry` does not fill disabled branches. `_jetPt` was missing, so the jet pT cut and the b-veto did nothing. Skims made before the fix (all of `/pnfs/.../K0sAnalysis/skim_results/run3/` as of 2026-11-12) effectively have **no b-veto**, and `_nimloth_nJets` is ~0 in them.
- `getYearFromFileName` must check sub-campaigns before their parent year (2022EE before 2022, 2023BPix before 2023). Data is matched on the era (Run2022E/F/G, Run2023D) and MC on the campaign name (`Run3Summer22EE`, `RunIII2024Summer24`, …). Unknown periods throw in `bTagSettings` instead of silently returning 0.
- Test the skimmer without real ntuples by compiling a harness that does `#define main skim_main` + `#include "skim_ztomumu.cc"`. Build it inside `cmsenv`: the system ROOT in `/usr/bin` lacks `Math/Vector4D.h`.

## Skim job submission (changed 2026-11-24)
- Each skim job (`skim_submit.py::make_job_commands`) runs the per-file skims and the `hadd` in a temporary directory on the worker's local disk. Only the merged `skimmed_N.root` is copied to `/pnfs`, then read back and compared (`cmp`), with up to 3 tries. The job runs under `set -e`: any failure exits non-zero and writes `###error###` to stderr.
- Why: writing per-file skims straight to `/pnfs` and having `hadd` read them back occasionally gave unreadable files (`R__unzip: error -5`, ~0.14% of files). With 100–300 files per job that broke ~27% of jobs in the 2026-11-23 round. The old scripts then deleted the per-file skims anyway and exited 0, leaving ~850-byte stub outputs.
- Keep `-n` small (5–10). Each file takes ~5 min, so 300-file jobs run more than a day.
- All jobs from one call go into one `condor_submit` (`condortools.submitCommandsAsCondorJobsBatched`). Condor still gives each job its own ClusterId, because `executable` changes between jobs.
- `skim_resubmit.py --since <date> [--dryrun] [-n 10]` resubmits failed jobs. It skips jobs still in `condor_q`, deletes the broken output, and splits large jobs into `skimmed_N-k.root`. It also marks the old executable `# skimjob-superseded` so that a rerun doesn't resubmit the same job again. Run it from `skimming/`. New-style jobs count as failed only if condor's `_log_` file shows a non-zero return value; stderr is not scanned, because a stage-out copy can crash and then succeed on retry. Old-style jobs (exit code always 0) still use the text checks in `tools/jobcheck.py`.
- `skim_submit.py` lists the input files when it is run. CRAB files that arrive later are not picked up. In the 2026-11-23 round, 1,768 such files (18 tasks, mostly `Muon1/Run2025G`) were skimmed afterwards as `skimmed_topup-<k>.root` in the existing output folders. Check that the CRAB tasks are complete before submitting.
- 2026-11-25: worker node `192.168.36.9` threw bus errors at ~09:32 (crashes in `cp`, in the skimmer, and missing `libGenVector.so`). 23 jobs failed; jobs that finished on that node before or despite the errors were checked and are intact.

## Known issues (open)
- No official BTV loose WPs exist yet for 2025 (any tagger), for DeepJet/PNet/RobustParT in 2024, or for UParT in 2022/2023. These are `PLACEHOLDER`s in `bTagSettingsTable`. The default veto taggers only use placeholders in 2025.
- In the `normmode=='range'` block of `mcvsdata_fill.py`, `file_total` is computed but no longer used (it was a leftover diagnostic).
- `mcvsdata_submit.py` hard-codes `CMSSW = '/user/jbierken/CMSSW_15_0_17'`.
- 2026 is wired in, but `reweighting/pileup/2026/Pileup_ratio.root` still needs the 2026 data pileup histogram (`reweighting/pileup/data/run3/dataPileupHistogram-2026pp_Golden-66000ub.root`, then `make_pileup.py -e 2026`). 2026 simulation cannot be reweighted until then.

## Conventions
- Tree names follow Tolkien naming: `laurelin` (K0S), `telperion` (Λ), `nimloth` (event-level), `celeborn` (muon branches).
- The repo is checked out under both `/ada_mnt/ada/user/jbierken/...` and `/user/jbierken/...`; these are the same files.
- `results/` holds committed example configs (Run 2 UL). Generated outputs (root/png) are not committed.
- There is no test suite or linter. `python3 -m py_compile` is the minimal check; real runs need the CMSSW environment and access to `/pnfs`.
