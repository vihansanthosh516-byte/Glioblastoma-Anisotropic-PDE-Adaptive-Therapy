> **HISTORICAL BACKGROUND (corrected 2026-10-19):** Dated old-GBM audit record, not current instructions or authorization to run analyses. Current goal/status: START_HERE.md. Simulation settings are assumptions, not measured treatment effects. Original approval, chronology and process claims are author-reported unless supported by retained logs.

# Audit and fix run, 2026-10-18: index

This folder holds the files cited by ledger entries 2bi, 2bj, 2bk, 2bm, 2bn and 2bo in `vault/PAPER_FINDINGS_LEDGER.md`, plus the plan and review record for the run. Nothing here changes a result: nulls stay null, and the new simulations are exploratory with an assumed kill-rate scale.

## Author-reported historical changes (exact commit membership unverified here)
| File | State |
|---|---|
| `README_STATUS.md` | new in git (was untracked); corrected after a 35-claim number audit |
| `REPO_MIND_MAP.md` | new in git (was untracked); stale rows and claims corrected, section 13b added |
| `vault/PAPER_FINDINGS_LEDGER.md` | additions only: 3 dated annotations and entries 2bh to 2bo |
| `CLAUDE.md` | line 79 (test count and module list) and line 123 (Requirements note) |
| `Requirements.txt` | merge conflict resolved by keeping the union of both sides |
| `vault/audit_2026-10-18/` | this folder |
Pre-edit copies of the two files that had no git history are in `originals/`. The other three files are recoverable from git history. `.gitignore` was not edited.

Originals notice (2026-10-19 audit correction): `originals/README_STATUS_original_2026-10-17.md` and `originals/REPO_MIND_MAP_original_2026-10-14.md` are unchanged dated snapshots (ARCHIVE in `START_HERE.md` §6), kept byte-identical on purpose and not current guidance. Where they differ from `README_STATUS.md`, `README_STATUS_change_log.md` and `vault/PAPER_FINDINGS_LEDGER.md` win. Known superseded items in the README snapshot: the 190+/150 script counts (current: 202 `.py`, 154 numbered); `failed_hypotheses.md` path (now `vault/failed_hypotheses.md`); adaptive-therapy "54/61 less resistance" (corrected: MTD ends with more resistant cells in 54/61, 7 ties, 0 reversed, simulation with assumed kill rate, `output/emax_sweep_summary.json`); inverse-estimation 124/154 (corrected: D at bound 147/154, rho at lower bound 124/154); "25% win" (simulated cohort64, n=20); H-3 C2 +0.031 (displayed +0.030, dev, interim); the PDE H1 interval (earlier 2,000-draw bootstrap, superseded by 10,000 draws); the 3.8 GB size (unreproduced); pytest status; the H-4 candidate list (2026-10-18 desk review confirmed none eligible); and the "locked" labels on the H-2/H-3 sets, which are spent/post hoc development data after Amendment 12. Values were cross-checked only against the later status files and change log, not re-extracted from outputs.

Mind-map snapshot notice (2026-10-19): `originals/REPO_MIND_MAP_original_2026-10-14.md` is a dated 2026-10-14 snapshot, byte-identical on purpose, with no banner of its own; the maintained copy is `archive_old_gbm/REPO_MIND_MAP.md`. Its file counts, storage sizes, git state, test count (16 files/101 functions, "not run") and `Requirements.txt` merge-marker note are stale; markers were resolved and current counts are in `README_STATUS.md`/`START_HERE.md`. Its per-script result rows were compiled from vault notes and a few JSONs without rerunning scripts, so treat each as an unverified index entry. Known corrections: adaptive-therapy source is `output/emax_sweep_summary.json` (MTD more resistant in 54/61, 7 ties, 0 reversed, simulation with assumed kill rate, adaptive final tumour mass higher), not `adaptive_cohort_summary.json`; the .787 random-cell C-GAT result is leak-prone and the leak-free value is .638 +/- .133 (`gat_pca_clean`); core forecast +.026 is post hoc (shrinking cores -.039, n=50); MPC dose sparing 68.9% came with higher variance (the -31.8% is an increase); FNO >1000x speed is withdrawn; rho>.02 10/10 vs 0/51 is post hoc exploratory. Rows for H3 exact-solver, LUMIERE, script 60 and spatial ablation were not re-extracted from JSON here and stay unverified.

## Files here
| File | What it is | Ledger |
|---|---|---|
| `readme_status_audit.md`, `.csv` | the 35-claim audit of the old README_STATUS.md against the ledger and `output/` | 2bj |
| `README_STATUS_change_log.md` | old and new text for each README edit made after the audit | 2bj |
| `kill_rate_grid_*.md/.csv/.png` | E_MAX ratio 100 to 8000, 61-patient adaptive-therapy simulation: declaration, report, per-ratio table, figure | 2bi |
| `P4_threshold_sweep_*` | controller-threshold sweep (THRESHOLD_OFF x THRESHOLD_ON, 8 cells): declaration written before the runs, report, table, figure | 2bn |
| `P4_ON_inertness_check.csv` | per-OFF count of patients whose results differ between ON 1.0 and ON 1.2 | 2bn |
| `raw_metrics/` | per-patient metrics and cohort summaries from every new run, unedited (JSON, LFS-tracked: present locally, not in this commit; see the LFS note below) | 2bi, 2bn |
| `literature_check.md` | three papers read in full, and what they imply for this project | 2bk |
| `H4_dataset_desk_review.md` | A13.4 eligibility desk check of three candidate locked-test datasets (metadata only; nothing downloaded; none chosen) | 2bm |
| `FIX_PLAN_v1_initial.md`, `FIX_PLAN_v2_1_approved.md`, `plan_review_record.md` | the plan, its two review rounds and the revisions they caused | |
| `proposed_gitignore_additions.txt` | a proposed `.gitignore` patch, NOT applied | |
Ledger entry 2bo (pytest) is self-contained: 101 of 101 tests passed, but in a copy with partial data and not in one single run.

## Reproducing a simulation run
`src/44_adaptive_therapy.py` is unchanged. From the repo root, with output outside the repo:
`PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 ADAPTIVE_E_MAX_RATIO=<ratio> ADAPTIVE_THRESHOLD_OFF=<x> ADAPTIVE_THRESHOLD_ON=<y> ADAPTIVE_OUT_DIR=<dir> python src/44_adaptive_therapy.py`
Author-reported timing: about 27 minutes alone or 38 minutes with seven runs on 8 cores; dated timing log/hardware environment missing. The default setting (ratio 1000, OFF 0.80, ON 1.0) reportedly reproduced the committed `output/emax_sweep/r1000` exactly (max absolute difference 0.0).

## Things to know
- Author-reported LFS history (exact failed-add transcript and historical commit SHA missing; current presence does not prove prior membership): `*.png` and `*.json` are Git LFS files in this repo (`.gitattributes`). `git add` of those files failed in the environment used for this run (the LFS filter could not run), and the repo's LFS rules were deliberately not bypassed. So `kill_rate_grid.png`, `P4_threshold_sweep.png` and everything under `raw_metrics/` are in the working tree but NOT in the commit. From a machine where git-lfs works, run `git add vault/audit_2026-10-18` and commit; until then the two report figures show as broken images on GitHub and `raw_metrics/` is missing there.
- The two plan reviews were run as single-shot model calls using the reviewer profile's system prompt, because the delegation tool was switched off; they were not run by a separate agent. See `plan_review_record.md`.
- Threshold sweep: the progression pattern is sensitive to the thresholds, and the declared expectation that THRESHOLD_ON is inert was not met (it changes results for 12 to 26 of 61 patients, all faster growers).
- Author-reported date/clock explanation (independent clock record missing): entries are dated 2026-10-18 (session date); the machine clock read late 2026-10-17 local time, so the 2026-10-17 audit and these entries are the same working session.
- Open decisions for the author: apply the proposed `.gitignore` patch; whether to keep the large untracked run kits out of git; the H-4 dataset choice (an A13.4 amendment, before any label is read).

Declaration preservation note: `P4_threshold_sweep_declaration.md` and `kill_rate_grid_declaration.md` are immutable historical declarations; their lack of an inline BACKGROUND banner is resolved here rather than changing their bytes. Settings are model assumptions; before-run timing and reported default reproduction need independent timestamp/log evidence. P4 expected ON inertness was falsified; see the report and ledger 2bn. Current retained P4 declaration SHA-256 is 66d7b40b9cb76c6203a8ea8b1cad89cc67537951533796a7fa2d2568cde713d8; it differs from the originally reported artifact-version hash. This verifies local bytes only, not pre-run identity.
