# Fix plan v1 (before review): finish everything still open on the GBM project

Author: analysis assistant. Date: 2026-10-18. Status: DRAFT for review by PLAN_REVIEWER_BIOSTAT. Nothing below has been started.

## Context (verified earlier in this project)
- README_STATUS.md was corrected after a 35-claim audit (25 match, 2 minor, 7 fixed or flagged, 1 unverifiable) and written into the repo with the user's approval. Original bytes are saved as an artifact.
- A kill-rate sensitivity grid (E_MAX ratio 100 to 8000, 61-patient adaptive-therapy simulation) was run exploratorily. Ratio-1000 re-run reproduced the committed output exactly.
- Open items from the audit: ledger entries that still carry the old source or stale CI; Requirements.txt merge-conflict markers; CLAUDE.md says 56 tests, tests/ has 101 `def test_`; REPO_MIND_MAP row 34 stale; script numbers reused (15, 16, 22, 50-56); untracked run kits (2.0 GB untracked, 4.87 GB in kaggle_run); H-4 step 3 blocked (no eligible locked test set); controller thresholds (ADAPTIVE_THRESHOLD_OFF 0.80 / ON 1.0) never varied in the sweep.
- Study rules (analysis_plan_v1.md): each test declared before data is read; a locked test set is used once; negatives are reported. H-2 and H-3 primaries (enhancing recurrence) are null; all-recurrence is a secondary endpoint and positive.

## Guardrails (apply to every step)
G1. No result is relabelled. No null is reworded as a positive. Nulls stay where they are.
G2. No locked-test run. No label of any candidate H-4 dataset is read or downloaded. Dataset review is a desk check of public metadata only.
G3. Every repo edit uses edit_file after saving the exact original bytes as an artifact. The ledger is append-only: corrections are dated annotations; original text is never deleted.
G4. No commit, push, rename or deletion unless the user asks.
G5. Everything new is labelled exploratory or post hoc.
G6. The threshold sweep (Phase 4) is declared in writing before any run. All settings are reported; none dropped.
G7. If a step fails twice, stop and report instead of improvising.

## Phases
P0 Baseline. Record `git status --short -unormal` and SHA-256 of every file to be edited; save backups as artifacts.
P1 Ledger (vault/PAPER_FINDINGS_LEDGER.md), append-only.
 1.1 Row B1: annotate that 54/61 is verified in output/emax_sweep_summary.json, ratio 1000 (7 ties, 0 reversed); keep original text.
 1.2 Entry 2w: annotate that its PDE CI is superseded by 2y.
 1.3 Row C1: state both numbers explicitly: PPO 25% on cohort64 (simulated, n=20); 20/21 on real_test day-90 volume with assumed kill rates; no TTP gain.
 1.4 New entries: retrain finished (training only, no test data read; SHA-256, loss); kill-rate grid; README audit summary; literature check (three papers, limits).
P2 Repo hygiene (no deletions, no renames).
 2.1 Requirements.txt: resolve merge markers. Rule: pyproject.toml wins; any package present in either side but absent from pyproject.toml is listed for the user, not silently dropped.
 2.2 Tests: run pytest on a COPY of src/, tests/, configs/ and small output files in the workspace (not in the repo), with the cache provider off, because tests might write into output/. Update the CLAUDE.md test counts to what is observed (collected, passed, failed). Failures are reported, not fixed, unless trivial and the user agrees.
 2.3 REPO_MIND_MAP: update row 34 to the H-3 UCSF locked test result (ledger 2bd); add a table of the colliding script numbers and their files instead of renaming.
 2.4 Untracked kits: re-measure sizes; write a proposed .gitignore patch as a separate artifact. Do not edit .gitignore.
P3 H-4 blocker desk review. For BraTS 2025, UCSD-PTGBM, Duke longitudinal (Zenodo 20548528): check each against A13.4 using public metadata only (pre-op MRI; recurrence segmentation with separate enhancing label; at least 40 patients; research-use licence; never scored by this project). Output a table of yes / no / unknown with a source for each cell. Do not choose a dataset: A13.4 requires the user's amendment, chosen on eligibility and logistics only.
P4 Controller-threshold sweep (exploratory). E_MAX ratio fixed at 1000. Declare first. New settings (THRESHOLD_OFF, THRESHOLD_ON): (0.50, 1.0), (0.65, 1.0), (0.90, 1.0), (0.80, 0.90). Default (0.80, 1.0) already reproduced. Same metrics as script 92. No p-values; counts and means only. Four parallel runs, about 30 min.
P5 Verification. Re-run the 35-row audit against the final README and ledger. List every repo file that changed versus the P0 baseline and confirm only intended files changed. Produce one final summary artifact.

## Stop conditions
Stop and ask the user if: any step would touch a test set or labels; pytest needs the user's data in a way that risks overwriting outputs; the Requirements.txt conflict cannot be resolved by the stated rule; any verification row fails.

## What I want reviewed
1. Does anything in this plan compromise the study's pre-declaration, locked-test or negative-reporting rules?
2. Is the P4 design (four settings, no inference, declared first) adequate for an exploratory sensitivity analysis, or should it change?
3. Is the P3 desk review safe with respect to A13.4 (no outcome information used to choose a dataset)?
4. Is adding ledger entries and annotations acceptable practice for a vault that is called "the single source for every number in the paper"?
5. What is missing, or what would you cut?
