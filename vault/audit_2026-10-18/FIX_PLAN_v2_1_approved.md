# Fix plan v2.1 (revised after review rounds 1 and 2): APPROVED TO START

Author: analysis assistant. Date: 2026-10-18. Status: reviewed twice by PLAN_REVIEWER_BIOSTAT (round 2: approve with changes, may start; changes below applied). Nothing below has been started at the time of writing.

## What changed from v1 and why
R1 (F1) Row C1 annotation is split into three separately tagged sub-entries.
R2 (F2) P4 is now a full crossed 4 x 2 design with a pre-specified interpretation rule. A code reading also changed the design: script 44 restarts dosing when `mass > THRESHOLD_ON * baseline OR mass > 0.9 * peak` (hard-coded 0.9), and stops dosing when `mass < THRESHOLD_OFF * baseline OR mass < THRESHOLD_OFF * peak`. Because peak >= baseline, the 0.9*peak clause fires first in most cases, so THRESHOLD_ON is expected to be near-inert. The ON axis is kept as a declared expectation (identical results to ON=1.0) and not as a lever. Script 44 stays unchanged.
R3 (F3) P3 gets an operational exclusion list and an 'unknown' rule.
R4 (F4) Every ledger addition starts with a fixed header: [DATE | STEP-ID | TYPE: measured / assumed / exploratory / correction].
R5 (F5) The P4 declaration is saved as a timestamped artifact, with its SHA-256, before any P4 run; its artifact id is cited in P4 outputs.
R6 (F6) P5 has a regression stop condition.
R7 (F7) Ledger additions state that the H-2 and H-3 test sets were reused as dev under Amendment 12 and are spent for any future confirmatory use.
R8 (Q5) The retrain's training-only metrics are labelled descriptive; no confirmatory claim may be made from them. P0 verifies backup checksums.

## Guardrails (apply to every step)
G1. No result is relabelled. No null is reworded as a positive. Nulls stay where they are.
G2. No locked-test run. No label of any candidate H-4 dataset is read or downloaded.
G3. Every repo edit uses edit_file after saving the exact original bytes as an artifact. The ledger is append-only: corrections are dated annotations with the R4 header; original text is never deleted.
G4. No commit, push, rename or deletion unless the user asks.
G5. Everything new is labelled exploratory or post hoc.
G6. The P4 declaration, including the interpretation rule, is saved as an artifact before any P4 run. All cells are reported; none dropped.
G7. If a step fails twice, stop and report instead of improvising.

## Phases
P0 Baseline. Record `git status --short -unormal` and SHA-256 of every file to be edited. Save backups as artifacts and confirm each artifact checksum equals the local SHA-256 (backups live outside the repo).
P1 Ledger (vault/PAPER_FINDINGS_LEDGER.md), append-only, every addition with the R4 header.
 1.1 B1 annotation: 54/61 verified in output/emax_sweep_summary.json ratio 1000 (7 ties, 0 reversed); original text kept.
 1.2 Entry 2w annotation: PDE CI superseded by 2y.
 1.3 C1 annotation as three sub-entries: (a) simulated PPO 25% on cohort64, n=20, TYPE assumed (simulated kill rates), not confirmatory; (b) real_test 20/21 on day-90 volume, TYPE assumed (same assumed kill rates for all 21 patients), measured growth rates, not a real-drug test and not confirmatory; (c) no time-to-progression gain, and a one-line efficiency rule matches the policy.
 1.4 New entries: retrain finished (training only, no test data read, SHA-256, loss; descriptive only); kill-rate grid (E_MAX 100 to 8000); P4 threshold results; README audit summary; literature check (three papers, limits); note on script 44 (stale comment says 75%, default is 0.80; hard-coded 0.9*peak clause).
 1.5 Every new entry says: H-2 and H-3 test sets were reused as dev (Amendment 12) and are spent for any future confirmatory analysis.
P2 Repo hygiene (no deletions, no renames).
 2.1 Requirements.txt: resolve merge markers. pyproject.toml wins; packages present in either side but absent from pyproject.toml are listed for the user, not dropped.
 2.2 Tests: run pytest on a COPY of src/, tests/, configs/ and small output files in the workspace, cache provider off. Update the CLAUDE.md test-count line in place (overwrite, one line only) to the observed collected / passed / failed counts, and record before and after in a dated ledger entry with the R4 header (round-2 F10). Failures are reported, not fixed.
 2.3 REPO_MIND_MAP: update row 34 to the H-3 UCSF locked test result (ledger 2bd); add a table of colliding script numbers and their files instead of renaming.
 2.4 Untracked kits: re-measure sizes; write a proposed .gitignore patch as a separate artifact; do not edit .gitignore.
P3 H-4 blocker desk review. Candidates: BraTS 2025, UCSD-PTGBM, Duke longitudinal (Zenodo 20548528). Check each against A13.4 (pre-op MRI; recurrence segmentation with a separate enhancing label; at least 40 patients; research-use licence; never scored by this project).
 Allowed sources: dataset description pages, licence terms, cohort and annotation documentation.
 Excluded: any page, table or section reporting model performance, benchmark results or leaderboard rankings for any candidate dataset. If encountered incidentally it is noted as excluded and not used.
 Output: yes / no / unknown per criterion with a source per cell. An 'unknown' is never inferred; an unknown on any criterion blocks selection until the user resolves it with the data provider. I do not download anything and I do not choose a dataset (A13.4 requires the user's amendment).
P4 Controller-threshold sweep (exploratory, descriptive, no p-values). E_MAX ratio fixed at 1000. Crossed design: THRESHOLD_OFF in {0.50, 0.65, 0.80, 0.90} x THRESHOLD_ON in {1.0, 1.2}. The cell (0.80, 1.0) is the default, already reproduced exactly. Seven new runs in parallel (about 30 to 55 min; memory monitored). The simulation is deterministic (ratio-1000 re-run differed by 0.0), so no replicates are needed; spread is reported across the 61 patients.
 Reported per cell: MTD more resistant / tied / adaptive more resistant; median, IQR and range of per-patient VR (adaptive / MTD final mass) and of drug % of MTD; TTP adaptive earlier / later / same; the rho > 0.02 split (earlier counts above and below); fraction of time drug is on; both-censored count.
 Pre-specified interpretation rule (fixed before any run):
  C-A 'adaptive never selects more resistance' is called robust only if adaptive-more-resistant = 0 in all 8 cells; otherwise sensitive, naming the cells.
  C-B 'dose-sparing' is called robust only if mean drug % of MTD is below 50 in all 8 cells (default is 32.4); otherwise sensitive. Rationale: a halving of the MTD dose is the smallest saving I would call substantial; this is a convention, not a clinical standard, and passing it is a modest form of 'robust'. (Round 2 asked to tighten the first draft cutoff of 80; 50 is fixed here before any P4 run.)
  C-C 'cost' (VR) is described with numbers only; no robust/sensitive label.
  C-D 'progression pattern' (10 of 10 above rho 0.02 earlier, 0 of 51 below) is called robust only if both counts equal the default in all 8 cells; otherwise sensitive, naming the cells.
  Deviation reporting (round-2 F9): for every cell, the write-up gives the difference from the default cell (counts, VR median and IQR, drug %), whether or not the cell is labelled sensitive. Labels are not to be read without these magnitudes.
  C-E 'ON is inert' (declared expectation): for each OFF value, ON = 1.2 and ON = 1.0 give identical per-patient metrics. Any difference is reported with the patients affected.
 The declaration (settings, metrics, interpretation rule, fixed E_MAX ratio) is saved as an artifact and its SHA-256 recorded before the first run.
P5 Verification. Re-run the 35-row audit against the final README and ledger. List every repo file changed versus the P0 baseline and confirm only intended files changed. Produce one final summary artifact.
 Regression stop condition: if any previously matching audit row fails, halt further edits and report to the user.

## Stop conditions
Stop and ask the user if: any step would touch a test set or labels; pytest needs the user's data in a way that risks overwriting outputs; the Requirements.txt conflict cannot be resolved by the stated rule; any verification row fails; any P4 run does not complete cleanly twice.

## Not in scope
Running the H-4 locked test; choosing the H-4 dataset; fixing failing tests; renaming scripts; committing or pushing; editing script 44.

## What I want reviewed (round 2)
1. Do R1 to R8 resolve each round-1 finding? Name any that do not.
2. Is the P4 design and the pre-specified interpretation rule adequate? Are the thresholds in C-B declared reasonably (no look at results)?
3. Any new blocker introduced by v2?
