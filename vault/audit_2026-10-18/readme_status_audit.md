# README_STATUS.md number audit

Read-only audit of `README_STATUS.md` (written 2026-10-17) against `vault/PAPER_FINDINGS_LEDGER.md` and the files in `output/`, `configs/`, `kaggle_run/`, `src/`, `tests/`. No script was re-run and nothing in the repo was changed. Values below were read from the JSON/CSV on disk or recomputed from them (UCSF 32/49 from `ucsf_test_pairs.csv`; the retrain SHA-256 from the `.pt` file; test, script and size counts from the file tree).

## Result

35 claims checked: **25 match**, 2 minor, 7 need a fix or a re-check, 1 could not be checked from the repo.

Every headline number in section 2 (positive results) matches its file, except that the adaptive-therapy count is cited to the wrong file. The problems are concentrated in section 3 (negative results).

## Items that need action (most serious first)

| # | Status | Claim as written | Value found in files | Note |
|---|---|---|---|---|
| 20 | WRONG | RL vs Stupp on real patients: 25% win | 25% is PPO on cohort64 (simulated, n=20) in output/rl_kill_conditioned/evaluate.json (25.0%). On real_test (n=21) the same PPO arm wins 95.2% (rl_kill_conditioned) / 95.2% (rl_equal_budget) | Wrong patient set and wrong file. On the 21 real fitted-growth-rate patients PPO beats Stupp on day-90 volume in 20/21, using assumed kill rates. The ledger caveat is that it never wins on time to progression and a one-line rule matches it. |
| 13 | DISCREPANCY | PDE beats persistence (H1): +0.0130, CI -0.0043 to +0.0302 | GBM-only (n=105) aniso r=1: +0.0130, CI -0.0038 to +0.0292. All-eligible (n=133): +0.0088, CI -0.0060 to +0.0236 | Mean matches. CI does not: README repeats ledger 2w (older bootstrap); ledger 2y re-ran at 10,000 draws and the file now says -0.0038 to +0.0292. Correct CI is -0.0038 to +0.0292. README also does not say the figure is GBM-only (n=105). |
| 6 | SOURCE ERROR | Adaptive therapy, 61 patients: 54/61 less resistance | Not in cited file (it has n=61, 10/61 earlier, 39/61 censored, VR 2.39). 54 is in output/emax_sweep_summary.json ratio 1000: n_mtd_more_resistant=54, ties=7, reverse=0 | Number is correct; the cited file is wrong. 7 of 61 are ties, so "54/61 less resistance" is 54 less, 7 equal, 0 more. Ledger B1 still lists 54/61 as "vault only"; entry 2n later confirmed it. |
| 18 | MISLABEL | Inverse estimation of D: 124/154 at bound; D not identifiable | n_rho_at_lower_bound=124/154; n_D_at_bound=147/154; D_identifiable_from_volumes=False | 124/154 is growth rate rho at its lower bound, not D. D is at a bound in 147/154. Conclusion (D not identifiable) is right. |
| 16 | CITATION | H-3v: dev null, locked LUMIERE worse than persistence | Content matches ledger 2ba (locked: worse by +0.37 and +0.33 log error). 2ay is the H-2 sensitivity entry | Cite ledger 2ba (and 2bc if RANO is meant), not the range 2ay-2bc. |
| 34 | PATH | failed_hypotheses.md at repo root | File is at vault/failed_hypotheses.md; not present at the root | Fix the path in sections 1 and 7. |
| 27 | MINOR | "190+ numbered scripts"; track ranges 01-41, 42-48, 49-65 | src/ has 198 top-level .py files (207 with subfolders) but only 150 begin with a number (135 distinct numbers, max 136). Numbers 1-41: 41 present; 42-48: 7; 49-65: 17 | "190+ scripts" is right, "190+ numbered scripts" is not (150). MIND_MAP says 196 script rows. |
| 23 | MINOR | H-3 dev (137 TUM): C1 +0.030, C3 +0.027, C2 +0.031 | C1 0.02983; C3 0.02743; C2 0.03048; n=137 | C2 is +0.0305, which rounds to +0.030, not +0.031 (likely rounded twice). |
| 31 | UNVERIFIED | About 3.8 GB untracked run kits and caches | Not reproduced. Untracked, not ignored (git status): 2.0 GB, of which kaggle_run 1.73 GB; pde_manifest cache_* about 2 MB. All of kaggle_run on disk incl. gitignored dirs: 4.87 GB | Size changes as files are added; re-measure. 3.8 GB matches neither slice. |
| 35 | UNVERIFIABLE | BraTS 2025, UCSD-PTGBM, Duke longitudinal (Zenodo 20548528) status | Not checked: external datasets, no repo file to compare against | Would need a web lookup. |

## What I would change in the README

1. **Row 20 (RL vs Stupp).** Replace "RL vs Stupp on real patients: 25% win" with what the files show. Either report PPO on the simulated `cohort64` set (25%, n=20, `output/rl_kill_conditioned/evaluate.json`), or report the real-patient set honestly: PPO wins 20/21 on day-90 volume under assumed kill rates, never wins on time to progression, and a one-line efficiency rule matches it. As written, a positive result on the real-patient set is listed as a negative.
2. **Row 13 (PDE H1).** Change the CI to -0.0038 to +0.0292 and add "GBM-only, n=105". The conclusion (every CI includes 0) is unchanged.
3. **Row 6 (adaptive 54/61).** Cite `output/emax_sweep_summary.json` (ratio 1000), and write "54 less, 7 equal, 0 more resistant". Consider adding the cost already in the ledger: final tumour mass is higher under adaptive (VR 2.39) and 10/61 progress earlier.
4. **Row 18 (inverse estimation).** 124/154 is rho at its lower bound; D is at a bound in 147/154.
5. **Rows 16, 23, 27, 34.** Cite ledger 2ba for H-3v; C2 is +0.030; say "about 150 numbered scripts" (or "190+ scripts"); `failed_hypotheses.md` is in `vault/`.

## Things that match but need a qualifier

- Row 5: the 0.88 noise floor is the middle volume tertile; the smallest tertile is 1.03, the same as the assumed 1.02.
- Row 7: elongation is a secondary endpoint; the primary fractal-dimension endpoint moved the opposite way (-0.018).
- Row 15: p 0.0003 is the HD-GLIO-AUTO segmenter only; DeepBraTumIA gives p 0.004.
- Row 21: name the variant (`gat_pca_clean`); the other two leak-free-labelled variants reach 0.742 and 0.758 but still carry leaks (ledger 2o).
- Row 3: the core forecast loses to no-change on the 50 shrinking or stable cores (-0.039).

## Not checked

- Row 35: dataset status claims (BraTS 2025, UCSD-PTGBM, Duke/Zenodo 20548528) are external.
- Row 31: the 3.8 GB figure is not reproduced (2.0 GB untracked by git status; 4.87 GB for all of `kaggle_run` including gitignored dirs).
- Pytest was not run, so the test count is a grep count only.
- The "H-4 step 3 blocked: no eligible dataset" statement is a judgment about the outside world, not a number.
- I did not audit the ledger itself. Where the ledger and a JSON file disagreed (row 13), I report both and say which one is current.

## Side finding for task 3 (kill-rate sensitivity)

`output/emax_sweep_summary.json` already holds a kill-scale sweep at E_MAX ratio 500 / 1000 / 2000 for the 61-patient cohort (ledger 2n). Across it, "MTD more resistant" is 57 / 54 / 51 of 61, and final-mass ratio rises 1.73 / 2.39 / 2.89. A finer or wider grid would extend that sweep; ratios outside 500-2000 were not run.
