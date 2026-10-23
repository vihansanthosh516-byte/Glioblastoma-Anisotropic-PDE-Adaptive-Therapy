# GBM computational oncology project: status README

**This file records earlier GBM Tracks A/B/C results. GBM is now the main science-fair project; current goals are in START_HERE.md and science_fair_gbm/README.md.** Written 2026-10-17. Corrected the same day after a number audit against the ledger and `output/`. All numbers come from `vault/PAPER_FINDINGS_LEDGER.md` (entries cited as "ledger 2xx") or files in `output/`. I did not re-run any script for this file. Interim numbers are marked interim. README.md now gives the GBM project overview. This file preserves older platform results and plans; audit corrections dated 2026-10-19 do not reactivate those plans. Historical run-once/freeze compliance needs the named execution records.

## 1. What this project is

A glioblastoma (GBM) research platform with four tracks and 217 recursive `.py` files in `src/` (208 top-level, 159 numeric-prefixed top-level files, maximum prefix 145; counted 2026-10-19). This is a file inventory, not a runnable-pipeline count.

| Track | Scripts | Domain |
|---|---|---|
| A | `src/01-41` | Single-cell multi-omics, GRN, drug screen |
| B | `src/42-48` and later | Anisotropic PDE, cohort forecast, adaptive therapy |
| C | `src/49-65` and later | Digital twin, inverse estimation, MPC, RL |
| C-Extended | `src/neural_pde/`, `rl/`, `sensing/`, `hil/` | FNO, PPO, virtual sensors, pump interface |

Within the older GBM platform, the recurrence-target study records the main historical result (PREDICT-GBM, H-2 to H-4). Prospective rules require declarations before data access, single-use locked tests and negative-result reporting; blanket historical compliance is not established. H-2/H-3 are now spent/post hoc development sets under Amendment 12 (`analysis_plan_v1.md`, `vault/failed_hypotheses.md`).

## 2. What works (positive results)

| Result | Number | Source | Weakest point |
|---|---|---|---|
| All-recurrence coverage, H-2 former test, now spent | +0.081 (CI +0.065 to +0.097), n=101 | `output/predictgbm/final/test_results.json` | Secondary endpoint |
| All-recurrence coverage, H-3 former test, now spent (UCSF, new centre) | +0.085 (CI +0.059 to +0.115), n=59 | `output/predictgbm/h3_final/ucsf_test_results.json`, ledger 2bd | Secondary endpoint |
| Forecast on core target vs no-change | Dice +0.026, Holm p 0.0015 | `output/forecast_labels_core/results.json` | Target chosen post hoc. No DTI benefit detected; contrasts in `primary_tests_all` and `descriptive_shrank_or_same.tests`. Loses to no-change on shrinking cores (−0.039, n=50). |
| Age association with OS in adjusted CGGA model | adjusted HR 1.011, p 0.018 | `output/cgga_adjusted_age_mgmt.json:Y5_age_adjusted.terms.age` | Observational association, not predictive validation; about 1/3 the TCGA size |
| Two-tool disagreement proxy | floor 0.88 (middle volume tertile) vs 1.02 assumed | `output/seg_noise_floor.json` | Two automatic tools, not test-retest. Floor by tertile: 1.03 / 0.88 / 0.35. |
| Adaptive therapy, 61 patients (simulation) | MTD ends with more resistant cells in 54/61 (7 tied, 0 reversed), E_MAX ratio 1000 | `output/emax_sweep_summary.json` | Assumed kill rates. Adaptive ends with higher final tumour mass (VR 2.39). 10/61 progress earlier. Kill-scale sweep, ratios 100 to 8000: adaptive never selects more resistance (0/61 at all 7), the final-mass cost rises from 1.00 to 3.17, and the progression pattern is identical only at 500 to 2000. Controller-threshold sweep (8 cells, exploratory): still 0/61 more resistant and under 50% of the MTD dose in every cell; the progression pattern changes at OFF 0.50 and at ON 1.2 with OFF 0.80 or 0.90 (`vault/audit_2026-10-18/P4_threshold_sweep_results.csv`; ledger 2bn). |
| Elongation from anisotropy | +0.178, p 4e-19 (Wilcoxon, n=103) | `output/fractal_aniso_vs_iso.json` | Simulation only. Secondary endpoint. The primary endpoint (fractal dimension) went the other way (−0.018). |
| Released PREDICT-GBM U-Net rebuilt | Author-reported reproduction check, unverified without comparison artifact | `src/125_published_unet_inference.py` | The reported eight-patient comparison lacks a retained result citation; not independent validation |
| Process | Freeze files with hashes (`configs/h2_frozen.yaml`, `h3_frozen.yaml`, `h3v_frozen.yaml`), tags `h2-frozen`, `h3-frozen`, `h3v-frozen` | `git tag --list *frozen*`; configs | Tags exist. H-2/H-3 config hashes match saved score JSON; tag/hash presence alone does not prove one-shot compliance. |

## 3. What does not work (negative results, kept on purpose)

| Result | Number | Source |
|---|---|---|
| **Enhancing recurrence, H-2 former test (primary; now spent)** | +0.0088, CI −0.0115 to +0.0307, n=98 | `output/predictgbm/final/test_results.json` |
| **Enhancing recurrence, H-3 former test (primary; now spent)** | +0.0044, CI −0.0182 to +0.0369, n=49 | ledger 2bd |
| H-3 at 90% volume | −0.0053. Non-inferiority at −0.01 not shown | ledger 2bd |
| PDE beats persistence (H1), GBM-only n=105, aniso r=1 | +0.0130, CI −0.0038 to +0.0292 (all eligible, n=133: +0.0088) | `output/pde_manifest/results.json`, ledger 2y |
| DTI direction on forecast (aniso r=1 vs iso_same, GBM-only n=105, exact solver, sensitivity analysis) | −0.0004 | `output/pde_manifest_monotone/results.json`, ledger 2ak |
| Growth-only forecast vs no-change (LUMIERE) | worse (p 0.0003, HD-GLIO-AUTO; DeepBraTumIA p 0.004) | `output/lumiere_volume_forecast.json` |
| H-3v volume forecast | dev null. Locked LUMIERE worse than persistence | ledger 2ba |
| Stratified forecast | AUC 0.46, n=108 | `output/forecast_stratified/results.json:classifier.oof_auc` |
| Inverse estimation of D | D at a bound in 147/154. Growth rate rho at its lower bound in 124/154. D not identifiable | `output/inverse_est_metrics.json:summary` (in this fitted-volume setup) |
| Per-patient treatment effect | MGMT TTP p 0.60. Resistance gain 0 of 6 | `output/mgmt_resistance/results.json:clinical.logrank_p` (n=130); `output/plastic_resistance/results.json:verdict_real_test_seed59` (0/6 cells per AT50/AT80 arm, assumed-parameter simulations) |
| RL vs Stupp (simulated cohort64, n=20) | PPO wins 25%. On the 21 real fitted-growth-rate patients (real_test, assumed kill rates) PPO wins 20/21 on simulated day-90 volume (not a real-drug or confirmatory test), never on time to progression, and a one-line efficiency rule matches it. | `output/rl_kill_conditioned/evaluate.json`, ledger C1 |
| Track A leak-free classifier | acc 0.638 (gat_pca_clean, patient-level folds; random split 0.787 was leaking) | `output/cgat_leak_free.json` |

**Plain reading:** The two historical evaluations report positive secondary all-recurrence coverage differences near +0.08; both sets were later reused as post hoc development under Amendment 12. It does not beat the 15 mm plan on enhancing recurrence (null, twice). Growth forecasting and per-patient treatment effect are negative. Do not call any of this a "digital twin" (no validated forecast, no two-way data flow).

**Main limit on enhancing recurrence:** In 32 of 49 UCSF patients the 15 mm plan already covers at least 0.999 of enhancing tissue (`output/predictgbm/h3_final/ucsf_test_pairs.csv`, primary time1_n_enh>0 subset, time1_standard_enh>=0.999). There is little room to gain.

## 4. Historical H-2 to H-4 status snapshot

The following statuses describe the old platform. H-2/H-3 cannot become fresh locked tests. Current work follows `START_HERE.md`.

| Step | What | Status |
|---|---|---|
| H-2 | Train on TUM + MU, lock-test on LUMIERE + RHUH (n=101) | Historical result; sets now spent. Primary null, secondary positive |
| H-3 dev | Arms C1, C2, C3, U-Net compared on 137 TUM dev patients | Done. C1 +0.030, C3 +0.027, C2 +0.030 (dev, interim) |
| H-3 | C1 frozen, locked test on UCSF (n=59) | Historical result (ledger 2bd); set now spent. Primary null, secondary positive |
| H-4 step 1 | Plan rules on pooled dev (G10 rule picked) | Done, **post hoc** (ledger 2be). enh +0.0223, all +0.0869 |
| H-4 step 2 | Pooled retrain, 636 samples, 15 epochs, seed 20261009 | **Training done 2026-10-17** (attempt 1 crashed on disk, ledger 2bg). Model in `kaggle_run/out_h4_retrain/`. SHA-256 `b9846e1f…` checked. Training loss 0.7994 to 0.7164. No test run. |
| H-4 step 3 | New locked test set | **Historical desk-review state: no candidate confirmed eligible.** Desk check 2026-10-18: no candidate confirmed (see 5a.1). |

Retrain caveat: attempt 2 is a second launch of a run the plan calls one-shot. Ledger 2bg discloses it.

## 5. Historical proposed work (superseded as active next steps)

### 5a. Historical H-4 proposal (not the current project)
1. Find a new locked test set that meets A13.4: pre-op MRI, recurrence segmentation with a separate enhancing label, at least 40 patients, usable licence, never scored. Checked and rejected: CFB-GBM v2.0 (no pre-op MRI, one merged GTV), Burdenko-GBM-Progression (post-op baseline, restricted licence, no enhancing label). Desk-checked 2026-10-18 from public descriptions only (no download, no labels opened): BraTS 2025, UCSD-PTGBM, Duke longitudinal (Zenodo 20548528). None is confirmed eligible. Duke metadata now confirms the label map, automated labels and restricted/request-based access with CC BY-NC-ND 4.0 or a custom agreement; see ledger 2bm for the dated primary-source check. This does not establish eligible case count or approve a new test. UCSD-PTGBM is CC BY 4.0 but is described as post-treatment imaging. BraTS 2025 has unknown licence, per-patient pairing and overlap with sets already used. An unknown on any criterion blocks selection (A13.4). Others in `paper/research_prompt_v2.md` not checked.
2. If no set is found: apply A13.7. H-4 stays "post hoc dev" and the headline must say so.
3. If a set is found: write the amendment naming it before any label is read (A13.4), then create `configs/h4_frozen.yaml` and tag `h4-frozen` (A13.6), then run the one allowed test.
4. Ledger entry for the finished retrain (training only, no result): written 2026-10-18 as entry 2bh.
5. Commit `kaggle_run/out_h4_retrain/h4_retrain_train.json` and the log. Decide about the 5.6 MB `.pt` file.

### 5b. Research (prompts ready)
- Run `paper/research_prompt_v2.md` in a deep-research session. Topics: why enhancing recurrence is hard, post-op MRI, extra inputs (diffusion, perfusion, PET), small-data tricks, better loss, new test sets, growth forecast, treatment effect.
- `paper/literature_search_digital_twin.md` is a first pass from search snippets only. Verify each link before citing. Unverified: SAILOR access, any growth model beating persistence on real GBM, NASEM validation wording.

### 5c. Paper and review
- The paper is not written. The old figure-count claim is unverified without a defined inventory; `vault/PAPER_FIGURES.md` is a historical plan, not proof of a current RANO figure set.
- References are unchecked.
- Outside review (ChatGPT): check each claim against files, record in `paper/chatgpt_review_response.md`. Label interim numbers.
- Author-reported change history (not a fresh status check): `archive_old_gbm/REPO_MIND_MAP.md` row 34 was stale (said the UCSF test was not run). Updated 2026-10-18 with the locked-test result, together with the RL, adaptive-therapy and size lines the audit corrected.

### 5d. Repo hygiene
- `Requirements.txt` had merge-conflict markers; resolved 2026-10-18 by keeping the union of both sides. `pyproject.toml` does not list SALib, nibabel, dipy, scikit-learn or h5py, so add them there before dropping `Requirements.txt`.
- Script numbers reused (15, 16, 22, 50-56). `archive_old_gbm/REPO_MIND_MAP.md` section 13b lists the files behind each number. No file was renamed.
- Untracked run kits and caches (`kaggle_run/*`, `output/pde_manifest/cache_*`). Historical storage measurements lack a retained filesystem manifest and are not current facts. Re-measure before deciding what to ignore.
- Ledger 2bo records author-reported historical split scratch runs totaling 101/101, not one invocation. Current static inventory is 18 test files and 117 test functions. The prior audit reports 108 passed/13 warnings in a Git-backed scratch clone (`paper/codex_audit_report.md`, section E); the latest audit-fix scratch run passed 119 test cases with 13 warnings in 27.74 s (`paper/codex_audit2/fix_progress.md`).
- Run `graphify update .` after code changes.

### 5e. Long-term options (all unproven)
- Part 2 (growth forecast): needs 3+ time points per patient. LUMIERE and MU have few. Hierarchical or Bayesian calibration is untested.
- Part 3 (treatment effect): needs per-patient dose and TMZ dates. Candidates: CFB-GBM v2.0 dose maps, Burdenko (restricted). Test-first, simulation-second.

## 6. Historical old-work rules (current rules: START_HERE.md)
- Every number cites a file in `output/` or the ledger.
- Name the weakest point next to each result.
- Words to use: "simulation", "assumed", "post hoc", "exploratory", "interim".
- Never hand-edit output JSON. Fix code and re-run.
- Commit format: `Script N: brief description`. Check `git status --short` before every commit. Push after every commit.

## 7. Historical old-work index (current index: START_HERE.md)
| File | Use |
|---|---|
| `analysis_plan_v1.md` | Amendments through A21; prospective rules, with chronology requiring per-run evidence |
| `vault/PAPER_FINDINGS_LEDGER.md` | Every result with source file |
| `vault/failed_hypotheses.md` | Negatives |
| `archive_old_gbm/REPO_MIND_MAP.md` | Full inventory with status tags |
| `configs/h*_frozen.yaml` | Freeze files with hashes |
| `paper/research_prompt_v2.md` | Next research prompt |
| `paper/literature_search_digital_twin.md` | First-pass literature notes |
