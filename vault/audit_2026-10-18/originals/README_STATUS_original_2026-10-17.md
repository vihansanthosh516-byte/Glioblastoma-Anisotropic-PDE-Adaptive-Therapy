# GBM computational oncology project: status README

Written 2026-10-17. All numbers come from `vault/PAPER_FINDINGS_LEDGER.md` (entries cited as "ledger 2xx") or files in `output/`. I did not re-run any script for this file. Interim numbers are marked interim. The old `README.md` is unchanged.

## 1. What this project is

A glioblastoma (GBM) research platform with four tracks and 190+ numbered scripts (see `REPO_MIND_MAP.md` for the full inventory).

| Track | Scripts | Domain |
|---|---|---|
| A | `src/01-41` | Single-cell multi-omics, GRN, drug screen |
| B | `src/42-48` and later | Anisotropic PDE, cohort forecast, adaptive therapy |
| C | `src/49-65` and later | Digital twin, inverse estimation, MPC, RL |
| C-Extended | `src/neural_pde/`, `rl/`, `sensing/`, `hil/` | FNO, PPO, virtual sensors, pump interface |

The main science result now sits in the recurrence-target study (PREDICT-GBM, H-2 to H-4). Rules for the whole study: each test is declared before data is read, a locked test set is used once, negative results are reported (`analysis_plan_v1.md`, `failed_hypotheses.md`).

## 2. What works (positive results)

| Result | Number | Source | Weakest point |
|---|---|---|---|
| All-recurrence coverage, H-2 locked test | +0.081 (CI +0.065 to +0.097), n=101 | `output/predictgbm/final/test_results.json` | Secondary endpoint |
| All-recurrence coverage, H-3 locked test (UCSF, new centre) | +0.085 (CI +0.059 to +0.115), n=59 | `output/predictgbm/h3_final/ucsf_test_results.json`, ledger 2bd | Secondary endpoint |
| Forecast on core target vs no-change | Dice +0.026, Holm p 0.0015 | `output/forecast_labels_core/results.json` | Target chosen post hoc. DTI adds nothing. |
| Age predicts survival (TCGA to CGGA) | adjusted HR 1.011, p 0.018 | `output/cgga_adjusted_age_mgmt.json` | About 1/3 the TCGA size |
| Measured segmentation noise | floor 0.88 vs 1.02 assumed | `output/seg_noise_floor.json` | Two automatic tools, not test-retest |
| Adaptive therapy, 61 patients | 54/61 less resistance | `output/adaptive_cohort_summary.json` | Simulation, assumed kill rates |
| Elongation from anisotropy | +0.178, p 4e-19 | `output/fractal_aniso_vs_iso.json` | Simulation only |
| Released PREDICT-GBM U-Net rebuilt | reproduces released predictions | script 125 | Checked on 8 dev patients |
| Process | Freeze files with hashes (`configs/h2_frozen.yaml`, `h3_frozen.yaml`, `h3v_frozen.yaml`), tags `h2-frozen`, `h3-frozen`, `h3v-frozen` | git | — |

## 3. What does not work (negative results, kept on purpose)

| Result | Number | Source |
|---|---|---|
| **Enhancing recurrence, H-2 locked (primary)** | +0.0088, CI −0.0115 to +0.0307, n=98 | `output/predictgbm/final/test_results.json` |
| **Enhancing recurrence, H-3 locked (primary)** | +0.0044, CI −0.0182 to +0.0369, n=49 | ledger 2bd |
| H-3 at 90% volume | −0.0053. Non-inferiority at −0.01 not shown | ledger 2bd |
| PDE beats persistence (H1) | +0.0130, CI −0.0043 to +0.0302 | `output/pde_manifest/results.json` |
| DTI direction on forecast | −0.0004 | `output/pde_manifest_monotone/results.json` |
| Growth-only forecast vs no-change (LUMIERE) | worse (p 0.0003) | `output/lumiere_volume_forecast.json` |
| H-3v volume forecast | dev null. Locked LUMIERE worse than persistence | ledger 2ay-2bc |
| Stratified forecast | AUC 0.46 | `output/forecast_stratified/` |
| Inverse estimation of D | 124/154 at bound. D not identifiable | `output/inverse_est_metrics.json` |
| Per-patient treatment effect | MGMT TTP p 0.60. Resistance gain 0 of 6 | `output/mgmt_resistance/`, `output/plastic_resistance/results.json` |
| RL vs Stupp on real patients | 25% win | `output/rl_equal_budget/evaluate.json` |
| Track A leak-free classifier | acc 0.638 (random split 0.787 was leaking) | `output/cgat_leak_free.json` |

**Plain reading:** The model finds the whole recurrence area better than the 15 mm plan (+0.08, replicated on two centres). It does not beat the 15 mm plan on enhancing recurrence (null, twice). Growth forecasting and per-patient treatment effect are negative. Do not call any of this a "digital twin" (no validated forecast, no two-way data flow).

**Main limit on enhancing recurrence:** In 32 of 49 UCSF patients the 15 mm plan already covers at least 0.999 of enhancing tissue. There is little room to gain.

## 4. The H-2 to H-4 study, step by step

| Step | What | Status |
|---|---|---|
| H-2 | Train on TUM + MU, lock-test on LUMIERE + RHUH (n=101) | Done. Primary null, secondary positive |
| H-3 dev | Arms C1, C2, C3, U-Net compared on 137 TUM dev patients | Done. C1 +0.030, C3 +0.027, C2 +0.031 (dev, interim) |
| H-3 | C1 frozen, locked test on UCSF (n=59) | Done (ledger 2bd). Primary null, secondary positive |
| H-4 step 1 | Plan rules on pooled dev (G10 rule picked) | Done, **post hoc** (ledger 2be). enh +0.0223, all +0.0869 |
| H-4 step 2 | Pooled retrain, 636 samples, 15 epochs, seed 20261009 | **Training done 2026-10-17** (attempt 1 crashed on disk, ledger 2bg). Model in `kaggle_run/out_h4_retrain/`. SHA-256 `b9846e1f…` checked. Training loss 0.7994 to 0.7164. No test run. |
| H-4 step 3 | New locked test set | **Blocked. No eligible dataset found.** |

Retrain caveat: attempt 2 is a second launch of a run the plan calls one-shot. Ledger 2bg discloses it.

## 5. What is left to do

### 5a. Must do next (H-4)
1. Find a new locked test set that meets A13.4: pre-op MRI, recurrence segmentation with a separate enhancing label, at least 40 patients, usable licence, never scored. Checked and rejected: CFB-GBM v2.0 (no pre-op MRI, one merged GTV), Burdenko-GBM-Progression (post-op baseline, restricted licence, no enhancing label). Not yet checked: BraTS 2025, UCSD-PTGBM, Duke longitudinal (Zenodo 20548528, auto labels), others in `paper/research_prompt_v2.md`.
2. If no set is found: apply A13.7. H-4 stays "post hoc dev" and the headline must say so.
3. If a set is found: write the amendment naming it before any label is read (A13.4), then create `configs/h4_frozen.yaml` and tag `h4-frozen` (A13.6), then run the one allowed test.
4. Write the ledger entry for the finished retrain (training only, no result).
5. Commit `kaggle_run/out_h4_retrain/h4_retrain_train.json` and the log. Decide about the 5.6 MB `.pt` file.

### 5b. Research (prompts ready)
- Run `paper/research_prompt_v2.md` in a deep-research session. Topics: why enhancing recurrence is hard, post-op MRI, extra inputs (diffusion, perfusion, PET), small-data tricks, better loss, new test sets, growth forecast, treatment effect.
- `paper/literature_search_digital_twin.md` is a first pass from search snippets only. Verify each link before citing. Unverified: SAILOR access, any growth model beating persistence on real GBM, NASEM validation wording.

### 5c. Paper and review
- The paper is not written. Figures exist (14 in `paper/figures`, made by `make_paper_figures.py`).
- References are unchecked.
- Outside review (ChatGPT): check each claim against files, record in `paper/chatgpt_review_response.md`. Label interim numbers.
- `REPO_MIND_MAP.md` row 34 is stale (says UCSF test not run). Update it with section 4 above.

### 5d. Repo hygiene
- `Requirements.txt` has merge-conflict markers. Use `pyproject.toml`.
- Script numbers reused (15, 16, 22, 50-56).
- About 3.8 GB untracked run kits and caches (`kaggle_run/*`, `output/pde_manifest/cache_*`). Decide what to ignore.
- Pytest not run recently. CLAUDE.md says 56 tests. Grep finds 101 `def test_`. Fix the count.
- Run `graphify update .` after code changes.

### 5e. Long-term options (all unproven)
- Part 2 (growth forecast): needs 3+ time points per patient. LUMIERE and MU have few. Hierarchical or Bayesian calibration is untested.
- Part 3 (treatment effect): needs per-patient dose and TMZ dates. Candidates: CFB-GBM v2.0 dose maps, Burdenko (restricted). Test-first, simulation-second.

## 6. Rules to keep
- Every number cites a file in `output/` or the ledger.
- Name the weakest point next to each result.
- Words to use: "simulation", "assumed", "post hoc", "exploratory", "interim".
- Never hand-edit output JSON. Fix code and re-run.
- Commit format: `Script N: brief description`. Check `git status --short` before every commit. Push after every commit.

## 7. Key files
| File | Use |
|---|---|
| `analysis_plan_v1.md` | All amendments (A1 to A13). Pre-declared rules |
| `vault/PAPER_FINDINGS_LEDGER.md` | Every result with source file |
| `failed_hypotheses.md` | Negatives |
| `REPO_MIND_MAP.md` | Full inventory with status tags |
| `configs/h*_frozen.yaml` | Freeze files with hashes |
| `paper/research_prompt_v2.md` | Next research prompt |
| `paper/literature_search_digital_twin.md` | First-pass literature notes |
