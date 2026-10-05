# Implementation plan for GBM_ISEF_PHD_Level_Masterplan.md

Source: `C:\Users\vihan\Downloads\GBM_ISEF_PHD_Level_Masterplan.md` (sections §1–§210).
Scope: all 37 priority items (Tier S, A, B, §200) in build order (§201).
Rules for every step: patient is the unit (§76), numbers cite `output/` files, negative results are kept, no tuning on LUMIERE after freeze (§79).
New scripts start at `src/94_*`. Reuse existing scripts where listed. "Reuse" names come from `src/` file names; check each one before building on it.

## Gates (a phase does not start until the gate before it passes)
- G1 end of Phase 1: leakage tests green, manifest exists.
- G2 end of Phase 2: solver tests green, baseline table exists.
- G3 end of Phase 3: `v1.0-development-frozen` git tag set BEFORE the external run.
- G4 end of Phase 4: calibration plot exists, selector beats fixed choice or the negative result is logged.

## Phase 0 — Lock the plan (§107, §210.1) — 1–2 days
| Step | Output | Done when |
|---|---|---|
| 0.1 | `analysis_plan_v1.md`: Q1–Q4, primary endpoint (next-scan spatial forecast vs persistence), secondary and exploratory lists (§72), H1–H12 table (§119), subgroup defs, correction rule (Holm / BH, §75) | committed to git with date |
| 0.2 | `vault/claim_ledger.md` (§95), `vault/learned_vs_assumed.md` (§96, §40), `vault/failed_hypotheses.md` (§195) | seeded from `PAPER_FINDINGS_LEDGER.md` |
| 0.3 | ISEF record: current-period work log, Form 7 continuation note, AI-use table (§167–169) | `docs/isef_compliance.md` |

## Phase 1 — Data and leakage (§173, §103, §29–30) — Month 1
| Step | Output | Reuse |
|---|---|---|
| 1.1 | `src/94_patient_manifest.py` → `data/manifests/patient_manifest.csv` (columns from §173) plus dataset manifests (§172) | `70_load_mu_glioma`, `81_label_correct_forecast` cache |
| 1.2 | Missingness + CONSORT flow table + cohort-shift table (§138–140) | new |
| 1.3 | `tests/test_no_leakage.py`, `tests/test_split.py`: disjoint patient IDs, duplicate scans, temporal leakage, fit-on-test (§103) | `86_patient_level_cv` |
| 1.4 | Track A: patient-held-out graph built from train patients only; patient-ID / batch / cell-count baselines; label permutation (§29–30, §66, §128) | `88_cgat_leak_free` |
| 1.5 | Experiment manifest writer (seed, versions, commit, data hash) used by every script (§104, §106) | new `src/common/manifest.py` |
**G1.**

## Phase 2 — Verify solver, build baseline ladder (§7–9, §105, §129–130) — Month 1–2
| Step | Output |
|---|---|
| 2.1 | `tests/test_solver.py` + `src/95_solver_verification.py`: manufactured solution, grid and dt convergence, positivity, conservation, boundary test, tolerance sensitivity |
| 2.2 | Baseline ladder M0–M1 plus extrapolation baselines: persistence, geometric (4 variants), volume-matched, last-rate, linear, exponential, Gompertz (§130) |
| 2.3 | Model ladder M2–M3: isotropic, anisotropic (reuse `74_fractal_aniso_vs_iso`, `run_improved_aniso.py`) |
| 2.4 | DTI ablation set A–E (§87) and the four "why DTI fails" tests A–D (§86); keep E only if external gain |
| 2.5 | Forecastability ceiling + SNR / impossibility test (§181–182) |
| 2.6 | Volume vs shape scored separately (§132); metric suite: Dice, Jaccard, HD95, ASD, volume, centroid, SPE (§22) |
| 2.7 | Statistics module: patient-level bootstrap, paired permutation, Wilcoxon, effect size, per-patient ΔDice, catastrophic-error rate, power sim (§23–24, §73, §76–77) |
**G2.**

## Phase 3 — Rolling-origin, external, failure atlas (§16–21, §18, §27) — Month 2–3
| Step | Output | Reuse |
|---|---|---|
| 3.1 | `src/96_rolling_origin.py`: filter-only forecast, patient-clustered CIs, fixed horizon bins set in advance (§18–19) | `79`, `81` |
| 3.2 | Subgroups from pre-forecast info only: growth / stable / shrink, early post-RT vs late (§20, §186–187) | `79_prespecified_stratified_forecast` |
| 3.3 | Segmentation perturbation (erode/dilate ±1,2, jitter) and registration sensitivity (rigid/affine/deformable) (§13–14) | new |
| 3.4 | Compartment forecasting: enhancing / necrotic / FLAIR / cavity; cavity as explicit state (§82–84) | `70_load_mu_glioma` |
| 3.5 | Negative controls 1–6 (§128) | new |
| 3.6 | Residual maps + residual-vs-covariate analysis (§98, §133) | new |
| 3.7 | Residual-driven model fixes only (§99). Treatment-aware M4 if timing explains error | new |
| 3.8 | **Freeze:** `configs/development.yaml`, checksum, git tag `v1.0-development-frozen` (§101, §192) | — |
| 3.9 | `evaluate_external.py --dataset LUMIERE` loads frozen config only; optional UCSF as "distribution-shift" cohort (§79–80) | `91_lumiere_volume_forecast` |
| 3.10 | Protocol / scanner / interval shift vs performance (§15, §140); RANO-aware labels (§21); "why persistence wins" explanation (§89–90) | new |
| 3.11 | Failure atlas, 13 modes (§27) | new |
**G3.**

## Phase 4 — Uncertainty, sequential update, selector (§5, §10–12, §25–26, §91–93, §146–150) — Month 4
| Step | Output |
|---|---|
| 4.1 | Synthetic recovery of ρ, D vs interval, scan count, noise, size; profile likelihood; parameter correlation; Fisher information (§11). Report identifiable combos if ρ, D not separable |
| 4.2 | Bootstrap posterior → sequential update: ML refit → EnKF → (stretch) particle filter (§5). Reuse `51_inverse_parameter_estimation` |
| 4.3 | Observation model H(u,v,n,e) with segmentation noise (§6) |
| 4.4 | Uncertainty budget: measurement / parameter / model-form / shift (§12) |
| 4.5 | Calibration: coverage at 50/80/90/95, CRPS, NLL; patient-level conformal for scalars only (§22, §25, §93) |
| 4.6 | Selector: start 2-state growth vs non-growth (§148); learn on train, freeze, test external (§149). Compare Always-PDE / Always-persistence / Adaptive (§92). "Growth information content" analysis (§90) |
| 4.7 | Model-form variants: logistic, Gompertz, two-state, heterogeneous D; AIC/BIC + predictive log score (§52, §127) |
| 4.8 | Ablation table: anisotropy, DTI, treatment, molecular prior, uncertainty, selector, sequential update (§100) |
**G4.**

## Phase 5 — Treatment control (§39–58, §121–126) — Month 5
| Step | Output | Reuse |
|---|---|---|
| 5.1 | Latent per-patient kill rate k ~ p(k), resistant fraction f_r ~ prior; observed/derived/inferred/assumed table (§40, §122, §125) | `93_real_rho_kill_sensitivity`, `82/83` |
| 5.2 | Heterogeneity 2×2: growth × treatment response (§123) | extend `93` |
| 5.3 | Belief-state (POMDP) wrapper; sequential posterior over k (§41–42) | `75_probe_and_paced_policies` |
| 5.4 | Policies: standard, MTD (`89_dual_mtd_baseline`), paced, MPC, robust MPC (`52_robust_mpc_controller`), RL (`58_rl_adaptive_steering`, `src/rl/`), oracle RL labelled upper bound (§43, §57) | — |
| 5.5 | Multi-objective utility + Pareto frontier (§45) | new |
| 5.6 | Replace 0.02/day rule: continuous growth rate or train-locked threshold with bootstrap CI (§46–47, §141) | `85_rho_threshold_check` |
| 5.7 | Stress tests: noise 0–20%, monitoring 1–4 wk, parameter misspecification ±10/20/30%, model mismatch (§48–50, §54) | new |
| 5.8 | Simulator holdout: train generator A, test generator B; worst-case regret (§51–53) | new |
| 5.9 | Value of information + stopping rule for probing (§55–56) | `75` |
| 5.10 | Cluster bootstrap / patient-level aggregation (§74) | stats module |
| 5.11 | Optional two-state → S/I/R phenotype model only if identifiable (§124–126) | new |
Treatment is labelled "simulation", never causal (§184).

## Phase 6 — Molecular / spatial integration (§28–38, §62–70, §32) — Month 6
| Step | Output |
|---|---|
| 6.1 | Gene-set programs (pre-declared sources/versions), pathway scores, pseudo-bulk patient-level tests (§31, §33, §68–69) |
| 6.2 | Stability under leave-one-patient-out; batch-correction comparison (§67) |
| 6.3 | Spatial validation on `multiomic-gbm`: Moran's I, neighbourhood enrichment; continuous tumor-state score s∈[0,1] vs zones (§63–64) |
| 6.4 | Heterogeneity metrics: entropy, gradient, core-periphery contrast (§65) |
| 6.5 | Survival as downstream check: continuous predictors, bootstrap, C-index, calibration, events-per-variable cap (§37, §70, §142–144) |
| 6.6 | Molecular→PDE hierarchical prior: Stages A–D (§34, §36). Link is a "transferable prior", not matched patients (§35). If no gain: label exploratory |
| 6.7 | Known limit from vault: zone-expression files likely synthetic (L32) → do not use as evidence |

## Phase 7 — Benchmark, paper, ISEF (§108, §136–138, §156–174, §193–194, §210.17–20) — Month 6+
| Step | Output |
|---|---|
| 7.1 | `gbmbench` CLI: `run --model all --dataset MU`, `--model frozen --dataset LUMIERE`, `--policy all --scenario stress`; outputs metrics, CIs, figures, config, hash (§137) |
| 7.2 | Repo layout + `README`, `LICENSE`, `CITATION.cff`, `environment.yml`, `configs/`, clean-machine run test, second person reruns (§102, §174, §193) |
| 7.3 | Stability stress-test of headline result: seeds, folds, preprocessing (§194) |
| 7.4 | Literature protocol §177 + novelty matrix §178 (real papers only); TRIPOD+AI / PROBAST+AI checklist (§108, §177) |
| 7.5 | Paper from results: strongest result first, failure + external sections, validation matrix (§94), claim hierarchy L1–L4, banned words (§109–110, §154) |
| 7.6 | Board (3–4 findings, full-loop figure §157), notebook (§170), interview drills §158–164 |

## Tier B stretch (only after G4, one at a time)
Neural residual (§59), FNO surrogate with its own error budget (§60–61, reuse `src/neural_pde/`), Bayesian model averaging, target-trial emulation if metadata allow (§185), conformal sets beyond scalars.

## Known risks
- Parts of Track B depend on LUMIERE and MU data already on disk (`data/external/`). Re-check each path before Phase 1.
- Several existing numbers are interim or post hoc. Mark them so in the ledger.
- Scope is about 6 months of work. Phases 0–4 give the strongest science. Phases 5–6 can be cut to "exploratory" without hurting Phases 0–4.
- ISEF 2027 is May 8–14, 2027 (§167). Today is 2026-10-04. Work backward: external run done by ~Jan 2027, paper by ~Mar 2027.
