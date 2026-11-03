# Track B LUMIERE amendment A-2026-11-02-T: time-aware required baselines (DRAFT for review, 2026-11-02)

**Status: DRAFT. Not adopted. Not frozen.** It changes nothing in the frozen runner, the frozen recipes or the MU/Duke protocol. It does not authorize the MRI-QC preflight, voxel reads or scoring. LUMIERE results stay exploratory.

Applies to: the LUMIERE scoring draft `paper/track_b_lumiere_scoring_amendment_draft_2026-11-01.md` (also a draft). Prototype code and synthetic tests: commit `ff4b9f0` (`src/fk_time_aware_baselines.py`, `tests/test_fk_time_aware_baselines.py`). Not wired into any runner.

## 1. Change

The LUMIERE required mask baselines become:

| # | Required comparator | Status vs the frozen design |
|---|---|---|
| 1 | Copy-last, clipped to support | unchanged |
| 2 | **Time-aware fixed growth** (§2b) | **replaces** gap-blind fixed growth |
| 3 | PDE-only | unchanged |
| 4 | **FiLM time-aware U-Net** (§2a) | **replaces** the gap-blind 3-channel U-Net |

Gap-blind fixed growth and the gap-blind U-Net are not required LUMIERE comparators under this amendment. They may be reported as secondary, labelled "gap-blind". The count model, volume-trend, geometric train-rate and ablations are unchanged secondary comparators. The MU/Duke 90-day design keeps its own frozen baselines.

**Same for every arm:** the same QC-passing patients, the same fold per patient (scoring draft §3b: sorted by Δw, round-robin into folds 0–4), the same origin-only support S, the same 2 mm grid and the same target (DeepBraTumIA labels 1|2 at visit 3). All fitting uses outer-training patients only. No arm receives future tumor size.

## 2. Recipes (exactly as implemented in `ff4b9f0`)

### 2a. FiLM time-aware U-Net

- **Backbone:** the frozen `UNet3D` and `UNET_RECIPE`: widths 16/32/64/128; 3 image channels (visit-1 core, visit-2 core, origin domain); Adam 1e-3, betas (0.9, 0.999), weight decay 0; 300 epochs; 3 seeds (20261102 + fold + 0/1000/2000); soft Dice + BCE; flips and axial 90° rotations; threshold 0.5; clip to the origin domain.
- **Interval inputs:** z_f = log(Δw(v2→v3) / 13) and z_h = log(Δw(v1→v2) / 13), from published week-label differences. The constant 13 is fixed (`INTERVAL_SCALE_WEEKS`). Nothing is fitted in the scaling.
- **FiLM network:** Linear(2→16) → LeakyReLU(0.01) → Linear(16→1,408). It gives a separate γ and β for each channel of each of the 14 InstanceNorm layers (704 channels). The FiLM network has 23,984 parameters; the full model has 1,426,817 (`output/track_b_timing_pilot_20261031/film_unet_timing_record.json`).
- **Applied** after each InstanceNorm and before LeakyReLU: x ← IN(x) × (1 + γ) + β.
- **Zero init:** the last FiLM Linear starts at zero, so the untrained model equals the frozen U-Net with the same backbone weights.
- **Training:** FiLM and backbone are trained jointly on outer-training patients only.

### 2b. Time-aware fixed growth

- **Interval:** t_i = 7 × Δw_i nominal days, Δw_i = v2→v3 week-label difference. This is an approximation, not an exact interval.
- **Rate grid:** k ∈ {0.00, 0.01, …, 0.19} mm per nominal day (20 values, `RATE_GRID_MM_PER_DAY`).
- **Radius:** r_i = k × t_i mm. At Δw = 15 the largest radius is 19.95 mm. The code refuses any grid that exceeds the frozen 20 mm cap at 15 weeks.
- **Mask:** F_i = { v ∈ S : Euclidean distance(v, last core) ≤ r_i }, as in the frozen fixed-growth spec.
- **Selection:** one k per outer fold, from that fold's outer-training patients only; best mean patient-macro plain Dice; ties within 1e-9 go to the smaller k; no inner folds.
- **Flags:** k = 0 → "equals copy-last"; k = 0.19 → "boundary hit". The grid is not widened after scoring starts.

## 3. Synthetic verification done (code behaviour only; not patient evidence)

Scratch copy, Python 3.14: 45 passed (`tests/test_fk_time_aware_baselines.py` 6, plus the learned-baselines and fixed-growth suites).
- Zero-init FiLM with identical backbone weights: logits, loss and every backbone gradient match the frozen U-Net by `torch.equal` (exact element equality) and by byte comparison, for three interval pairs.
- With fixed non-zero FiLM weights and masks held fixed, changing only the forecast or only the history interval changes the output.
- Rate selection uses outer-training patients only; the 20 mm cap is enforced; ties go to the smaller k.

Not done: the 300-epoch synthetic response studies and their negative controls (not adoption gates by user decision).

## 4. Pass rule (replaces the comparator list in the scoring draft §3a)

> The PDE candidate passes only if, on patient-macro plain Dice, every max-T-adjusted paired 95% patient-bootstrap interval for the candidate versus each of the four required comparators — copy-last, time-aware fixed growth (§2b), PDE-only and the FiLM time-aware U-Net (§2a) — lies strictly above zero (appendix §1.3). If any comparator is unavailable for any patient, the comparison is incomplete, not passed.

**Claim if it passes:** "On LUMIERE (exploratory; at most 28 patients; automated DeepBraTumIA masks; 11–15 week-label window), the PDE candidate outperformed the four specified comparators on patient-macro plain Dice."

**Not claimed:** that the comparison isolates a causal PDE effect; superiority over other time-aware or learned methods; agreement with true tumor; generalisation, novelty or clinical value.

## 5. Run time (wall-clock hours on this PC, one job at a time; not CPU-core-hours)

All four analyses (three prior supports + central-gap refit), Amendment R reuse = 2 FiLM training sets of 100,800 steps each.

| Scenario | FiLM U-Net | PDE (scenario projection) | Total |
|---|---|---|---|
| S1: FiLM step from this session, 4.141 s | 232.0 h | 323.4–466.2 h | **555.4–698.2 h** |
| S2: earlier pilot frozen step 7.973 s × measured FiLM/frozen ratio 1.038 | 463.6 h | 323.4–466.2 h | **787.0–929.8 h** |

Sources: `film_unet_timing_record.json` (FiLM 4.141 s, frozen 3.989 s same session, ratio 1.038), `lumiere_timing_record.json` (frozen 7.973 s), scoring draft §5a and `lumiere_history_gap_stress.json` (PDE scenarios).

**Caveats:**
- **Neither range is a bound.** Both are projections.
- The absolute U-Net step time varied about 2× between sessions on this PC (3.989 s vs 7.973 s for the same frozen model).
- Time-aware fixed growth rate selection (20 rates), **measured (synthetic)** at 91 × 109 × 91: 19.4, 20.9 and 23.7 s for 22, 23 and 22 training cases (`output/track_b_timing_pilot_20261031/timed_fixed_growth_timing_record.json`). **Extrapolated:** about 5 folds × 21 s ≈ 105 s per analysis, about 0.12 h over four analyses. Negligible next to the ranges above; not added to them.
- PDE inner cross-fitting is **excluded**; the real PDE cost can be higher.
- The central-gap training set is counted at full size.

## 5a. Empty visit-2 core rule (draft, 2026-11-02)

- A patient whose visit-2 (origin) core is empty is **excluded from every arm**. This check uses the visit-2 history mask only.
- Order: registration QC exclusions, then empty visit-2 exclusions, **before any target mask is read**. Both counts are reported.
- If fewer than **25** patients remain after all frozen eligibility and QC rules: **stop**. No target is read, nothing is scored, no re-folding, no relaxed rule.
- Remaining patients keep their fold (scoring draft §3b).
- Code: `src/fk_lumiere_runner.py:eligibility` (synthetic tests only). This replaces the earlier gap note about the prototype raising on an empty last core: the runner excludes such patients before rate selection.

## 5b. Integration status (synthetic only, 2026-11-02)

- `src/fk_lumiere_runner.py:run_lumiere_cv`: guarded LUMIERE path with the four proposed comparators plus the PDE candidate; common eligible patients, folds, support and targets; time-aware fixed growth rate selected per outer fold on training patients only; FiLM U-Net trained per outer fold on training patients only; nominal days 7 × Δw passed to the PDE callables (`nominal_days`); lazy target reads after the stop check.
- **Fail closed:** real cases are refused unless all 11 `LUMIERE_GATES` are recorded True; synthetic and real cases cannot be mixed.
- **FiLM reuse fingerprint:** `film_training_fingerprint` = the frozen U-Net fingerprint (inputs, targets, fold, recipe, environment) + FiLM config + each patient's (forecast, history) weeks. Reuse only on a full match.
- **Not defined here:** the PDE-only and PDE-candidate masks. They are injected callables; the LUMIERE PDE-only definition and the PDE wiring remain open.
- MU/Duke runner unchanged.
- Tests (scratch copy, Python 3.14): `tests/test_fk_lumiere_runner.py` 5 new; with `test_fk_time_aware_baselines.py`, `test_fk_learned_baselines.py` and `test_fk_mu_pipeline_components.py`: all passed (44 + 19 in two runs; 3 first-run failures were a missing `paper/` folder in the scratch copy, fixed and re-run).
- Not added: bootstrap intervals and the max-T pass test in this runner.

## 5c. PDE-only comparator: proposed definition and cutoff check (2026-11-03; STOPPED)

**Proposed definition.** Same history-fitted joint posterior, parameter draws, forecast interval (7 × Δw nominal days), solver, domain and seeds as the candidate, from the **same solver run**. For each raw concentration draw u: mask = {u ≥ c} ∩ domain, **before** volume matching and boundary perturbation. Combine draws with the candidate's inclusion rule (mean inclusion ≥ 0.5, ties included). Clip to the same support. No test target or test volume is used. It shares the candidate's posterior calibration but **omits its forecast-volume model and boundary corrections**. The matched-volume no-PDE ranker stays a separately named ablation. Code: `src/fk_hybrid_forecast.py`, option `pde_only_concentration_cutoff` (default `None` leaves the output unchanged; tests `tests/test_fk_pde_only_cutoff.py`, 5 passed; existing hybrid and MU pipeline tests 31 passed, scratch copy). Not wired into the LUMIERE runner.

**Cutoff check with c = 0.5: FAIL** (record `output/track_b_pde_only_cutoff_check_20261103.json`; synthetic only).
- Predeclared rule: PASS if the concentration scale is [0, 1] (carrying capacity 1, initial core = 1) **and** in every synthetic posterior case ≥ 90% of draws give a non-empty, non-full 0.5 mask.
- Scale: concentrations 0.0–1.0 over 16 (ρ, D) pairs × 3 horizons. The scale part passes. But 7 of 48 fields had **no voxel ≥ 0.5**: all at D = 1.0 mm²/day with ρ = 0.001–0.03/day (diffusion-dominated; peak u 0.08–0.51).
- Posterior draws (32³ grid, 40 draws, reduced posterior settings): empty 0.5 masks in **11/40** (slow case), **7/40** (moderate) and **23/40 plus 3 full-domain** (fast) draws. In the fast case the combined PDE-only mask was **empty**.
- Reading: the joint posterior includes high-D, low-ρ draws on the ρ–D ridge, where the field spreads below 0.5. A fixed 0.5 cutoff then calls the tumor gone. 0.5 is not a meaningful cutoff for these draws.
- Weak points: synthetic spheres; small grid; reduced posterior settings.

**Not done at the time of the check:** wiring PDE-only into the runner, and any cutoff change. The 0.5 cutoff was **not** silently changed.

### 5d. Proposed rule after the failed check: training-selected cutoff (2026-11-03; not adopted)

**Arm name and meaning.** PDE-only = a PDE posterior forecast with a training-selected concentration cutoff, **without** forecast-volume matching or boundary perturbation. It shares the candidate's posterior draws and solver run. The matched-volume no-PDE ranker stays a separately named ablation. No relative half-peak cutoff is used.

**Rule.**
- Grid: c ∈ {0.05, 0.10, …, 0.50} (`PDE_ONLY_CUTOFF_GRID`).
- Per outer fold, choose one c from **outer-training patients only**: mean patient-macro plain Dice of their PDE-only masks ((inclusion probability at c ≥ 0.5) ∩ support) against their targets. Ties within 1e-9 go to the **larger** c.
- Apply that c to each held-out patient's raw concentration draws; combine with the candidate's inclusion vote (≥ 0.5); clip to the common support.
- A valid **empty** PDE-only forecast stays empty and is scored by the frozen empty-mask rule. It is not replaced by persistence.
- A PDE `ForecastFailure` uses the frozen persistence fallback for the candidate and PDE-only, and is recorded per patient.
- Recorded per fold: chosen c, grid-boundary hit (0.05 or 0.50), training Dice per c, training patients whose calibration run failed. Per held-out patient: empty draws at c, empty PDE-only mask, fallback flag.

**Code path check and correction.** The earlier §5c sentence ("the inner calibration runs already compute their concentration fields, so masks for all 10 cutoffs could be recorded there") needed a correction. The saved calibration-only output (`pre_boundary_probability`) is **volume-matched**, so it does **not** give raw cutoff masks. But the calibration-only loop does compute the raw concentration field of every unique posterior draw before volume matching (`src/fk_hybrid_forecast.py`, `_batched_parameter_concentrations` loop). The new option `pde_only_cutoff_grid` thresholds that raw field at every grid cutoff inside the same loop, in both full and calibration-only runs. So **no extra solver run** is needed.

**Is the reuse safe?**
- No target and no held-out patient enters it: in `src/fk_lumiere_runner.py`, a training patient's calibration call receives only the other outer-training patients; held-out targets are never passed. A test checks that one outer fold is always absent from every calibration call.
- **Mismatch (state with results):** in the frozen MU design the inner calibration runs use the initial boundary covariance and a fixed 90-day horizon. For LUMIERE, each training patient's calibration run must use its own 7 × Δw nominal horizon (`nominal_days`), which this runner passes. The posterior of a training patient's calibration run may use a different boundary covariance from the held-out forecast (the posterior takes the covariance as an input), so the cutoff is chosen on slightly different posteriors from the ones it is applied to.
- Default output unchanged: with `pde_only_cutoff_grid=None` the output has no new keys; with the grid, candidate `mask_draws` and calibration `pre_boundary_probability` are identical (tests).

**Cost (synthetic; `output/track_b_pde_only_grid_timing_20261103.json`).**
- Measured: 0.038 s per unique raw field for all 10 cutoffs at the LUMIERE 91 × 109 × 91 grid.
- End-to-end calibration-only runs (32³, 40 draws): 23.9 s and 16.6 s without the grid, 20.6 s and 15.9 s with it. The overhead is below run-to-run noise.
- Extrapolated upper figure: ≤ 200 unique draws × 0.038 s ≈ 7.6 s per forecast; 140 forecasts per analysis ≈ 0.30 h; ≈ 1.2 h over four analyses.

**Revised projection (wall-clock hours on this PC, four analyses; projections, not bounds).** S1 555.4–698.2 h, S2 787.0–929.8 h (§5), plus about 0.12 h time-aware fixed growth and up to about 1.2 h PDE-only cutoff bookkeeping: **about 556.7–699.5 h (S1) or 788.3–931.1 h (S2)**. The PDE part still excludes inner cross-fitting; the FiLM step time still varied about 2× between sessions.

**Verification (scratch copy, Python 3.14):** `tests/test_fk_lumiere_runner.py` 7 and `tests/test_fk_pde_only_cutoff.py` 6, with the hybrid forecast, full-hybrid gate, hybrid stress and MU pipeline-component suites: 44 passed. Covered: fold isolation of calibration and held-out calls, ties to the larger cutoff, empty PDE-only stays empty and scores 0 against a non-empty target, ForecastFailure fallback recorded, unchanged default candidate output.

## 5e. Paired bootstrap and max-T intervals: implementation choices (frozen in this draft 2026-11-03, before coding; not adopted)

Appendix §1.3 says: paired 95% cluster-bootstrap intervals, whole patients resampled, max-T adjustment across all required comparator contrasts. It does not give the statistic, the number of resamples or the seed. Chosen now, before any patient result exists:
- **Contrasts:** k = 1..4, d_ik = Dice(candidate, i) − Dice(comparator k, i) for copy-last, time-aware fixed growth, PDE-only and FiLM U-Net, over the same scored patients.
- **Resampling:** B = 10,000 resamples of patient indices with replacement, the **same** indices for all four contrasts (paired). Generator `numpy.random.default_rng(20261110)`.
- **Statistic:** single-step max-T on studentized means. se_k = standard deviation (ddof = 1) of the B bootstrap means of contrast k. T_b = max_k |mean_b,k − mean_k| / se_k. c = 95th percentile of T_b (`numpy.quantile`, default linear method).
- **Adjusted interval:** mean_k ± c × se_k.
- **Pass:** every adjusted lower bound > 0 (strict). **Incomplete:** any required arm missing for any scored patient, fewer than 2 patients, or any se_k = 0 / non-finite. Incomplete is not a pass.
- Weak points: percentile choice and B are conventions; with ≤ 28 patients the bootstrap distribution is coarse.

Code: `src/fk_lumiere_stats.py` (synthetic tests only).

## 5f. Runner completion on synthetic data (2026-11-03; not adopted)

- **PDE provider** `src/fk_lumiere_pde_provider.py:LumierePDEProvider`. It follows the frozen MU outer-nuisance rules (the MU runner is unchanged): volume nuisance from outer-training patients only; inner patient-held-out boundary calibration on outer-training patients only; an inner fallback or fewer than two valid residuals makes the fold not ready, so every held-out forecast uses persistence; expected model failures and empty history masks raise `ForecastFailure`. LUMIERE days are nominal (0, 7 × history weeks, + 7 × forecast weeks), so the inner calibration forecasts each training patient at its own nominal horizon (`forecast_horizon_days=None`). Every call carries the PDE-only cutoff grid, so the inner calibration runs give the training rows for the cutoff selection.
- **One shared solver call** per held-out patient gives the candidate (`mask_draws`), PDE-only (cutoff grid) and the separately named **matched-volume no-PDE ablation** (`no_pde_ranker_mask_draws`). The ablation is scored but is not a required comparator.
- **Count arm** (frozen Poisson GLM, `fit_count_model`): fitted per outer fold on outer-training patients only; counts = `new_label_pieces(visit 2, visit 3)` inside support; offset log(7 × Δw); test μ = rate × 7 × Δw. Scored by Poisson deviance, absolute count error and Brier score of P(any new piece). It is never a Dice mask.
- **Pass test:** §5e max-T on the four required contrasts; missing arms give "incomplete".
- **Tests** (scratch copy, Python 3.14): `tests/test_fk_lumiere_runner.py` 9 and `tests/test_fk_lumiere_stats.py` 4, 13 passed. Covered: fold isolation (the held-out fold is absent from every fit view; each held-out patient is forecast once); identical patient sets for all arms; shared PDE call; empty PDE-only stays empty; fallback recorded with reason for candidate, PDE-only and ablation; count arm separate from Dice; pass, not passed and incomplete logic; seeded max-T. One test runs the real hybrid forecast on a 16³ grid: all 5 folds were model-ready, no fallbacks, chosen cutoff 0.10 in each fold; the pass test was "incomplete" because the 1-epoch tiny U-Net was unavailable.

**Remaining ambiguities (recorded, not chosen):**
1. `calibration_minimum_patients`: the MU runner requires at least 5 and uses 8 by default; the LUMIERE value is not set. With about 20–23 QC-passing training patients per fold, 8 is the MU default.
2. Boundary-calibration settings for LUMIERE (inner folds 5, ℓ candidates 5/10/20 mm, initial covariance, seeds): copied from the MU runner; not separately frozen for LUMIERE.
3. Count-arm voxel volume uses the 2 mm grid spacing; the "new piece" margin and minimum size are the frozen MU values (2 mm, 2 voxels).
4. The count arm has no pass rule; it is reported only.
5. Prior-sensitivity runs and the central-gap refit are not wired in the LUMIERE runner.

## 6. Gates (all open or blocked; nothing passes by this draft)

**Source gates**
- LUMIERE derived-output license terms: not stated. OPEN.
- DeepBraTumIA version: not stated (ARTORG page checked; NITRC page unreachable). OPEN.
- DeepBraTumIA label 0 meaning: not stated. OPEN. No cavity label (stated limit).

**Approval gates**
- School written approval (date, approver, conditions): not saved. OPEN.
- User approval of the exact 84-file MRI-QC preflight read scope: not given. **Preflight BLOCKED.**
- User adoption of this amendment and of the scoring draft: not given. OPEN.

**Cohort gates**
- Header inventory: 28 eligible, 84 visits; voxel content unchecked. Headers only.
- Real-patient registration QC: BLOCKED until preflight approval.
- Real-patient support check (non-empty, on grid): OPEN.
- Minimum of 25 QC-passing patients with stop rule and no re-folding: proposed, not adopted.
- Frozen cohort list and hash: not done.
- History-gap outliers (Patient-086, -070, -078): keep-and-flag with predeclared sensitivity proposed, not adopted.

**Protocol gates**
- Pass-rule reading (no numeric margin; adjusted intervals strictly above zero): needs user confirmation.
- FiLM U-Net and time-aware fixed growth: wired into the guarded LUMIERE runner (synthetic only, §5b). Not adopted.
- FiLM reuse fingerprint: done (§5b).
- Empty visit-2 core: rule written and coded (§5a); not adopted.
- Nominal-day step (7 × Δw): passed to the PDE callables; the count arm uses the nominal offset (§5f).
- LUMIERE PDE-only mask: 0.5 cutoff FAILED (§5c); training-selected cutoff rule proposed and wired in the guarded runner with synthetic tests (§5d); not adopted. Real PDE provider: wired (§5f).
- Bootstrap intervals and max-T pass test: implemented with choices frozen in §5e (synthetic only); not adopted.
- Real PDE provider and count arm: wired (§5f, synthetic only); remaining ambiguities listed in §5f.
- Time-aware fixed growth timing: done (§5).
- Protocol freeze with SHA-256 of code, cohort, thresholds and records: not done.
- **Patient scoring: BLOCKED.**

No voxels were read, no one was contacted, nothing was downloaded, and no patient was scored for this draft.
