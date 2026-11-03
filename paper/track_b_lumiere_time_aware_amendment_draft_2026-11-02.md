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

**Not done (stopped as instructed):** wiring PDE-only into the runner, and any cutoff change. The cutoff is **not** silently changed.

**Proposal for review (not adopted):** select c per outer fold from outer-training patients only.
- Grid: c ∈ {0.05, 0.10, …, 0.50} (10 values).
- Score: mean patient-macro plain Dice of the PDE-only mask against training targets.
- Ties within 1e-9 go to the **larger** c (closest to the prespecified 0.5).
- Cost: needs PDE-only masks for training patients. The inner calibration runs already compute their concentration fields, so masks for all 10 cutoffs could be recorded there without a new solver run; this is not implemented or timed.
- Limit: draws whose peak concentration is below 0.05 stay empty at every c.
- Alternative (not training-selected): a per-draw relative cutoff, u ≥ 0.5 × max(u). Needs its own decision.

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
- Nominal-day step (7 × Δw): passed to the PDE callables; the count arm is not wired in the LUMIERE runner.
- LUMIERE PDE-only mask: definition proposed (§5c); 0.5 cutoff FAILED the synthetic check; cutoff rule OPEN; real PDE wiring OPEN.
- Bootstrap intervals and max-T pass test in the LUMIERE runner: OPEN.
- Time-aware fixed growth timing: done (§5).
- Protocol freeze with SHA-256 of code, cohort, thresholds and records: not done.
- **Patient scoring: BLOCKED.**

No voxels were read, no one was contacted, nothing was downloaded, and no patient was scored for this draft.
