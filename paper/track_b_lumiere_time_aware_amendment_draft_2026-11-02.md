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
- Time-aware fixed growth (20 rates) has **not** been timed.
- PDE inner cross-fitting is **excluded**; the real PDE cost can be higher.
- The central-gap training set is counted at full size.

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
- FiLM U-Net and time-aware fixed growth are not wired into the runner.
- The U-Net reuse fingerprint must include the FiLM recipe and the interval inputs before reuse applies to FiLM models. Not done.
- Empty last core in rate selection: the frozen spec leaves such training patients out and counts them; the prototype raises an error instead. Must be aligned before wiring.
- Nominal-day time step (7 × Δw) for PDE and count arms: proposed, not wired.
- Time-aware fixed growth timing: not done.
- Protocol freeze with SHA-256 of code, cohort, thresholds and records: not done.
- **Patient scoring: BLOCKED.**

No voxels were read, no one was contacted, nothing was downloaded, and no patient was scored for this draft.
