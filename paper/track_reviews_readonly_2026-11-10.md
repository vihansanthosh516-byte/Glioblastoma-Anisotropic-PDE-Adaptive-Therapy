# Read-only track reviews (2026-11-10)

Three independent read-only reviewer agents (Track A, B, C). They edited no files, ran no tests or analyses, and opened no patient data. The lead checked the findings listed as "verified" against the files. Findings not marked verified are reviewer claims that were not re-checked. No frozen result changes.

## Track A (A5, Yu 2020)

Overall: code matches the frozen protocol on the points checked; no donor leakage found; START_HERE numbers match `output/trackA_yu2020/result.json`.

| # | Sev | Finding | Evidence | Lead check |
|---|---|---|---|---|
| A-1 | Med | Within-cell-type check (PROTOCOL §8-9) is in result.json but is never mentioned in START_HERE or gap_status; only Myeloid scored | result.json ~7198 | not re-checked |
| A-2 | Med | Mean AURC gap driven by two donors (GS1 about −0.10, GS5 about −0.05); GS3, GS13 favor logistic | per-donor rows in result.json | not re-checked |
| A-3 | Low | Manifest `gbmap_version` is a hard-coded string, not derived from the snapshot | `src/189_trackA_yu2020_lodo.py:296` | **verified** |
| A-4 | Low | AURC integrates coverage from ceil(N/2)/N to 1, not normalized; values are not on the usual 0-1 scale. Matches protocol; docs need a scale note | `189_...:486-493` | **verified** |
| A-5 | Info | 5-of-6 rule implemented correctly; 4 of 6 donors | result.json 7144-7148 | **verified** |
| A-6 | Low | `ablation_no_gate` sets lambda_1 = 0 but still trains the gate; name misleads | `189_...:643-644` | not re-checked |

## Track B

Overall: `screen.json` matches the docs for the growth and shrinkage numbers quoted. Completed Duke result not affected by any finding (private records not opened).

| # | Sev | Finding | Evidence | Lead check |
|---|---|---|---|---|
| B-1 | High (future runs) | Scoring gates come from a caller-supplied dict and are not sealed | `src/fk_track_b_run_guard.py:227` | **verified** |
| B-2 | High (future runs) | `qc_passed` is a free case flag, not tied to `registration_qc` | `fk_lumiere_runner.py:62, 98`; `fk_longitudinal_prep.py:164` | not re-checked |
| B-3 | High (prototype) | Voxel spacing hard-coded at 2 mm with no input check | `src/track_b_transition_estimand.py` `sampling=2.0` | **verified** (line 32) |
| B-4 | Med | Empty previous FLAIR silently makes all risk voxels "new"; should raise | `track_b_transition_estimand.py:31-33` | **verified** |
| B-5 | Med | `test_track_b_future_blind.py` tests `165_signed_transition_cv.py`, not the transition estimand; README lists it as an estimand check | `tests/test_track_b_future_blind.py:13` | **verified** |
| B-6 | Med | gap B8 omits the plateau regimes and the other grid arms; growth and shrinkage differences are not significant (z about −1.1, −1.5) | `screen.json` | not re-checked |
| B-7 | Low | Plan says registration QC is frozen; code labels it "Draft (2026-10-31)" | `fk_longitudinal_prep.py:1` | not re-checked |

## Track C

| # | Sev | Finding | Evidence | Lead check |
|---|---|---|---|---|
| C-a | — | Reviewer said the f0 prior is too narrow in the integrated path | `216_...:147-148` default | **Rejected for the integrated benchmark**: `211_...:993` passes `INTEGRATED_BELIEF_RANGES` (f0 0.000497-0.9) to every MPC arm. True only for direct default construction of `SignalAmbiguityBeliefMPC` |
| C-b | Med | B0 generator reaches 0.65; filter prior 0.05-0.5; sim docstring assumes B0 ≤ 0.5 | `trackc_sim.py:15-16, 43` | **verified**; covered by the prior amendment draft |
| C-c | — | Config q grid 0.5-0.9 vs world levels 0.6-0.9 | config lines 29-36, 95-100 | **Not a drift**: 0.5-0.9 is the controller model set (`q_signal_model_set`); world draws use `primary_q_signal_sampling_levels` 0.6-0.9. Already covered by C9 |
| C-d | Med | Liu-West resample-move has no Metropolis-Hastings correction; it is an approximation not labelled as one | `trackc_sim.py:256-258` | **verified** (code has no accept step) |
| C-e | Low | Linear parameters use a logit clipped at ±20, not a reflection | `trackc_sim.py:197-199` | **verified**; recorded in amendment rev 2 |
| C-f | Med | Signal-model residual is N(0,1) only under the fixed f0 envelope, not under each world's f0 law | `trackc_integrated_signal_model.py:18, 130-147` | not re-checked |
| C-g | Med | GE-betting bisection assumes monotone accept sets; not proven in code | `trackc_ge_betting.py:140-171` | not re-checked; belongs in the independent statistical review (C3) |
| C-h | Info | GE support [−790, 790] is correct from utility [−60, 730]; it is loose against a world SD near 30 days, which is the root of C4 | config 380-383 | **verified** |
| C-i | Low | gap C3 wording "A valid candidate now exists" overstates a synthetic screen without review | gap_status C3 | wording issue; propose "passed a synthetic screen" |
| C-j | Info | Signal likelihood conditional variance and molecular flip are correct; noise CRN holds inside `simulate`; cross-arm seed reuse in 211 not verified | `trackc_integrated_signal_model.py:224, 303-308` | not re-checked |

## Not covered

Full test suite; hash checks; private Duke records; GE-betting paper equations; 208 and world runner/summary in full; cross-arm CRN in 211.
