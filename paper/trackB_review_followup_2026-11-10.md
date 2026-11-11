# Track B review follow-up (2026-11-10)

No patient data, images, masks, or private Duke records were opened. The completed Duke 10-mm analysis stays exploratory and previously used; it is not independent validation.

## Code fixes (commit a765b00; synthetic tests only)

| Review item | Fix | Test |
|---|---|---|
| B-1 gates caller-supplied, unsealed | `build_seal(..., gates=)` records the gate decisions; `run_sealed_protocol` refuses real cases unless the entry gates equal the sealed gates and all are True | `tests/test_fk_track_b_run_guard.py::test_seal_records_gates_and_real_run_requires_same_gates`, `test_seal_without_passed_gates_refused_for_real_cases` |
| B-2 `qc_passed` not tied to QC | `LumiereCase.qc_record` (output of `registration_qc`); real cases are refused without it or if `qc_passed` disagrees; the record's hash is bound into the eligibility manifest that the seal checks | `test_real_case_without_qc_record_refused`, `test_qc_flag_disagreeing_with_record_refused`, `test_qc_record_is_bound_into_manifest` |
| B-3 spacing unchecked | `construct_history_states` requires `voxel_spacing_mm` and refuses anything but 2 mm on all axes | `tests/test_track_b_transition_estimand.py::test_non_2mm_spacing_refused` |
| B-4 empty previous FLAIR | now raises `empty previous FLAIR mask` | `test_empty_previous_flair_and_empty_core` |
| B-5 future-blind test mismatch | new signature and source check for the estimand; README now states that `test_track_b_future_blind.py` covers script 165 only | `test_history_constructor_is_future_blind` |

Scratch run: 105 passed (8 Track B test files), scratch files byte-identical to the worktree. Existing seals created before this change lack `gates` and are refused for real cases, which is intended (fail closed).

## Verified documentation findings

- **B-6 (verified).** `output/track_b_two_history_screen_20261102/screen.json` has three regimes, 27 synthetic cases each, and 21 arms. Gap B8 quoted only growth and shrinkage for the default control. Full picture (Dice with true volume, paired difference vs last-mask, SE):
  - Growth: default control −0.0071 (SE 0.0065). Six `tau28` grid arms were higher than last-mask (best +0.0109, SE 0.0060); `tau7_cap1` +0.0061. These are post hoc choices among 18 grid arms on the same 27 cases, with no multiplicity control.
  - Plateau: all control and grid arms 0.9970-0.9972; copy-last was highest (0.9986, +0.0014).
  - Shrinkage: no arm beat last-mask; default −0.0358 (SE 0.0244), best grid arm −0.0122.
  - Reading: on this synthetic screen the two-history control is not a demonstrated improvement over last-mask extrapolation. Any grid arm would need selection on separate cases before a claim.
- **B-7 (verified).** `src/fk_longitudinal_prep.py` line 3 labels the registration QC rule "Draft (2026-10-31)". Plan text that calls the QC frozen describes a future requirement, not the current state.

## Not changed

The 10-mm Duke result, its private records, and the 2-mm prototype's separate estimand. No Duke-derived output is shared; written terms are still not on file (U1-U3).
