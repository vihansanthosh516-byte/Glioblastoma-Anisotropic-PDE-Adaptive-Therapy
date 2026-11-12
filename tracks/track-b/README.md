# Track B — longitudinal MRI

## Current role

Track B contains longitudinal MRI forecasting and transition analyses. A private Duke patient-mask analysis was completed on 2026-11-04; its cohort and outcomes are exposed exploratory evidence, not an untouched test set. The integrated phase will use only aggregate summaries that the exact access agreement permits. It will not reopen patient images, masks, ADC maps, or patient-level records.

LUMIERE forecast scoring ended in a final no-go on 2026-11-11 under the frozen protocol: 7 of 28 patients passed registration QC, below the minimum of 25. No coverage count, tumor label, target mask, or score was used ([record](../../paper/trackB_lumiere_registration_nogo_2026-11-11.md)). The coverage-aware runner code below is kept but has no active cohort.

## Primary records

- [ACRIN 6677 area forecast result and sRCBV linkage no-go](../../paper/trackB_acrin6677_advanced_mri_linkage_result_2026-11-11.md): the area forecast failed its frozen rule; the sRCBV audit found 15 eligible patients against a minimum of 25 and stopped before measuring-value analysis.
- [ACRIN sRCBV protocol](../../paper/trackB_acrin6677_advanced_mri_protocol_2026-11-11.md) and [aggregate audit output](../../output/track_b_acrin6677_advanced_mri_linkage_20261111/linkage.json).
- [Duke readiness and source record](../../paper/track_b_duke_readiness_2026-11-03.md)
- [Duke source ledger](../../paper/track_b_duke_source_ledger_2026-11-02.md)
- [FLAIR transition prior-art update](../../paper/track_b_flair_transition_prior_art_update_2026-11-07.md)
- [Duke go/no-go packet](../../paper/track_b_duke_go_no_go_packet_2026-11-02.md)
- [Integrated plan (active, revised 2026-11-08)](../../paper/trackABC_integrated_research_plan_revised_2026-11-08.md)

## Main implementation and checks

- `src/track_b_transition_estimand.py`, `src/track_b_run_provenance.py`, `src/fk_duke_archive_inventory.py` — transition and provenance prototypes.
- `src/fk_longitudinal_prep.py`, `src/fk_track_b_run_guard.py` — registration preparation and guarded-run components.
- `tests/test_track_b_transition_estimand.py`, `tests/test_track_b_run_provenance.py`, `tests/test_track_b_manifest_triplets.py`, `tests/test_track_b_geometry_contract.py`, `tests/test_track_b_future_blind.py` (covers `src/165_signed_transition_cv.py` only), `tests/test_fk_track_b_run_guard.py` — focused checks. The transition estimand's own future-blind check is in `tests/test_track_b_transition_estimand.py`. Since 2026-11-10 the LUMIERE forecast runner scores only on geometry-derived common scan coverage (`src/track_b_coverage_scoring.py`; synthetic-tested), and the guard seals gate decisions, binds each real case's saved registration QC record, and the estimand refuses non-2-mm spacing and an empty previous FLAIR mask.
- `output/track_b_feasibility_pilot_20261102/` and `output/track_b_two_history_screen_20261102/` — synthetic development records only.
- Private Duke files remain in access-controlled local storage and are not listed, copied, hashed, or committed here.

## Interpretation and access boundary

The private Duke result has already been examined. A saved-record audit does not validate image registration, target coverage, labels, or clinical meaning. The newer transition-estimand prototype is not a replication of the completed 10-mm analysis. The 2026-11-09 review confirmed from the private run records that target masks were read and the primary endpoint scored on 2026-11-04 (33 patients, 101 windows, 10-mm rule), so any further analysis of this cohort is exploratory. One user-authorized case was timed for runtime only (one rigid registration about 13 s; no target file read, no score computed). ACRIN 6677 later supplied a fresh public clinical-table path: the area forecast failed, and the sRCBV exact-link audit found 15 of the required 25 patients and stopped before measurement magnitudes. No active patient-level Track B forecast path remains. The proposed two-history control is not a proven strong comparator: in the synthetic screen its default setting scored below the simpler last-mask control in the growth and shrinkage regimes (`output/track_b_two_history_screen_20261102/screen.json`). The accepted Duke licence option, the school approval, and the sharing terms are not on file; the user must supply them (see the [open-gap record](../../paper/gap_status_2026-11-09.md)). Check the agreement before using or sharing any derived summary; if permission is unclear, omit Duke-derived summaries and use the plan's standard-normal synthetic MRI signal fallback.

- [ACRIN enhancing-area forecast](../../paper/trackB_acrin6677_result_2026-11-11.md) failed the frozen rule in the amended exploratory analysis.
