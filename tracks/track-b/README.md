# Track B — longitudinal MRI

## Current role

Track B contains longitudinal MRI forecasting and transition analyses. A private Duke patient-mask analysis was completed on 2026-11-04; its cohort and outcomes are exposed exploratory evidence, not an untouched test set. The integrated phase will use only aggregate summaries that the exact access agreement permits. It will not reopen patient images, masks, ADC maps, or patient-level records.

## Primary records

- [Duke readiness and source record](../../paper/track_b_duke_readiness_2026-11-03.md)
- [Duke source ledger](../../paper/track_b_duke_source_ledger_2026-11-02.md)
- [FLAIR transition prior-art update](../../paper/track_b_flair_transition_prior_art_update_2026-11-07.md)
- [Duke go/no-go packet](../../paper/track_b_duke_go_no_go_packet_2026-11-02.md)
- [Integrated plan (active, revised 2026-11-08)](../../paper/trackABC_integrated_research_plan_revised_2026-11-08.md)

## Main implementation and checks

- `src/track_b_transition_estimand.py`, `src/track_b_run_provenance.py`, `src/fk_duke_archive_inventory.py` — transition and provenance prototypes.
- `src/fk_longitudinal_prep.py`, `src/fk_track_b_run_guard.py` — registration preparation and guarded-run components.
- `tests/test_track_b_transition_estimand.py`, `tests/test_track_b_run_provenance.py`, `tests/test_track_b_manifest_triplets.py`, `tests/test_track_b_geometry_contract.py`, `tests/test_track_b_future_blind.py`, `tests/test_fk_track_b_run_guard.py` — focused checks.
- `output/track_b_feasibility_pilot_20261102/` and `output/track_b_two_history_screen_20261102/` — synthetic development records only.
- Private Duke files remain in access-controlled local storage and are not listed, copied, hashed, or committed here.

## Interpretation and access boundary

The private Duke result has already been examined. A saved-record audit does not validate image registration, target coverage, labels, or clinical meaning. The newer transition-estimand prototype is not a replication of the completed 10-mm analysis. The 2026-11-09 review confirmed from the private run records that target masks were read and the primary endpoint scored on 2026-11-04 (33 patients, 101 windows, 10-mm rule), so any further analysis of this cohort is exploratory. One user-authorized case was timed for runtime only (one rigid registration about 13 s; no target file read, no score computed). No fresh labelled longitudinal cohort was found. The proposed two-history control is not a proven strong comparator: in the synthetic screen its default setting scored below the simpler last-mask control in the growth and shrinkage regimes (`output/track_b_two_history_screen_20261102/screen.json`). The accepted Duke licence option, the school approval, and the sharing terms are not on file; the user must supply them (see the [open-gap record](../../paper/gap_status_2026-11-09.md)). Check the agreement before using or sharing any derived summary; if permission is unclear, omit Duke-derived summaries and use the plan's standard-normal synthetic MRI signal fallback.
