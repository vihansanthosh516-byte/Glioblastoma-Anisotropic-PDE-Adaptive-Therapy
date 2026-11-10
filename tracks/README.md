# Project map: Tracks A, B, and C

This folder is the navigation layer for the research project. It does not relocate source code, datasets, or outputs: the existing `src/`, `tests/`, `data/`, `output/`, and `paper/` paths remain canonical so older scripts and provenance links keep working.

## Start here

1. Read [`START_HERE.md`](../START_HERE.md) for the latest dated status, data boundaries, and project rules.
2. Read the [integrated research plan master](../paper/trackABC_integrated_research_plan_master_2026-11-07.md) for the complete proposed method. The [submission copy](../paper/trackABC_integrated_research_plan_submission_2026-11-07.md) shortens only the Procedure section.
3. Use the track hubs for the current scope, primary records, code, tests, outputs, and next gates.
4. Read [`shared/HANDOFF_AND_REPRODUCIBILITY.md`](shared/HANDOFF_AND_REPRODUCIBILITY.md) before trying to reproduce or run anything.

## Plan status (checked 2026-11-09)

Labels come from each file's own status line or from the links in `START_HERE.md`. "Historical" means an older project or direction kept for the record. The open-gap record is [`paper/gap_status_2026-11-09.md`](../paper/gap_status_2026-11-09.md).

| Plan file (in `paper/` unless noted) | Status |
|---|---|
| `trackABC_integrated_research_plan_master_2026-11-07.md` | **Active** complete plan (unfrozen; not an evaluation authorization) |
| `trackABC_integrated_research_plan_submission_2026-11-07.md` and `.docx` | **Active** short submission copy |
| `trackABC_integrated_research_plan_revised_2026-11-08.md` | **Status undecided**: linked from no hub; the user must label it |
| `Track A5 Research Plan Science Fair MD Final.md` | Completed and frozen; A5 failed its 5-of-6 rule |
| `trackc_world_benchmark_protocol_draft_2026-11-07.md` | Proposed draft; evaluation stopped (see Track C hub) |
| `trackABC_integrated_research_plan_submission_2026-11-04.md`, `..._submission_short_2026-11-05.md`, `trackABC_integrated_virtual_patient_research_plan_draft_2026-11-04.md` | Superseded by the 2026-11-07 pair |
| `trackABC_execution_plan_2026-10-31.md`, `trackAB_science_and_execution_plan_2026-11-01.md` | Superseded by the 2026-11-07 master plan |
| `trackB_research_plan_draft_2026-11-03.md`, `track_b_next_phase_plan_2026-11-02.md` | Proposed drafts; no Duke protocol is frozen |
| `track_b_codex_direction_and_execution_plan_2026-11-02.md` | Superseded on 2026-11-03 (own header) |
| `track_b_no_contact_dataset_plan_2026-10-31.md` | Historical LUMIERE exploratory plan |
| `research_plan_trackC_draft.md` | Superseded on 2026-11-03 (own header) |
| `research_plan_trackC_v2_2026-11-04.md`, `trackc_restart_plan_2026-11-03.md` | Superseded by the integrated plan and the 2026-11-06 calibration correction |
| `trackA_alternate_dataset_plan_20261101.md`, `trackA_gse318562_no_contact_plan_20261101.md` | Metadata-only screens; no dataset selected |
| `research_plan_complete_project_form.md`, `research_plan_complete_project_appendix.md`, `research_plan_form_v1_2026-10-22.md`, `research_plan_form_v1_reproduction_appendix_2026-10-22.md`, `research_plan_full_story_2026-10-22.md`, `research_plan_next_phase.md`, `phd_level_plan_probabilistic_rano.md`, `../analysis_plan_v1.md` | Historical |

## Track hubs

- [Track A — expression and regional biology](track-a/README.md)
- [Track B — longitudinal MRI](track-b/README.md)
- [Track C — mathematical simulation and treatment policy](track-c/README.md)
- [Integrated phase — the synthetic bridge](integrated/README.md)

## Shared project rules

- The three tracks do not share patient identifiers. Their connection in the integrated phase is through explicitly assumed synthetic signals, not patient matching.
- Prior patient analyses and prior simulation outcomes stay labeled as exposed development evidence. They are not fresh confirmation.
- A folder name, hash, or successful synthetic check does not establish data-use permission, clinical validity, or independent replication.
- Keep restricted inputs outside Git. Stage named files only; never stage the entire worktree blindly.
