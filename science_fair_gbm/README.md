# GBM Science-Fair Project: File Guide

This folder is the entry point for the main GBM science-fair project. It organizes existing files without moving or rewriting historical code, datasets, or results.

## Start and plan

- [Current project instructions](../START_HERE.md)
- [Track A/B/C improvement backlog](../paper/GBM_TRACKS_IMPROVEMENT_BACKLOG.md)
- [Project status and historical results](../README_STATUS.md)

Prior drafts are preserved, not adopted as the current plan. The forward plan in [the full-story draft](../paper/research_plan_full_story_2026-10-22.md) centers a MU/LUMIERE growth-forecast study and does not set roles for all three tracks. The [form-fields draft](../paper/research_plan_form_v1_2026-10-22.md) and its [reproduction appendix](../paper/research_plan_form_v1_reproduction_appendix_2026-10-22.md) use a probabilistic RANO-derived framing. Use the improvement backlog when preparing the new GBM plan.

## Evidence and data

- [Findings ledger](../vault/PAPER_FINDINGS_LEDGER.md): result values and links to output files.
- [Data-source inventory](../vault/PAPER_DATA_SOURCES.md): dataset roles, recorded cohort counts, and known lineage gaps.
- [Old-track audit](../paper/codex_audit/old_tracks.md): documented leakage, overclaims, simulation limits, and code review boundaries.
- [Audit summary](../paper/codex_audit2/SUMMARY.md): later repository audit summary.
- [Failed hypotheses](../vault/failed_hypotheses.md): negative findings and prior tests.

## Historical Track A/B/C records

- Track A results: [archive_old_gbm/TRACK_A_RESULTS.md](../archive_old_gbm/TRACK_A_RESULTS.md) and [archive_old_gbm/vault/TRACK_A_RESULTS.md](../archive_old_gbm/vault/TRACK_A_RESULTS.md).
- Track B results: [archive_old_gbm/vault/TRACK_B_RESULTS.md](../archive_old_gbm/vault/TRACK_B_RESULTS.md).
- Track C results: [archive_old_gbm/vault/TRACK_C_RESULTS.md](../archive_old_gbm/vault/TRACK_C_RESULTS.md).
- Combined summary: [archive_old_gbm/paper/summary_all_tracks.md](../archive_old_gbm/paper/summary_all_tracks.md).
- Prior digital-twin plan: [archive_old_gbm/paper/project_plan_trusted_twin.md](../archive_old_gbm/paper/project_plan_trusted_twin.md). This is an unrun historical proposal, not the current plan.

## Code and datasets

Track code and outputs remain spread across src/, archive/track_a/, archive/track_b/, output/, and data/. Their presence does not show that a script is safe to run or that a dataset is eligible for a new analysis. Check the current plan, license, manifest, and ledger first.

## Project boundary

Tracks A, B, and C remain within the GBM project. Existing records do not establish patient-level links across their datasets. Keep molecular prediction, MRI forecasting, and treatment simulations distinct unless future records prove a valid join and outcome.
