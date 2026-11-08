# Project map: Tracks A, B, and C

This folder is the navigation layer for the research project. It does not relocate source code, datasets, or outputs: the existing `src/`, `tests/`, `data/`, `output/`, and `paper/` paths remain canonical so older scripts and provenance links keep working.

## Start here

1. Read [`START_HERE.md`](../START_HERE.md) for the latest dated status, data boundaries, and project rules.
2. Read the [integrated research plan master](../paper/trackABC_integrated_research_plan_master_2026-11-07.md) for the complete proposed method. The [submission copy](../paper/trackABC_integrated_research_plan_submission_2026-11-07.md) shortens only the Procedure section.
3. Use the track hubs for the current scope, primary records, code, tests, outputs, and next gates.
4. Read [`shared/HANDOFF_AND_REPRODUCIBILITY.md`](shared/HANDOFF_AND_REPRODUCIBILITY.md) before trying to reproduce or run anything.

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
