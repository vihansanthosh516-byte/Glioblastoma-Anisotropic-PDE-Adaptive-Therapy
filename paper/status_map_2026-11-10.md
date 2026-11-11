# Status map and backlog (2026-11-10)

Built from disk at commit e58c9d1 on `main`. No research job was running (only MCP server processes). Authoritative gap record stays `paper/gap_status_2026-11-09.md`; this file adds what changed since and a ranked backlog.

## State

| Item | State | Evidence |
|---|---|---|
| Track A A5 | Completed, failed frozen 5-of-6 rule (4 of 6) | `output/trackA_yu2020/result.json`; gap A1 |
| Track A external cohort | Blocked, none passes gates | gap A5 |
| Track B Duke 10-mm analysis | Completed, exploratory, cohort spent | gap B1 (private records) |
| Track B terms / school record | Blocked on user (U1-U3) | gap B6 |
| Track C confirmatory eval | Stopped; 14-day target infeasible at 7 worlds/family | gap C4; `paper/trackc_ge_betting_coverage_audit_result_2026-11-09.md` |
| GE-betting interval | Candidate, not adopted, no independent review | gap C3 |
| Science-corrected plan | Proposed, not adopted | `paper/trackABC_integrated_research_plan_submission_science_revision_2026-11-09.md` |
| PF prior vs generator support | **New check, 2026-11-10.** Gap confirmed | below |

## New finding: particle-filter prior does not cover the generator (outcome-blind)

- Analytic generator support (from `src/trackc_integrated_sampler.py` world axes): r about 0.00050 to 0.123 per day; B0 0.02 to 0.65.
- Prior in `configs/trackABC_integrated_v2.json` and `src/trackc_sim.py` `RANGES`: r 0.001 to 0.05; B0 0.05 to 0.5. m, c, epsilon are covered.
- Synthetic draw (100 worlds per family, 12 patients each, 2400 patients, seeds 9,200,000+ / 9,300,000+): 246 patients (10.3%) outside the r prior, 324 (13.5%) outside the B0 prior.
- Pinned as two strict expected failures in `tests/test_trackc_prior_support.py`.
- Weakest point: support is analytic for the declared axes; no policy run measured the effect on utility, and none should before the protocol decides.
- Changing the prior changes a config setting. It needs a dated amendment and a user decision. Not changed.

## Ranked backlog (safe = no patient data, no policy outcome)

1. Decide the PF prior (user decision; amendment drafted only on approval).
2. Track C U7 decisions: 81 vs 135 models; 112-day forecast length (user).
3. Runtime config loader for the Track C driver (code; safe; blocked in value until C4 is resolved).
4. Independent statistical review of GE-betting (needs a reviewer outside this session).
5. CGGA MGMT association (C1 in proposed plan): check data terms and current code status before any run.
6. Documents only the user can supply: Duke terms, school approval, token revocation (U1-U4).

## Checks run

Scratch copy of `src/`, `tests/`, `configs/`; source vs scratch `src/` and new test identical (`diff -rq`).
`PYTHONPATH=<scratch>/src py -3.14 -m pytest tests/test_trackc_prior_support.py tests/test_trackc_integrated_sampler.py tests/test_trackc_config_parity_extended.py -q -s -p no:cacheprovider` → 22 passed, 5 expected failures, 2.11 s.
