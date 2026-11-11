# Track C prior-recovery diagnostic: result (2026-11-10, report-only)

Protocol: `paper/trackc_amendment_pf_prior_support_2026-11-10.md`, revision 3, section 4 (committed in e41f515 before the run). Code at commit 04c620f plus `src/221_trackc_prior_recovery_diagnostic.py` (identical content, first named `219_...`). Output: `output/trackc_prior_recovery_diagnostic_20261110/diagnostic.json` (sha256 aa188bf2de663b00d1c2b3d785f5c9d174069a6bf236fa1b6abf55c206a9367f) and `per_patient.jsonl` (400 rows). Run in a scratch copy whose `src/` matched the worktree byte for byte; outputs copied unedited. `git_head_at_run` is null in the JSON because the scratch copy has no `.git`. Seeds recorded in `src/trackc_seed_registry.py` (script 221).

**This is a synthetic filter-recovery diagnostic. No controller ran, no utility was computed, no policies were compared. Report-only: no pass threshold. The results did not and cannot change the V3 bounds.**

## Results (400 synthetic patients, 50 worlds per family, 500 particles, 5 scans)

Joint 90%-interval coverage (all six parameters covered):

| Family | Subset | n | Narrow V2 | V3 | Paired V3 − narrow, mean [world-cluster 95% CI] |
|---|---|---|---|---|---|
| single | all | 200 | 0.300 | 0.455 | +0.155 [0.090, 0.220] |
| single | truth outside narrow r/B0 | 47 | 0.000 | 0.447 | — |
| single | truth inside narrow | 153 | 0.392 | 0.458 | — |
| mixture | all | 200 | 0.360 | 0.445 | +0.085 [0.005, 0.170] |
| mixture | truth outside narrow | 44 | 0.000 | 0.614 | — |
| mixture | truth inside narrow | 156 | 0.462 | 0.397 | — |

Marginal coverage, all patients (narrow → V3): single r 0.710 → 0.795, B0 0.735 → 0.870; mixture r 0.730 → 0.865, B0 0.735 → 0.855. m, f0, c, eps changed by at most 0.05 and every one of their paired CIs includes 0.

Precision cost of V3: median log-width of r rose (single 1.206 → 1.510; mixture 1.053 → 1.405; paired mean +0.284 [0.217, 0.348] and +0.332 [0.280, 0.386]). m log-width rose (+0.237 and +0.305, CIs exclude 0). B0 width rose slightly (+0.016, +0.017).

Filter health: median per-patient minimum ESS single 57.4 → 52.3, mixture 51.3 → 54.0 (of 500). Resample events 1315 → 1336. Filter time 30.6 s → 30.8 s total; whole run 66.2 s.

## Reading, with limits

- Under the narrow prior, no patient whose true r or B0 lay outside its bounds had joint coverage (0 of 91). V3 restores coverage for those patients.
- V3 costs precision: wider r and m intervals for all patients. For mixture patients whose truth was inside the narrow bounds, joint coverage was lower under V3 (0.397 vs 0.462); no paired CI was computed for that subset.
- Neither prior reaches nominal 90% coverage (joint 0.30-0.46, marginal mostly 0.71-0.91). The filter is undercovering under both; this was expected and is why the rule is report-only.
- Weakest points: one fixed schedule; 4 patients per world; imaging channel only; the filter's resample-move step has no Metropolis-Hastings correction (review finding C-d), so its intervals are approximate. Says nothing about policy utility.
- Confirmatory Track C evaluation remains stopped (14-day precision gap C4 unchanged).
