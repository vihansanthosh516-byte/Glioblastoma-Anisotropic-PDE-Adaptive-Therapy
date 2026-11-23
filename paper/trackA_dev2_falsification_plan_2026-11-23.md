# Track A development falsification 2 (frozen before run)

2026-11-23. Generated data only. No real expression, labels or viability opened. Code: `src/trackA_dev2.py`. Output: `output/trackA_dev2/result.json` (new file; earlier outputs kept). Seed root 202611232. 20 fresh replicate draws per scenario.

## Regional: harm-calibrated override (candidate `override_calibrated`)

- Override is proposed for a pair when both reconstructed shared margins agree in sign and that sign differs from the nonspatial decision. Selection score = smaller absolute reconstructed margin.
- Threshold grid fixed: 0, 0.25, 0.5, ..., 3.0. On calibration split A choose the smallest threshold with at least 5 proposed overrides whose one-sided 95% Clopper-Pearson upper bound on the wrong-override rate is at most 0.35. Validate on independent split B with the same rule. If either step fails, no overrides (full fallback).
- Ablation `override_unchecked`: same overrides at threshold 0 without the check. Also reported: nonspatial, each spatial arm, molecule spatial, agreement only, prior bounded model is not refit here (kept in run 1).
- Scenarios: clean; shared error (both reconstructions shifted -1.8 toward the wrong side); divergent error; calibration mismatch (calibration clean, evaluation shared error); density-only (spatial carries no signal beyond density); missed contamination (shared error on 30% of pairs, calibration has none).
- Per-patient record: baseline, both reconstructions, final decision, changed, harmed.
- Pass rule (development only): candidate mean concordance minus strongest noncandidate arm > 0 with the 95% interval over replicates above 0, in clean; candidate harms per enrolled patient <= 0.02 in every adverse scenario; candidate beats `override_unchecked` in at least one adverse scenario. Otherwise the candidate fails.

## TMZ: replicated within-origin direction (candidate `replicated_direction`)

- Base = origin-balanced ridge (same penalty 0.2, same standardisation).
- Inside training: direction d = ridge on within-origin centered expression and viability (repeated origins only). Split training origins into two fixed halves. Admit d only if (a) cosine(d_A, d_B) > 0 and (b) adding the projection x.d_A to the base improves origin-mean MAE on half B, and the mirror holds. If admitted, final = origin-balanced ridge on [x, x.d]. If not, final = base.
- Ablations on the same splits: base; same with d from shuffled sister links (within stratum, whole rows kept); outcome-blind direction (first principal axis of within-origin centered expression); the existing one-slope contrast at weight 1.
- Scenarios: transfer (between-origin effect = within-origin effect); reversal (between = -within); unequal sisters (1 to 8 cultures); assay shift (evaluation expression feature mixing); null (no effect).
- Primary: per-origin absolute error of mean viability; paired difference base minus candidate.
- Pass rule (development only): in transfer, candidate beats base and the shuffled-link arm, with the replicate 95% interval of the paired difference above 0; in reversal and null, candidate is not worse than base by more than 0.01 MAE. A tie with base or with shuffled links means the sister-link mechanism is unproven.

Thresholds here are software-development settings, not scientific gates for an external test.
