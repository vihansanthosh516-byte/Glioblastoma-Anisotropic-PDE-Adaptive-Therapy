# Script 82: MGMT-calibrated resistance analysis (C3 on real data)

**Date:** 2026-10-04. **Script / outputs:** `src/82_mgmt_calibrated_resistance.py`, `output/mgmt_resistance/{results.json, km_mgmt_ttp.png, cells/}`. Hypothesis, endpoints, calibration rule and the "flips positive" rule are in the script header (written before the run). A smoke test of the clinical part (no design change) was run first.

## Primary: does MGMT status predict time to progression?
Cohort: Primary Diagnosis = GBM, MGMT in {0, 1}: 79 (MGMT 0) + 51 (MGMT 1) = 130 analysed, 4 excluded for no usable time. Events: 67 and 41. Non-progressors censored at the later of last MRI day and death day.

| | MGMT = 0 (unmethylated) | MGMT = 1 (methylated) |
|---|---|---|
| Kaplan-Meier median TTP, days (95% CI) | **169** (147-196) | **181** (108-207) |
| Median among progressors only | 161 | 152 |
| RMST to day 365 | 178.5 | 187.5 |

- **Log-rank p = 0.60.** Cox HR (MGMT 0 vs 1) 1.11 (95% CI 0.75-1.66, p = 0.59); adjusted for age 1.29 (0.85-1.97, p = 0.24, n = 130).
- RMST_365 ratio MGMT 0 / MGMT 1 = 0.95 (bootstrap CI 0.76-1.20).
- **The hypothesis is not supported:** the direction is as predicted (169 vs 181) but the difference is within noise. The earlier 167 vs 182 figures were all-grade, progressors-only medians; the pre-specified censoring-aware GBM-only values are 169 vs 181.

## Secondary: calibrated resistance and adaptive therapy
- Model RMST_365 under the Stupp schedule across f_r0 = 0.01 ... 0.6 changes by only 117.1 -> 114.8 days (ratio 1.00 -> 0.98). **The data cannot constrain f_r0**: the observed ratio 0.95 lies below the model's whole range, and the log-rank test is null. By the pre-specified clamp rule, the point estimate and lower-CI calibration therefore take the most resistant grid value, f_r0 = 0.6 (a clamp, not an estimate); the upper-CI calibration (ratio 1.20 > 1) equals the anchor 0.01 and is identical to script 76's primary cell. (Model TTP is the RANO criterion, ~117 d under Stupp; observed clinical TTP is ~180 d, so the two are not on the same definition.)
- Primary-cell settings (cost 0.25, eps 0, script-59 seed, 365-day window), AT arm minus same-window paced arm, days (95% CI) and resistance-attributable days (difference-in-differences vs the f_r0 = 0 control):

| f_r0 | Set | Arm | vs paced | Resistance-attributable | Median resistant fraction day 365 |
|---|---|---|---|---|---|
| 0.01 (script 76 primary) | real_test | AT50 | -141.7 (-192, -93) | +1.0 | 0.02 |
| | real_test | AT80 | -19.7 (-50, +2) | -21.6 | 0.02 |
| | cohort64 | AT50 / AT80 | +12.7 / +90.5 | +0.6 / -2.9 | 0.02 |
| 0.6 (clamped) | real_test | AT50 | **-35.0 (-58, -14)** | +107.8 | 0.99 |
| | real_test | AT80 | -6.6 (-15, +1) | -8.5 | 0.96 |
| | cohort64 | AT50 | **+22.0 (+20.9, +23.1)** | **+9.9** | 0.73 |
| | cohort64 | AT80 | +22.0 (+20.9, +23.0) | -71.3 | 0.73 |

- **Flips positive by the pre-specified rule? No.** On real_test (the only set with real patients) neither AT50 nor AT80 has a CI lower bound above 0 against the same-window paced arm at f_r0 = 0.6.
- **Where resistance does help adaptive therapy:** at f_r0 = 0.6, AT50 on cohort64 (a synthetic set) gains +22 days over paced, of which +9.9 days are resistance-attributable. This is the same corner as script 76 (resistance-driven gain in 2 of 384 rows): a very high resistant fraction. It is not supported by the clinical data (MGMT 0 does not progress faster) and is not on real patients.
- Compared with script 76: at the script-76 primary f_r0 the two analyses are identical by construction. The new calibration reproduces script 76's synthetic corner but adds no real-patient evidence for it.

## Verdict
C3 stays negative on real data, and the reason is now data-grounded: in this cohort MGMT status does not predict progression time, so there is no empirical basis for assigning a large resistant fraction to MGMT-unmethylated tumours, and the model's Stupp-TTP is insensitive to f_r0, so TTP data cannot calibrate it either. Report as a limitation: resistance-driven adaptive gains appear only at resistant fractions the clinical data do not support.

Caveats: model TTP (RANO volume criterion) and clinical TTP are different endpoints; observed RMST ratio CI is wide; the calibration is non-identifiable (clamped); only 21 real patients in `real_test`; MGMT is a coarse resistance proxy and resistance in GBM is partly plastic, which the fixed-clone Lotka-Volterra model does not represent.
