# P4: controller-threshold sweep of the 61-patient adaptive-therapy simulation

Status: EXPLORATORY, descriptive, no p-values. Simulation with an ASSUMED kill-rate scale (E_MAX ratio 1000, fixed). Not confirmatory.
Declared before any run: artifact `P4_threshold_sweep_declaration.md`, version id de0f987e-7125-46d6-baab-f826d425b53f, SHA-256 8febbe13fd18d2d09db86311b47b3d6c817a2e6bc5a44cc3376d4e25d2cda286. Design, metrics and interpretation rule below are exactly as declared.
Code: `src/44_adaptive_therapy.py`, unchanged, run from the repo with output to the workspace. Seven new runs (all exit code 0, about 38 min each in parallel); the cell OFF 0.80 / ON 1.0 is the default and comes from the committed `output/emax_sweep/r1000`, which a full re-run reproduced exactly earlier (max abs difference 0.0). The simulation is deterministic, so there are no replicates; spread is across the 61 patients.

![Threshold sweep](P4_threshold_sweep.png)

## Results per cell (61 patients; 10 with rho > 0.02, 51 with rho <= 0.02)
| OFF | ON | MTD more resistant | Tied | Adaptive more resistant | Mean VR (median, IQR) | Mean drug % of MTD | Drug-on fraction | TTP adaptive earlier / later / same | Earlier: 10 fast / 51 slower | Mean holidays |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.50 | 1.0 | 50 | 11 | 0 | 1.95 (1.25, 0.97 to 3.05) | 43.3 | 0.43 | 7 / 15 / 39 | 7 / 0 | 5.8 |
| 0.50 | 1.2 | 50 | 11 | 0 | 1.95 (1.25, 0.96 to 3.05) | 43.1 | 0.43 | 7 / 15 / 39 | 7 / 0 | 6.5 |
| 0.65 | 1.0 | 53 | 8 | 0 | 2.20 (1.83, 1.00 to 3.51) | 36.9 | 0.38 | 10 / 12 / 39 | 10 / 0 | 20.4 |
| 0.65 | 1.2 | 53 | 8 | 0 | 2.20 (1.83, 1.00 to 3.51) | 36.2 | 0.38 | 10 / 12 / 39 | 10 / 0 | 18.7 |
| 0.80 | 1.0 | 54 | 7 | 0 | 2.39 (1.99, 1.00 to 3.83) | 32.4 | 0.41 | 10 / 12 / 39 | 10 / 0 | 48.9 |
| 0.80 | 1.2 | 54 | 7 | 0 | 2.39 (1.99, 1.00 to 3.83) | 31.8 | 0.40 | 13 / 9 / 39 | 10 / 3 | 48.4 |
| 0.90 | 1.0 | 54 | 7 | 0 | 2.50 (1.99, 1.00 to 4.18) | 29.3 | 0.77 | 10 / 12 / 39 | 10 / 0 | 200.2 |
| 0.90 | 1.2 | 54 | 7 | 0 | 2.50 (1.99, 1.00 to 4.18) | 27.4 | 0.77 | 13 / 9 / 39 | 10 / 3 | 181.1 |

## Difference from the default cell (OFF 0.80, ON 1.0)
| OFF | ON | Δ MTD more resistant | Δ tied | Δ mean VR | Δ mean drug % | Δ drug-on fraction | Δ TTP earlier | Δ TTP later | Δ earlier (10 fast) | Δ earlier (51 slower) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.50 | 1.0 | -4 | +4 | -0.438 | +10.9 | +0.021 | -3 | +3 | -3 | +0 |
| 0.50 | 1.2 | -4 | +4 | -0.438 | +10.8 | +0.020 | -3 | +3 | -3 | +0 |
| 0.65 | 1.0 | -1 | +1 | -0.187 | +4.5 | -0.027 | +0 | +0 | +0 | +0 |
| 0.65 | 1.2 | -1 | +1 | -0.189 | +3.8 | -0.029 | +0 | +0 | +0 | +0 |
| 0.80 | 1.0 | +0 | +0 | +0.000 | +0.0 | +0.000 | +0 | +0 | +0 | +0 |
| 0.80 | 1.2 | +0 | +0 | -0.001 | -0.6 | -0.008 | +3 | -3 | +0 | +3 |
| 0.90 | 1.0 | +0 | +0 | +0.109 | -3.1 | +0.362 | +0 | +0 | +0 | +0 |
| 0.90 | 1.2 | +0 | +0 | +0.111 | -5.0 | +0.360 | +3 | -3 | +0 | +3 |

## Pre-specified labels
- **C-A: adaptive never selects more resistance. ROBUST.** Adaptive more resistant = 0 of 61 in all 8 cells. MTD more resistant ranges from 50 (OFF 0.50) to 54; ties from 7 to 11.
- **C-B: dose-sparing (mean drug below 50% of MTD in all cells). ROBUST, modestly.** Mean drug use is 27.4% to 43.3% across cells (default 32.4%). The 50% cutoff was a convention, so this is a modest form of robust. Note OFF 0.90 gives drug-on fraction 0.77 with 200 mean holidays: this is consistent with the controller chattering around the threshold, because it stops below 0.9 x baseline and restarts above 0.9 x peak (peak is close to baseline for most patients).
- **C-C: cost (final-mass ratio), numbers only.** Mean VR rises from 1.95 (OFF 0.50) to 2.50 (OFF 0.90); median from 1.25 to 1.99. The adaptive arm always ends with a larger tumour than MTD on average. Changing ON from 1.0 to 1.2 moves the mean VR by at most 0.003.
- **C-D: progression pattern (10 of 10 fast growers earlier, 0 of 51 slower earlier). SENSITIVE.** Deviating cells:
  - OFF 0.50, both ON values: 7 of 10 fast growers progress earlier (default 10); adaptive earlier / later / same = 7 / 15 / 39 (default 10 / 12 / 39).
  - OFF 0.80 and OFF 0.90 with ON 1.2: 3 of 51 slower growers progress earlier (default 0); earlier / later / same = 13 / 9 / 39 (default 10 / 12 / 39).
  - OFF 0.65 with either ON, and OFF 0.80 / 0.90 with ON 1.0, match the default pattern.
- **C-E: 'ON is inert' (declared expectation). NOT MET as stated.** ON changes the results for 12, 16, 16 and 26 patients at OFF 0.50, 0.65, 0.80 and 0.90.

| OFF | Patients with any difference (of 61) | TTP differs | VR differs (max abs ΔVR) | Drug % differs (max abs Δ) | Holidays differ | Lowest rho among affected | Affected with rho > 0.02 |
|---|---|---|---|---|---|---|---|
| 0.50 | 12 | 5 | 12 (0.004) | 10 (6.8) | 7 | 0.0121 | 7 |
| 0.65 | 16 | 6 | 16 (0.026) | 16 (8.0) | 14 | 0.0112 | 10 |
| 0.80 | 16 | 9 | 16 (0.022) | 16 (5.3) | 15 | 0.0125 | 10 |
| 0.90 | 26 | 14 | 26 (0.151) | 26 (8.4) | 26 | 0.0068 | 10 |

The code reading predicted that ON would matter only when a patient's peak mass exceeds about 1.11 x baseline (the 0.9 x peak restart clause otherwise fires first). That fits who is affected: the lowest rho among affected patients is 0.0068, so only faster-growing patients change. The strict expectation of identical per-patient metrics was wrong; ON is inert for most but not all patients. Largest TTP changes (adaptive TTP, ON 1.0 to ON 1.2, in steps): OFF 0.80: PatientID_0099 2654 to 824, PatientID_0136 3031 to 1388, PatientID_0260 3429 to 1679, PatientID_0076 820 to 458; OFF 0.90: PatientID_0136 3031 to 1102, PatientID_0260 3434 to 1679.

## What this does and does not show
- Within the 8 cells, the resistance direction and the dose saving do not depend on the controller thresholds; the size of the final-tumour cost does (1.95 to 2.50).
- The ledger's progression pattern (10 of 10 earlier above rho 0.02, 0 of 51 below) is specific to OFF 0.65 to 0.90 with ON 1.0, as well as to E_MAX ratios 500 to 2000 (entry 2bi). Quote it with those conditions.
- It does not validate any threshold or any kill rate. Not varied: fitness cost 0.25, the hard-coded 0.9 x peak restart clause, the 500-day horizon (39 of 61 patients censored in both arms in every cell), the cohort of 61 positive-rho patients.
- Eight cells is a small design: it cannot separate OFF and ON effects from interactions beyond what is tabulated.
