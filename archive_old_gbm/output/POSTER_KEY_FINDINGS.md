> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# POSTER KEY FINDINGS — Month 10 Master Cohort Synthesis

All numeric values are sourced from `output/master_cohort_summary.json` and computed at write time. Eight synthetic cases with assumed parameters; numeric sources are statistics/validation fields in that JSON. No clinical validation or non-inferiority is established.

- Anisotropic tensor growth yields fractal invasion fronts (Df 1.04-1.86), significantly higher than the isotropic/spherical baseline (paired t = 11.92, p<0.001, Cohen's d = 4.22).
- Stromal microenvironment coupling maintains tumor-GF front correlation in the range 0.935-0.952 across all 8 patients (hard floor 0.90; all patients clear the floor).
- In this assumed simulation adaptive TTP is much shorter than continuous MTD (means 498.75 vs 1752.5 steps; `statistics.ttp_non_inferiority`), despite lower exposure (paired t = 30.83, p<0.001; TTP ratio mean = 0.277, range 0.194-0.392) at 21-29% lower cumulative drug exposure (mean ± SD = 24.1 ± 2.5%).
- Exploratory inflammatory-score tiers give non-monotonic time-to-progression; the association is not significant into Low/Mid/High tiers (mean TTP MTD = 1950 / 1519 / 1711 steps; Pearson r(infl vs TTP) = -0.50, p=0.2113).
- Drug-toxicity reduction correlates with inflammatory score (Pearson r = 0.53, p=0.1814; Spearman rho = 0.60, p=0.1195), both correlations are nonsignificant in eight simulated cases; no patient-specific benefit established.
- Internal tensor-field and mass-conservation numerical checks pass, not empirical/clinical validation (symmetry residual = 0.0 < 1e-12; relative mass error = 1.7669748230352867e-16).
