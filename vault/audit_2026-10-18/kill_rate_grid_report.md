# Kill-rate sensitivity grid, 61-patient adaptive-therapy simulation

Status: exploratory sensitivity analysis of the script 44 simulation. Inputs are the project's own MU-Glioma-Post fitted growth rates; every number here is a simulation with an assumed kill-rate scale. Nothing in this analysis is a locked test.

## What was run
- Script `src/44_adaptive_therapy.py`, unchanged, run from the repo with `ADAPTIVE_E_MAX_RATIO` set and `ADAPTIVE_OUT_DIR` pointing to the workspace. Nothing was written into the repo.
- New ratios: 100, 250, 4000, 8000, plus a full re-run of 1000. Existing ratios 500, 1000, 2000 were read from `output/emax_sweep/` (ledger 2n).
- Grid and rules were written before the new runs (see the declaration below). Provenance note (2026-10-19 audit correction): that ordering, the "unchanged" code and the output location are author-reported; no timestamped commit or run log is retained here, and run duration and exit status are not recorded. The numbers below match `kill_rate_grid_results.csv`.

## Reproduction check (declared in advance): passed
The full ratio-1000 re-run matches the committed output exactly: maximum absolute difference 0.0 over 28 numeric metrics for all 61 patients; summary 10 earlier, 22 any progression, 39 both censored, drug 32.3551% of MTD, mean VR 2.3913. The new ratios are therefore comparable with the existing 500, 1000 and 2000 runs.

## Results (61 patients, same cohort at every ratio)

| Ratio | MTD more resistant | Tied | Adaptive more resistant | Mean drug, % of MTD | Mean final-mass ratio (adaptive / MTD) | TTP adaptive earlier / later / same | Of 10 patients with rho > 0.02: earlier | Of 51 with rho <= 0.02: earlier |
|---|---|---|---|---|---|---|---|---|
| 100 | 27 | 34 | 0 | 96.0 | 1.00 | 27 / 0 / 34 | 10 | 17 |
| 250 | 59 | 2 | 0 | 64.5 | 1.24 | 14 / 7 / 40 | 10 | 4 |
| 500 | 57 | 4 | 0 | 44.9 | 1.73 | 10 / 12 / 39 | 10 | 0 |
| 1000 | 54 | 7 | 0 | 32.4 | 2.39 | 10 / 12 / 39 | 10 | 0 |
| 2000 | 51 | 10 | 0 | 26.3 | 2.89 | 10 / 12 / 39 | 10 | 0 |
| 4000 | 49 | 12 | 0 | 23.3 | 3.08 | 9 / 13 / 39 | 9 | 0 |
| 8000 | 43 | 18 | 0 | 22.1 | 3.17 | 5 / 13 / 43 | 5 | 0 |

Median resistant fraction, MTD / adaptive: 100: 0.24 / 0.24; 250: 0.51 / 0.34; 500: 0.85 / 0.31; 1000: 0.99 / 0.31; 2000: 1.00 / 0.31; 4000: 1.00 / 0.31; 8000: 1.00 / 0.35.

## Reading
1. **Resistance direction is robust.** Adaptive therapy ends with more resistant cells than MTD in 0 of 61 patients at all 7 ratios. This is the one result that holds across the whole grid.
2. **Its size is not.** At ratio 100 the drug barely acts: the median resistant fraction is identical in both arms (0.24), 34 of 61 patients are tied, and adaptive therapy uses 96% of the MTD dose. From 250 upward, 59, 57, 54, 51, 49, 43 patients have more resistance under MTD; ties grow (2 to 18) as the MTD arm saturates.
3. **The cost grows with the kill scale.** Mean final tumour mass under adaptive therapy relative to MTD rises from 1.00 (equal) at 100 to 3.17 at 8000, while drug use falls from 96% to 22% of MTD. The trade-off (less drug and less resistance, more tumour) is present at every ratio above 100.
4. **The time-to-progression result is specific to 500-2000.** The ledger's pattern (all 10 patients with rho > 0.02 progress earlier under adaptive therapy, 0 of 51 below) is identical at 500, 1000 and 2000. It changes outside that range:
   - At 4000 and 8000, 9 and 5 of the 10 high-rho patients progress earlier.
   - At 250 and 100, 4 and 17 of the 51 low-rho patients progress earlier. At ratio 100, 27 patients progress earlier under adaptive therapy and 0 later.
   So "adaptive progresses earlier only in fast growers" should be stated for the 500-2000 range, not as a general result.
5. **No ratio was dropped and no threshold was changed after seeing results.**

## What this does and does not show
- It shows the resistance and dose-sparing findings do not depend on the exact kill scale within 250-8000, and that the cost does depend on it.
- It does not validate any kill rate. E_max = ratio x patient growth rate is an assumption; real per-patient drug effect was not measured. This stays a simulation result.
- Not varied: fitness cost 0.25, THRESHOLD_OFF 0.80 and ON 1.0, 500-day horizon (39 of 61 patients are censored in both arms at every ratio >= 250), the 61 positive-rho patients (42 negative-rho responders excluded).
- One more axis worth declaring before running it: the controller thresholds (ADAPTIVE_THRESHOLD_OFF / ON), since the 0.80 setpoint pins adaptive nadirs by design.

## Declaration written before the new runs
# Kill-rate sensitivity grid: declaration (written before any new run)
Date: 2026-10-17. Status: EXPLORATORY sensitivity analysis of the simulation in script 44. Not a locked test.

Question: does the 61-patient adaptive-therapy result depend on the assumed kill scale E_MAX_RATIO (E_max = ratio x patient rho)?

Existing runs (read only, from output/emax_sweep/): ratio 500, 1000, 2000 (ledger 2n).
New runs (script 44 unchanged, env var ADAPTIVE_E_MAX_RATIO, output to workspace): ratio 100, 250, 4000, 8000.
Reproduction check: ratio 1000 re-run in this environment must match output/adaptive_cohort_summary.json
  (10 earlier, 22 any progression, 39 both censored, drug 32.36% of MTD, VR 2.391). If it does not, new ratios are reported as
  'not comparable' and the environment difference is described instead.
Metrics per ratio (same definitions as script 92): n MTD more resistant / adaptive more resistant / ties; median resistant fraction per arm;
  mean drug % of MTD; mean final-mass ratio (adaptive/MTD); TTP adaptive earlier / later / same; both-censored; rho>0.02 split.
Rules: no ratio is dropped, no threshold or conclusion is changed after seeing results; ratios where the result weakens or reverses are reported as such.
Not varied: fitness cost 0.25, THRESHOLD_OFF 0.80, ON 1.0, horizon 5000 steps, cohort of 61 positive-rho patients.

