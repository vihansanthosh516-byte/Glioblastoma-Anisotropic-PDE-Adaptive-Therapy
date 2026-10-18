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
