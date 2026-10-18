# P4 declaration: controller-threshold sweep of the 61-patient adaptive-therapy simulation
Written before any P4 run. Date: 2026-10-18. Status: EXPLORATORY, descriptive, no p-values. Plan: FIX_PLAN_v2_1_approved.md (reviewed twice).

Code: src/44_adaptive_therapy.py, unchanged. Run from the repo with PYTHONPATH=src, PYTHONDONTWRITEBYTECODE=1, output to the workspace.
Fixed: E_MAX ratio 1000; resistance fitness cost 0.25; horizon 5000 steps (500 days); the 61 positive-rho patients; hard-coded restart clause (mass > 0.9 x peak).
Code reading, before running: stop dosing when mass < OFF x baseline OR mass < OFF x peak; restart when mass > ON x baseline OR mass > 0.9 x peak. Because peak >= baseline, the 0.9 x peak clause is expected to fire first, so ON is expected to be near-inert.

Design: THRESHOLD_OFF in {0.50, 0.65, 0.80, 0.90} x THRESHOLD_ON in {1.0, 1.2} = 8 cells. Cell (0.80, 1.0) is the default and was reproduced exactly earlier (max abs diff 0.0). 7 new runs.
Deterministic simulation: no replicates; spread is reported across the 61 patients.

Reported per cell: MTD more resistant / tied / adaptive more resistant (patients); median, IQR and range of per-patient VR (adaptive / MTD final mass) and of drug % of MTD; TTP adaptive earlier / later / same; earlier counts among rho > 0.02 (n=10) and rho <= 0.02 (n=51); fraction of time drug is on; both-censored count. Differences from the default cell are reported for every cell.

Pre-specified interpretation rule (fixed before any run):
C-A 'adaptive never selects more resistance': robust only if adaptive-more-resistant = 0 in all 8 cells; otherwise sensitive, naming the cells.
C-B 'dose-sparing': robust only if mean drug % of MTD is below 50 in all 8 cells (default 32.4). Convention, not a clinical standard; passing it is a modest form of robust.
C-C 'cost' (VR): described with numbers only; no label.
C-D 'progression pattern' (10 of 10 earlier above rho 0.02, 0 of 51 below): robust only if both counts equal the default in all 8 cells; otherwise sensitive, naming the cells.
C-E 'ON is inert' (declared expectation): for each OFF value, ON = 1.2 and ON = 1.0 give identical per-patient metrics. Any difference is reported with the patients affected.
Labels are not to be read without the deviation magnitudes.
All 8 cells are reported. No cell is dropped, and no threshold or label above is changed after the results are seen.
