# Script 83: plastic resistance and adaptive therapy

Model: script 76 TwoPopSim plus drug-induced S<->T switching (sig_i on drug, sig_r off drug). f_r0 = 0.001, cost 0.25, eps 0, window 365 d.
Validation: sigma = 0 reproduces script 76 control (max RMST diff 0.0 d).
Rule (pre-specified in script header): positive iff >= 3 contiguous of 6 primary cells (real_test, seed59 occupancy) have AT - paced CI lower > 0 and plasticity-attributable days > 0 for the same arm.

Result (real_test, seed59): 0 of 6 cells for AT50, 0 of 6 for AT80. `positive_overall` = false.
- AT50 - paced: -85 to -143 d, all CIs below 0.
- AT80 - paced: -26 to +1.5 d; only sig_i 0.005 / sig_r 0.05 has CI lower bound 0.0 (+1.3 d), with plasticity-attributable -0.6 d.
- Plasticity-attributable AT50 days are positive (+0.1 to +57 d) but AT50 stays far below paced, so criterion (a) fails.

Verdict: C3 stays negative under fixed clones (76), MGMT-calibrated clones (82) and plastic resistance (83).
Data: output/plastic_resistance/results.json, cells/.
