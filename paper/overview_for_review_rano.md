# Overview for an independent rating (probabilistic RANO), 2026-10-18

You know nothing about this project. Rate it honestly. I want hard feedback, not praise. Every number below comes from a file named in the last column; interim numbers are marked.

## What it is
A single-author science-fair project (Greater Austin regional fair, February 2027) in computational biology of glioblastoma. Question: RANO 2.0 calls tumour progression at a volume rise of 40% or more (partial response at a fall of 65% or more; progression is judged against the nadir; the 0.5 cm3 floor is our own rule, not part of RANO 2.0). NOTE: Result 1 compares CONSECUTIVE scans, not the nadir, so it is not the RANO rule itself. Volumes now come from automated segmentation. Are those lines above the measurement noise, and what does noise do to trial endpoints (progression-free survival, PFS)? All work is retrospective analysis of public data or simulation. Nothing was tested on live patients. No paper yet.

## Method rules
Each test is written into `analysis_plan_v1.md` and committed before the code runs; the code is committed before the run; each runs once; no re-runs with other settings.

## Result 1 (A15, real data, exploratory)
LUMIERE, two automated tools (HD-GLIO-AUTO, DeepBraTumIA). 430 consecutive post-operative scan pairs; 309 (71.9%) excluded (a volume missing, zero or below 0.5 cm3); 121 eligible pairs, 43 patients. The RANO class differs between the tools in 54.5% of pairs (95% CI 45.6% to 63.8%). Against expert ratings (106 pairs) sensitivity for expert progression is 0.55 and 0.51, specificity 0.61 and 0.74. Minimal detectable change: smallest size band inside the noise; the two larger bands undetermined because the noise lower bound is degenerate. Source: `output/rano_noise/rano_tool_dependence.json`, ledger 2br.
Weakest points: the two tools use different label definitions and are not a scan-rescan pair, so their disagreement overstates noise; the lower noise bound (expert correction of an AI mask) understates it and is near zero for larger tumours; most pairs were lost; one centre; expert ratings are relative to baseline or nadir.

## Result 2 (A16, simulation, assumed parameters)
Two-arm trial simulation, 96 cells, 2000 trials each. Noise in the progression call moves the control median PFS earlier by 15 to 23 weeks (no-noise median 38.9 weeks) and cuts power (N=100, 8-week scans: 0.97 with no noise, 0.29 at change-SD 0.25, 0.05 at 1.0). Two-scan confirmation reduces but does not remove the bias. Source: `output/trial_sim/trial_sim_summary.json`, ledger 2bs.
Weakest points: effect sizes rest entirely on an assumed noise model (independent per-scan noise, nadir rule); noise at 0.5 already equals the +40% line; real noise is probably correlated and reviewed centrally; growth, spread and the 8-month median are assumed.

## Negative or failed findings
- Declared prediction H-S3 (confirmation delays PFS at zero noise) NOT supported: with no noise and monotone growth confirmation changes nothing; the prediction was flawed.
- Declared prediction H-S4 (type I error within 0.04 to 0.06 in all cells) NOT supported: 4 of 48 cells at 0.039 to 0.0605, inside the Monte Carlo error but outside the rule.
- The older GBM platform in this repository (growth forecasting, RL therapy) did not beat a no-change baseline on real data; it is archived in `archive_old_gbm/`.

## Not done
Identifiability proof, calculator, patient timeline viewer, paper, poster, clinician review. Novelty not verified (about 15 snippet-level searches).

## Please rate
1. Is the question real and important? 2. Is the evidence enough for a regional fair? For ISEF? 3. What is the biggest flaw? 4. What would a judge ask that I cannot answer? 5. What one thing should I do next?
