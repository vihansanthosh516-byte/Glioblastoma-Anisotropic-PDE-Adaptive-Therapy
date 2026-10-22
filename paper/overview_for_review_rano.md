# Overview for an independent rating (probabilistic RANO), 2026-10-18

You know nothing about this project. Rate it honestly. I want hard feedback, not praise. Each result gives its source; interim numbers are marked. This is a dated review snapshot, updated for audit corrections on 2026-10-19; START_HERE.md governs current status.

## What it is
A single-author science-fair project (Greater Austin regional fair, February 2027) in computational biology of glioblastoma. Question: RANO 2.0 calls tumour progression at a volume rise of 40% or more (partial response at a fall of 65% or more; progression is judged against the nadir; the 0.5 cm3 floor is our own rule, not part of RANO 2.0). NOTE: Result 1 compares CONSECUTIVE scans, not the nadir, so it is not the RANO rule itself. Volumes now come from automated segmentation. Are those lines above the measurement noise, and what does noise do to trial endpoints (progression-free survival, PFS)? All work is retrospective analysis of public data or simulation. Nothing was tested on live patients. No paper yet.

## Scope and statistics
Only the VOLUMETRIC component of RANO 2.0 is modelled (not new lesions, steroids, clinical status or human approval of segmentations). The two-scan confirmation rule is our simplification, not the RANO 2.0 algorithm. Confidence intervals for the real-data results are patient-cluster bootstrap (10,000 draws), because pairs from one patient are not independent. The four A19 predictions are four readings of ONE simulator, not four independent confirmations. The A19 change-SD rests on an ASSUMED conversion from a published repeatability coefficient (the paper gives no formula or scale); Amendment 21 audits it.

## Method rules
Each test is written into `analysis_plan_v1.md` and committed before the code runs; the code is committed before the run; each runs once; no re-runs with other settings.

## Result 1 (A15, real data, exploratory)
LUMIERE, two automated tools (HD-GLIO-AUTO, DeepBraTumIA). 430 consecutive post-operative scan pairs; 309 (71.9%) excluded (a volume missing, zero or below 0.5 cm3); 121 eligible pairs, 43 patients. Consecutive-scan ratio classes inspired by RANO thresholds differ between tools in 54.5% of pairs (95% CI 45.6% to 63.8%); this is not full RANO 2.0 classification. Against expert ratings (106 pairs), sensitivity for expert progression is 0.55 and 0.51, specificity 0.61 and 0.74. Both proxy-based MDC95 estimates exceed ln(1.40) in the smallest size band. In each larger band, proxy estimates fall on opposite sides of the line. Source: `output/rano_noise/rano_tool_dependence.json`, ledger 2br.
Weakest points: neither proxy is a scan-rescan measure or formal bound on pure measurement noise. Cross-tool disagreement assumes independent tool errors for its conversion and can include changing label-definition effects; AI-anchored expert correction is not guaranteed to be a lower bound. Most pairs were lost; one centre; expert ratings are relative to baseline or nadir.

## Result 2 (A16, simulation, assumed parameters)
Two-arm trial simulation, 96 cells, 2000 trials each. Noise in the progression call moves the control median PFS earlier by 15 to 23 weeks (no-noise median 38.9 weeks) and cuts power (N=100, 8-week scans: 0.97 with no noise, 0.29 at change-SD 0.25, 0.05 at 1.0). Two-scan confirmation reduces but does not remove the bias. Source: `output/trial_sim/trial_sim_summary.json`, ledger 2bs.
Weakest points: effect sizes rest entirely on an assumed noise model (independent per-scan noise, nadir rule); noise at 0.5 already equals the +40% line; real noise is probably correlated and reviewed centrally; growth, spread and the 8-month median are assumed.

## Result 3 (A18, sensitivity after the first outside review)
Corrections first: the 0.5 cm3 floor is our own rule (not in RANO 2.0), and Result 1 used consecutive scans, not the nadir. Sensitivity runs: in the simulation the earlier-PFS bias and power loss already appear at a change-SD of 0.10 (single-scan: -4.7 weeks, power change -0.054); correlated noise (AR(1), rho 0.8) cuts the bias by about half to two thirds but does not remove it; adding death keeps the direction. On the real data the tool flip rate is 62.6% with any positive volumes, 58.1% / 54.5% / 52.9% at floors 250 / 500 / 1000 mm3, and 54.4% (CI 37.7% to 66.3%) against each tool's own nadir, but that last analysis keeps only 68 pairs from 20 patients. Eligible pairs are later and more often expert-rated PD than excluded pairs. Source: `output/rano_sensitivity/rano_sensitivity_summary.json`, ledger 2bt. Weakest point: the real noise level in each tool is still unknown, so the simulation shows what noise of a given size would do, not what these tools do.

## Result 4 (A19, simulation anchored to a published study)
The change-SD was taken from a published scan-rescan study (Abu Khalaf 2021, AJNR; 20 recurrent glioblastoma patients scanned about 2 days apart; enhancing-tumour repeatability coefficient 39% to 46%, CI 28% to 67%), converted by an assumed rule s = RC%/1.96 to 0.143 to 0.342. In the same trial model (N=100, 8-week scans) single-scan PFS comes out 7.6 to 20.4 weeks earlier and power at the assumed true effect falls from 0.97 to 0.77 at the mildest anchor and lower at the others (power change -0.20 to -0.82); the two-scan rule roughly halves the bias. All four declared predictions held. Source: `output/anchored_sim/anchored_sim_summary.json`, ledger 2bu. Weakest points: one study, one scanner, large tumours only, a different tool, assumed conversion, and the 2-day gap may include real change and contrast-timing effects.

## Result 5 (A21, audit of the A19 conversion)
The published repeatability coefficient can be turned into a change-SD in three ways because the paper does not say which scale it used: 0.079 to 0.342 overall (linear 0.143 to 0.342, log back-transform 0.126 to 0.262, absolute mL 0.079 to 0.202). Simulated at seven selected change-SDs, not all 17 converted values (N=100, 8-week scans): single-scan PFS is earlier under every reading (point estimates -5.8 to -15.0 weeks) and power falls (-0.07 to -0.63 at the point estimates; -0.70 at s=0.262); the SIZE depends on the reading by about 2.5 times; at the lowest reading the two-scan rule keeps the bias under a month (-0.6 to -3.3 weeks) and only from s=0.168 does it exceed one month. All four declared predictions held on that selected simulation grid. Source: `output/conversion_audit/conversion_audit_summary.json`, ledger 2bx. Not identifiable from the published numbers: the scale, any removal of the later-scan rise, the error shape, volume dependence.

## Negative or failed findings
- Declared prediction H-S3 (confirmation delays PFS at zero noise) NOT supported: with no noise and monotone growth confirmation changes nothing; the prediction was flawed.
- Declared prediction H-S4 (type I error within 0.04 to 0.06 in all cells) NOT supported: 4 of 48 cells at 0.039 to 0.0605, inside the Monte Carlo error but outside the rule.
- Older PDE growth forecasts did not significantly beat persistence on real scans; the population-rate geometric baseline did show a positive gain (ledger 2w). RL therapy results concern simulations with assumed treatment effects, a different task. Older work is archived in `archive_old_gbm/`.

## Not done
Identifiability proof, paper, poster, clinician review. (Calculator `src/rano_calc.py` and timeline viewer `output/viewer_rano/viewer_rano.html` are built; arithmetic and display only, ledger 2bv and 2bw.) Novelty not verified (about 15 snippet-level searches).

## Please rate
1. Is the question real and important? 2. Is the evidence enough for a regional fair? For ISEF? 3. What is the biggest flaw? 4. What would a judge ask that I cannot answer? 5. What one thing should I do next?

## Audit corrections after this dated snapshot (2026-10-21)

- RANO 2.0 progression reference is context-dependent. Depending on setting and scan history, comparison may be to baseline or nadir; the project's running-minimum rule is a simplification. See https://pmc.ncbi.nlm.nih.gov/articles/PMC10860967/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC11493521/ .
- Keep A21 simulation ranges separate: seven A21 grid values give E1 shifts of -2.6 to -15.7 weeks; all 11 simulated A19+A21 values include -20.4 weeks at A19 s=0.342. The -5.8 to -15.0-week range refers to point-estimate conversion readings only, with upper endpoints from A19.
- Conversion sensitivity ratios are 2.05 within BraTumIA 1.2, 2.43 within BraTumIA 2.0, and 2.60 across compared endpoint readings. These describe modeled effect estimates, not the seven-value A21 grid. Source: `output/conversion_audit/conversion_audit_summary.json` and ledger 2bx.
- The 34.8-week simulation median remains assumed. The cited support in A16.2 is not supported by the reviewed text of PMC11241823. See the append-only erratum in `analysis_plan_v1.md`.
