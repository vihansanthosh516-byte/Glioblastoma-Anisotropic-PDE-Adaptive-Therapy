# Research plan (form text), Track C, version 2, 2026-11-04

Drafted with AI assistance (Claude). The author must read, edit, and own every claim, and state AI use wherever the venue requires it. Each field below is at most 2,000 characters. Status: DRAFT. It replaces `paper/research_plan_trackC_draft.md` (2026-10-27), which was written before the model was checked against published trial numbers. Numbers from the runs are in `paper/trackc_stage2_report_2026-11-04.md` and `vault/PAPER_FINDINGS_LEDGER.md` (rows C1 to C33); they are interim.

## Rationale

Glioblastoma (GBM) almost always returns after standard treatment, which includes temozolomide (TMZ) in repeated courses (Stupp et al., 2005). Adaptive therapy times drug courses to tumor burden instead of giving the most drug early; it can delay progression in prostate cancer (Zhang et al., 2017). Computer studies of dosing rules are often tested only inside the model that produced them. They rarely check that the model matches what real trials found. This project asks a narrow question about a simulated GBM tumor with sensitive and resistant cells: after the simulation is forced to match published trial numbers, can a dosing rule that reads scans beat a fixed schedule that uses the same amount of drug? The anchors are the median progression-free survival of the standard arm (5.5 months), no clear gain from dose-dense TMZ, and no clear gain from 12 instead of 6 cycles (Gilbert et al., 2013; Balaña et al., 2020). The work is simulation only. It uses no patient data and gives no medical advice. A negative answer is useful: it shows where dosing-rule research needs new measurements before more modeling. An earlier version of this project reported a gain for a belief-based rule in an uncalibrated toy world; this plan tests that kind of claim against trial numbers. No novelty is claimed. The prior-art search was not exhaustive.

## Research Question/Hypothesis(es)/Engineering Goal(s)

Research question: In simulated GBM tumors that match published trial numbers, how much time to progression can any dosing rule gain over the best fixed schedule at the same mean number of courses, and does the answer depend on how progression is called?

Hypotheses (declared before the runs; all exploratory):
- H1: Worlds that match the standard-arm median and both negative controls need a weak drug effect.
- H2: In those worlds no rule that reads scans beats the best fixed schedule at equal drug use.
- H3: Only a controller that knows each patient's parameters has room above fixed dosing, and the size of that room depends on the progression rule (RANO-like or baseline-referenced).
- H4: Knowing one patient-level quantity (growth rate, kill strength, initial burden) captures a large share of that room.

Expected outcomes: H1 and H2 are likely. H3 may fail for one progression rule. H4 is open. A result where nothing beats fixed dosing is valid and will be reported unchanged.

Engineering goal: a tested simulation pipeline in which every number can be traced to a file: seed registry, run manifests (git commit, code hash, seeds), and read-back checks.

## Materials

Computer: desktop PC, 8 logical cores, 31 GB RAM, Windows 11, no GPU. No purchases are needed.

Software: Python 3 (version recorded in each run manifest), numpy, scipy, pandas, pytest, git. Code: `src/trackc_*.py` and scripts 190 to 204.

Data: none from patients. Published numbers are read from the trial papers listed in the bibliography and recorded with page and level of confidence in `paper/trackc_calibration_targets_2026-11-03.md`.

Records: a git repository, a design file written before the runs (`paper/trackc_anchored_spec_2026-11-04.md`), a seed registry (`output/trackc_seed_registry.json`), one manifest per run folder, and a findings ledger (`vault/PAPER_FINDINGS_LEDGER.md`).

There are no chemicals, no organisms, no biological agents, and no contact with human participants or animals.

## Procedures

1. Declare. Commit the design and the code before a run. Draw seeds only from the registry. A run writes its manifest first (git commit, code hash, seeds, clean-tree flag).
2. Model. Sensitive and resistant cells compete (logistic growth). Drug kill is m x growth rate. Cycles last 28 days with a 5-day course; scans every 84 days; at most 12 courses. Progression is called by a RANO-like rule (true burden at least 1.4 times the smallest earlier scan) or by a baseline rule (1.2 times baseline). Both are reported.
3. Calibrate (test C). Sample 500 worlds from assumed ranges. Keep worlds whose simulated standard arm (6 courses) has a Kaplan-Meier median of 143 to 186 days.
4. Negative controls (test D). On fresh patients, dose-dense dosing (drug effect times 2) must give a hazard ratio of 0.75 to 1.00, and 12 versus 6 courses (landmark at the day-168 scan) 0.71 to 1.38. Surviving worlds are the controlled set.
5. Variants. A responder/non-responder mixture, a mixed course cap (3 or 6), an MGMT-stratified check, response depth, integration step.
6. Regime map (test E) and frontier (script 202). In the controlled worlds compare fixed schedules, scan-reading rules, response-guided stopping and an informed controller at equal mean courses.
7. Value of information (test H) and mechanism (test M).
8. Verify. Read each result back from its file; run tests in a scratch copy; use independent read-only checkers. No sealed test without a new declaration.

## Risk Assessment

Human participants: none. No patient data are used and no person is contacted.

Chemical, biological, electrical, and fire hazards: none. No lab work.

Computer use: runs take hours on a desktop PC; take breaks and keep the PC ventilated.

Scientific and misuse risk: simulation results could be misread as medical advice or as a validated "digital twin". Every figure and the abstract will say "simulation, assumed parameters, exploratory, not clinical". The study makes no claim of clinical benefit, no patient-specific claim and no digital-twin claim. Known limits are stated beside each result: one model class, assumed parameter ranges, a weak drug effect that the anchors imply, no toxicity or early death, true burden without measurement error, and many post hoc choices (listed in the report). Negative results are reported unchanged. Methods for AI use: Claude drafted code and text; the author checks every claim and number against the files.

## Data Analysis

Unit: a simulated world (a sample of an assumed family) with 1,000 simulated patients; patients are paired across arms (same patients in every arm).

Primary quantities: restricted mean time to progression at 365 and 730 days; mean courses used; the gap of an arm above the frontier of fixed schedules at the arm's own mean number of courses. Calibration uses the Kaplan-Meier median with linear interpolation between scan days. Hazard ratios use a Cox model with Efron ties.

Reporting: medians and ranges over worlds; counts of worlds above and below the frontier. Sign tests are descriptive and uncorrected; with 11 worlds the smallest possible p-value is 0.001, so no significance claim is made. Sensitivity runs: course cap, dose-dense multiplier, landmark estimate, progression rule, scan interval, integration step.

Checks: simulated outputs are re-read from disk; numbers in the report are compared with the files by independent checkers; tables list the weakest point beside each result.

## Bibliography

Balaña, C., Vaz, M. A., Sepúlveda, J. M., et al. (2020). A phase II randomized, multicenter, open-label trial of continuing adjuvant temozolomide beyond 6 cycles in patients with glioblastoma (GEINO 14-01). *Neuro-Oncology, 22*(12), 1851–1861. https://doi.org/10.1093/neuonc/noaa107

Blumenthal, D. T., Gorlia, T., Gilbert, M. R., et al. (2017). Is more better? The impact of extended adjuvant temozolomide in newly diagnosed glioblastoma: a secondary analysis of EORTC and NRG Oncology/RTOG. *Neuro-Oncology, 19*(8), 1119–1126. https://doi.org/10.1093/neuonc/nox025

Ellingson, B. M., Wen, P. Y., & Cloughesy, T. F. (2017). Modified criteria for radiographic response assessment in glioblastoma clinical trials. *Neurotherapeutics, 14*, 307–320. https://doi.org/10.1007/s13311-016-0507-6

Gilbert, M. R., Wang, M., Aldape, K. D., et al. (2013). Dose-dense temozolomide for newly diagnosed glioblastoma: a randomized phase III clinical trial. *Journal of Clinical Oncology, 31*(32), 4085–4091. https://doi.org/10.1200/JCO.2013.49.6968

Stupp, R., Mason, W. P., van den Bent, M. J., et al. (2005). Radiotherapy plus concomitant and adjuvant temozolomide for glioblastoma. *New England Journal of Medicine, 352*, 987–996. https://doi.org/10.1056/NEJMoa043330

Zhang, J., Cunningham, J. J., Brown, J. S., & Gatenby, R. A. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. *Nature Communications, 8*, 1816. https://doi.org/10.1038/s41467-017-01968-5
