# How noisy is the line? Measurement error and the reliability of volumetric progression criteria (RANO 2.0) in glioblastoma

PLAN, 2026-10-18. Nothing here is run. Novelty is UNVERIFIED (about 15 searches, snippet level). "PhD-level" here means method standards (formal model, partial identification, pre-declaration, simulation checks, sensitivity analysis). It does not mean the results are publishable; that depends on what the data show.

## 1. Problem and why it matters
- RANO 2.0 calls progression at a volume increase of 40% or more and partial response at a decrease of 65% or more, with a 0.5 cm3 measurability floor (https://pmc.ncbi.nlm.nih.gov/articles/PMC11630866/). Volumes now come from automated segmentation.
- Segmentation is noisy (inter-rater ICC 0.22 to 0.96 across studies); automated tools need expert supervision (https://pmc.ncbi.nlm.nih.gov/articles/PMC10513769/). No automated tool was formally validated against RANO 2.0 as of Nov 2024 (https://www.frontiersin.org/journals/neurology/articles/10.3389/fneur.2026.1797487/full).
- Progression-free survival (PFS) in glioblastoma changes with how progression is defined: median PFS 8.8, 9.5 and 11.8 months for three definitions in a 73-patient cohort (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11241823/). In other cancers, mismeasured progression biases median PFS by months in simulations (https://www.frontiersin.org/journals/drug-safety-and-regulation/articles/10.3389/fdsfr.2024.1423493/full).
- Who is affected: patients (treatment is changed or stopped on a progression call), trials (PFS is an endpoint; bias changes conclusions and sample size), and regulators and tool makers (AI volumes enter these decisions).

## 2. Position against the literature
Closest work found: a 2026 npj paper turns voxel-wise segmentation uncertainty into 95% credible intervals on tumour volume to separate noise from change (https://www.nature.com/articles/s41746-026-02902-0). It is model-internal uncertainty and does not test the RANO 2.0 lines or trial endpoints (from the snippet; paper not read in full). Not found: measured (tool-to-tool and expert-correction) noise put against the 40% / 65% lines; a minimal detectable change curve versus tumour size; a trial-level bias simulation for glioblastoma volumetric criteria. These are application-level gaps, not a new method.

## 3. Aims and hypotheses
- A1 Measurement model and identifiability. Scan i of patient p, tool k: y = ln V + b_k + e, e ~ N(0, s_k(V)^2), s_k decreasing with volume. With two tools the data identify only s_1^2 + s_2^2 (sum of variances of the change difference), not each. Prove this and derive a partial-identification bracket by adding a third measurement (expert correction of an AI mask, UPenn-GBM).
- A2 Misclassification and minimal detectable change. For each size band: MDC95 = 1.96 x s_change; probability that a true change of +40% is observed as progression, and that a true +10% is observed as progression. Closed forms from the log-normal model, checked by simulation.
- A3 Trial consequences. Simulate two-arm trials with an assessment schedule, true growth distribution and the noise bracket; outcomes: bias in median PFS, bias in hazard ratio, loss of power; compare single-scan calls with confirmed calls (mRANO-style).
- A4 Translation. An open-source function that returns, for two volumes and a tool, the probability of each RANO 2.0 class; a short reporting recommendation; a patient timeline viewer with uncertainty bands.

## 4. Declared analyses (first two, Amendment 15 in analysis_plan_v1.md)
Test 1: how often does the RANO 2.0 class differ between the two automated tools on real consecutive scans? Test 2: the MDC curve by size band against the progression line, with the noise bracket. Success and reading rules are fixed in the amendment before any run.

## 5. Possible conclusions, written before the data (each outcome has a meaning)
| Result | What it would mean | Recommendation we could support |
|---|---|---|
| Flip rate high and MDC above the 40% line for small tumours | The line is inside the noise for small tumours; calls near it depend on the tool | Report uncertainty; raise the measurability floor or require confirmation for small lesions; do not mix tools in one trial |
| Flip rate low and MDC below the line at all sizes | RANO 2.0 volumetric calls are robust to these tools | Supports use of automated volumes; limits stated |
| Result depends on where in the noise bracket the truth lies | Data cannot separate tool error from definition difference | State the bracket; call for scan-rescan data (the missing evidence) |

## 6. Public benefit (honest route, not a promise)
Better-calibrated progression calls can reduce premature treatment changes and false reassurance, and can make trial endpoints less biased. What this project can deliver: quantified error bounds, an open calculator, a reporting recommendation, and a clear statement of the data needed. It cannot show clinical benefit; that needs prospective data and clinical review.

## 7. Limits (stated now)
- Two tools with different label definitions are not a scan-rescan pair; their disagreement overstates noise. Expert correction of an AI mask understates it. Result = a bracket.
- LUMIERE is one centre; ratings cover selected follow-ups; ratings use 2D RANO.
- The closest literature may already cover parts of A2.
- Trial simulations depend on assumed growth and schedule; they are sensitivity analyses, not forecasts of real trials.

## 8. Timeline to February (fair date and forms unverified)
Weeks 1-2: scripts and tests 1-2 (declared). Weeks 3-4: A1 proof and simulation checks. Weeks 5-7: A3 declared then run. Weeks 8-9: calculator and timeline viewer. Weeks 10-13: paper, poster, mock interviews. Weeks 14-15: buffer and forms. Ask a faculty member or clinician to review the framing (not yet done).
