# Track B: prior-art check for the VIG hybrid (searched 2026-11-14 and 2026-11-15)

Scope: the two ideas inside VIG (`src/track_b_vig_lumiere.py`).
1. The time to the next scan carries information about what the tumor does (short gap = doctors expect change).
2. The last observed change partly reverses at the next scan, more after a fall than after a rise.

This file is a literature record. It contains no new result from our data.

## 1. How the search was done

- Tool: a paper index over PubMed, PubMed Central, bioRxiv, medRxiv and arXiv. About 16 queries, 12 to 40 results each.
- Query topics: next-scan glioma volume forecasting; scan interval as a predictor; informative visit or observation times; short-interval follow-up MRI in glioma; rebound and regression to the mean in serial tumor volumes; LUMIERE and MU-Glioma-Post forecasting; irregular visit times in health-record models.
- Depth: for four papers, passages of the full text were read (marked F). All others were read as abstract only (marked A).
- Not searched: Scopus, Web of Science, Google Scholar, conference-only proceedings (MICCAI workshops, SPIE), theses, patents, non-English papers. A citation-graph search on the closest paper returned nothing.

## 2. Idea 1: the scan interval

### 2a. The observation is already published (this is the most important finding)

| Paper | Depth | What it shows | Difference from VIG |
|---|---|---|---|
| Lindner & Fiehler 2026, BMC Med Imaging, PMC13421742, doi 10.1186/s12880-026-02612-9 | F | 60 patients from UCSF-ALPTDG matched to UCSF-PDGM (52 complete cases). Logistic model for "enhancing tumor increased": scan interval odds ratio 0.91 per day (95% CI 0.84 to 0.99), p = 0.023. Text: "The scan interval was inversely associated with progression risk". | The interval is an adjustment covariate. In the passages read, the authors do not explain it, do not call it informative scheduling, and do not use it to forecast volume. In their linear model for net volume change the interval was not significant (p = 0.454). Same cohort as one of ours. |

Result: **the sign of the interval effect (shorter gap, more growth) is not new.** It was reported on a subset of one of our own cohorts.

### 2b. The mechanism is known clinical practice

| Paper | Depth | What it shows |
|---|---|---|
| van Dijken et al. 2023, PMC10605145 | F | Unscheduled scans in glioma follow-up led to major treatment consequences in about 30%. |
| GIN CUP survey, PMC8043861 | F | "Additional short-interval imaging was carried out in cases of possible pseudoprogression in most centres (71%)". |
| BT-RADS papers, for example Kim 2022, pmid 36045502 | A | A structured report score sets the follow-up interval from the reader's suspicion. |
| Ji et al. 2021, PMC8099465; Kang et al. 2021, PMC8284622 | A | Recommended MRI intervals for glioma derived from progression hazard. |

### 2c. The statistical idea family is known

| Paper | Depth | What it shows |
|---|---|---|
| Buzkova 2010, pmid 21969983 | A | States that clinicians ask patients to return sooner based on the last exam, so observation times are informative (example includes tumor metastases). |
| Gasparini 2020, PMC6919310; Thomadakis 2025, PMC11911807; Stamoulopoulos 2026, PMC13454517; Coley et al., arXiv 1508.07511; Sun et al., arXiv 2004.12234; Lesperance 2015, pmid 26179520 | A | Joint models and weighting methods for informative visit processes. Their aim is unbiased estimation, not a better next-visit forecast. |
| Lipton 2016, arXiv 1606.04130; TA-RNN, PMC11211851; Multi-FIT, arXiv 1905.00125 | A | "Time since last observation" or the missingness pattern used as a model input in health-record prediction. |

### 2d. Glioma forecasting methods that use the interval in another way

| Paper | Depth | Use of the interval | Difference from VIG |
|---|---|---|---|
| OASA, Jung & Shin, arXiv 2609.12435 | F | Elapsed days / 180 as a constant input channel (time scaling). MU-Glioma-Post, 15 test triplets. Anchors the forecast to the last observation. Dice equal to persistence (0.6071 vs 0.6070). | Voxel maps and Dice, not volume error. Interval is elapsed time, not a scheduling signal. Shares the "stay close to persistence" idea. |
| STRIDE, arXiv 2610.08848 | A | Inter-scan interval conditions a latent transition. | Classifies stable / pseudo-progression / progression. Not a volume forecast. |
| TRACE, arXiv 2606.30313 | A | Scan interval as a pass-through concept. LUMIERE. | Response classification. |
| Weizman 2012, pmid 23286047; Glazar 2020, PMC7409184; Kebaili 2025, arXiv 2509.10824; Liu 2023, arXiv 2309.05406; PREDICT-GBM, PMC13558737; Pash 2025, PMC12132268 | A | Time is the integration horizon of a growth or generative model. | None treats the gap as information about the doctor's expectation. |

### 2e. What was NOT found for idea 1

No paper was found that does all of these: (i) forecasts next-scan glioma volume, (ii) uses the gap to the next scan as a weight between a growth law and a constant drift, and (iii) tests against persistence, population rate, two-point, random slope, Gompertz and gradient boosting.

## 3. Idea 2: the asymmetric reversal term

- Related biology: rebound growth after bevacizumab is stopped (Zuniga 2010, pmid 20151176, A); depth of response and time to regrowth after anti-VEGF therapy (Ellingson 2023, PMC10592195, A); regrowth after radiotherapy in low-grade glioma (Plaszczynski 2023, PMC10128962, A). These are treatment-specific growth models.
- Related statistics: regression to the mean and mean-reverting state-space forecasts are standard (for example arXiv 2002.11228, A; arXiv 2509.04718, A).
- Not found: a next-scan tumor-volume forecast with separate reversal weights for a prior fall and a prior rise.
- Honest reading: the term is regression to the mean with two slopes. It is a small step on a standard idea.

## 4. Verdict on novelty

| Part of VIG | Status |
|---|---|
| "Shorter scan gap goes with more enhancing growth" | **Not new.** Lindner & Fiehler 2026, same UCSF data. |
| "Visit times depend on the doctor's suspicion" | **Not new.** Known in practice and in biostatistics. |
| Interval as a model input | **Not new** as a general idea (health-record models, STRIDE, TRACE, OASA). |
| Interval as an informative-visit weight inside a growth law for next-scan volume | Not found in this search. Narrow. |
| Two-slope reversal term for volume forecasts | Not found in this search. Small step on regression to the mean. |
| Benchmark of all six standard baselines over four post-treatment cohorts, leave-one-cohort-out | Not found in this search. |

Plain statement: VIG is a new combination of known parts. It is not a new principle. A reviewer who knows Lindner & Fiehler 2026 can say the interval effect was already reported. Any text about VIG must cite that paper.

## 5. Weakest points of this check

- One index, about 16 queries. Novelty is "not found", not "proven absent".
- Most papers were read as abstract only. A methods section can contain an interval term that the abstract does not mention.
- For Lindner & Fiehler the full text was read in passages, not every page.
- No conference-only papers, theses or patents were searched.
- The check does not change the result status: the four-cohort result in `paper/trackB_code_audit_2026-11-14.md` stays exploratory and post hoc.

## 6. A design limit this search made clear

The gap to the next scan is set by the doctor after reading the current scan and the patient's state. A model that uses it borrows the doctor's judgement. This is allowed by the locked question (scan dates are inputs). It must be stated in every report, and VIG must not be described as predicting biology.
