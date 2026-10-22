# Does measurement noise matter for volumetric RANO progression calls? A technical note on Amendments 15-21

**DRAFT, updated 2026-10-21. Not peer reviewed. Not clinical evidence. Not for submission.**
Author of record: the project owner (personal research). Drafting assistant: Claude (Sonnet 5.5).

**How this draft was built.** Every project number below was read from an output file or a ledger entry in this session. The file or entry is cited next to it. Numbers I could check against a file are marked as checked in Appendix A. Anything I could not confirm is marked **[UNVERIFIED]**. Literature access is labelled in section 7.

---

## Abstract

**Question.** RANO 2.0 assesses volumetric progression against a context-specific reference: the pretreatment baseline or the smallest measurement after treatment starts, with special rules in some settings. This project's simulations use a simplified running-minimum nadir and a 40% increase threshold. How reliable is a call of this kind when an automated segmentation tool measures volume with error?

**What was done.** Four kinds of analysis were documented in amendments A15-A21. The commit history records an amendment, code and result sequence for the result-producing analyses, but does not establish when scripts were run. The real-data analysis is exploratory; simulation results depend on stated assumptions:
1. A comparison of two automated tools on real follow-up scans (LUMIERE, 121 eligible consecutive-scan pairs from 43 patients).
2. Trial simulations with *assumed* independent noise.
3. A simulation whose noise level was taken from *one* published 2-day repeat-scan study (Abu Khalaf 2021).
4. An audit of how that study's repeatability numbers can be converted to a noise level.

**What was found.**
- The two tools gave different consecutive-scan categories (PD-like, SD-like or PR-like ratio categories, not RANO calls) in 54.5% of eligible pairs (95% CI 45.6-63.8%). This is tool disagreement. It is not a noise estimate.
- In simulation, if per-scan noise is independent and Gaussian, single-scan nadir-based calls move median PFS earlier and cut power. The effect becomes at least one month from a change-SD of about 0.10 (log units).
- Three readings of the one published repeatability study, with their confidence limits, convert to a change-SD anywhere from 0.079 to 0.342 (the full converted range). **Not every value in that range was simulated.** A19 simulated four values (0.143, 0.199, 0.235, 0.342) and A21 a selected seven (about 0.08 to 0.262), 11 distinct values in all. The direction of the effect (earlier PFS, lost power) is the same at every simulated value; for single-scan calls it was material from 0.111 upward and not material at 0.08. The modeled shift ranges from 5.8 to 15.0 weeks across the lowest and highest point-estimate readings, a factor of about 2.6 across endpoints that combine different readings and tool versions. Within-version comparisons of R1 versus R3 give factors of about 2.0 (version 1.2) and 2.4 (version 2.0).
- A two-scan confirmation rule reduces the bias but does not remove it in most cells. At N = 100 and 8-week scans, at the simulated values up to 0.143 it keeps the median-PFS bias below the declared one-month materiality line (4.345 weeks); the power loss under that rule was still material from 0.126 upward.

**What was not measured.** The project has no scan-rescan noise estimate from its own pipeline. Every effect size rests on an assumed or converted noise value. Real noise may be correlated across scans, depend on volume, and include systematic shifts. These are only partly modelled.

**What this is not.** It is not a clinical result, a validated tool, a new theorem, or a claim of novelty. The calculation of false-progression risk against the number of scans is routine and was not pursued as a result.

---

## 1. Research question

RANO 2.0 sets a volumetric progression threshold of a 40% increase and a response threshold of a 65% decrease in tumour burden. The primary article states that these were proposed to be consistent with the 25% and 50% changes in area of the two-dimensional rule, that two-dimensional measurement remains the primary measure, and that studies to date have not shown a conclusive benefit of volumetric over two-dimensional analysis [Wen 2023, S9: article text fetched and summarized; its tables and data supplement were not readable]. The Ellingson guide says the 3D thresholds are "adopted from mRANO" and that the criteria do not clearly describe how to extend the 2D measurable-lesion rule to 3D [S1, summarizing fetch]; the derivation in the Wen data supplement was not read.

**The reference scan depends on context.**
- Progression is compared with the smallest tumour measurement at either the pretreatment baseline or after therapy began (the nadir). Response is compared with the setting-specific baseline described below [S9; S1 gives the same reference-scan distinction in its figure legend].
- The baseline itself depends on the setting: the post-radiotherapy MRI (about 21-35 days after radiotherapy) in newly diagnosed glioma treated with radiotherapy; the post-surgery, pretreatment MRI if there is no radiotherapy; the pretreatment MRI for recurrent glioma [S9].
- After a confirmed pseudoprogression the guide measures progression against that time point or the nadir after it, and special rules apply when the nadir scan has no measurable disease [S1].

**This project's simulations use a simplified reference:** progression is the first scan at least 1.40 times the running minimum of earlier *measured* volumes, starting from a baseline scan taken as the true nadir. They do not model pseudoprogression, a change of baseline, a response phase, or the baseline-timing rules. The real-data analysis A15 used consecutive-scan ratios (and A18 Part D tool-specific nadirs), not the full RANO 2.0 reference.

Automated segmentation tools produce volumes with some error. The project asked:

- **Q1.** Do two automated tools agree on the consecutive-scan volumetric category (PD-like, SD-like, PR-like) of the same scan pairs? (real data)
- **Q2.** If measurement noise has a stated size, what does it do to median PFS and trial power under a nadir-based progression rule? (simulation)
- **Q3.** What noise size does the one published scan-rescan study imply, and how sensitive are Q2's conclusions to how that number is converted? (simulation plus an assumption audit)

Not asked: how large the noise actually is in this project's data. The data to answer that were not available (section 4).

**Scope notes taken from the sources.**
- The 0.5 cm3 (500 mm3) floor used in A15-A21 is **the project's operational choice** [analysis_plan_v1.md A18.1 and the errata block]. The RANO 2.0 articles state a measurable-lesion rule in *diameters* (3D: at least 10 mm in all three perpendicular planes [S9, S1]). **This note does not claim the floor is equivalent to that rule, and the floor was not derived from it.** For scale only: a sphere 10 mm across has a volume of about 524 mm3 (pi/6 x 10^3). That is an illustration, not an equivalence; a lesion that meets the diameter rule can have a different shape and volume.
- A15's categories are *consecutive-scan* ratio categories, not RANO 2.0 calls. The "PR-like" category (ratio at or below 0.35, the -65% line) is **not** a RANO partial response: RANO response is judged against the setting-specific baseline and must be sustained for at least 4 weeks [S9]; the Ellingson guide explicitly requires response-confirmation scans after at least 4 weeks to demonstrate durable response [S1, Table 2 note]. A15 compared consecutive scans [analysis_plan_v1.md A18.1; ledger 2br and 2bt errata]. "PD-like" is likewise not a RANO progression call: RANO's reference for progression depends on context (the pretreatment baseline, or the smallest measurement after treatment begins, with special rules after a confirmed pseudoprogression [S9, S1]), whereas A15's PD-like category uses the ratio of consecutive scans.
- The simulation's "two-scan confirmation" is the project's simplification. For newly diagnosed glioma after radiotherapy, RANO 2.0 requires confirmation of progression within 12 weeks after radiotherapy completion. The Ellingson guide also uses "3 months after RT" in its Introduction, the Figure 3 legend and Table 2 footnote g, without reconciling that phrasing with "12 weeks" in its text, Table 1 and the Figure 4 legend (as seen in a summarizing fetch, 2026-10-21; Figure 3 uses "≥3 months" and footnote g ">3 months"). Taken literally, Table 1 makes confirmation mandatory at 12 weeks or less and the waiver starts after 3 months, which leaves an unaddressed window of about one week; the source does not say how to treat it. The guide describes other contexts where progression confirmation is optional, and requires response-confirmation scans at least 4 weeks later. The simulation models neither response nor these RANO confirmation rules [S1; analysis_plan_v1.md A21.5].

---

## 2. Methods and analysis history

All amendments are in `analysis_plan_v1.md`. Commit and code-commit identifiers are as recorded in the ledger. **Commit order was checked read-only against `git log` on 2026-10-21 (section 2.4).**

| Amendment | Ledger | Script | Output files | Type |
|---|---|---|---|---|
| A15 tool dependence and MDC | 2br | `src/140_rano_tool_dependence.py` | `output/rano_noise/rano_tool_dependence.json`, `rano_pairs.csv` | Exploratory, real data |
| A16 trial simulation | 2bs | `src/141_trial_sim.py` | `output/trial_sim/trial_sim_summary.json`, `trial_sim_cells.csv`, `trial_sim_E0.csv` | Simulation, assumed parameters |
| A18 sensitivity | 2bt | `src/142_rano_sensitivity.py` | `output/rano_sensitivity/rano_sensitivity_summary.json`, `rano_sensitivity_cells.csv`, `rano_sensitivity_pairs.csv` | Simulation (A-C) plus post hoc real data (D) |
| A19 anchored simulation | 2bu | `src/143_anchored_trial_sim.py` | `output/anchored_sim/anchored_sim_summary.json`, `anchored_sim_cells.csv` | Simulation anchored to one study |
| A21 conversion audit | 2bx | `src/145_conversion_audit.py` | `output/conversion_audit/conversion_audit_summary.json`, `conversion_audit_cells.csv`, `conversion_table.csv` | Simulation, assumption audit |
| A17 calculator | 2bv | `src/rano_calc.py`, `tests/test_rano_calc.py` | none (tool) | Arithmetic tool, no result |
| A20 timeline viewer | 2bw | `src/144_rano_timeline_viewer.py` | `output/viewer_rano/viewer_rano.html`, `viewer_rano_info.json` | Display, no result |

(Amendment 17 is listed after 21 only to keep the result-producing analyses together. A17 and A20 produced tools, not findings. The script and file names exist on disk; I confirmed that by listing them, not by running them.)

### 2.1 Data (real-data analyses)

LUMIERE longitudinal glioblastoma MRI with expert RANO ratings. Two automated segmentations were used, HD-GLIO-AUTO and DeepBraTumIA, taken as supplied in the dataset's feature tables, contrast-enhancing label [analysis_plan_v1.md A15.2]. **The project did not run either tool.** The dataset paper reports 91 patients, 638 study dates and 2,487 images [S3; Europe PMC record and abstract counts, 2026-10-21]. These are dataset-wide counts, and this project's 91 patients and 430 pairs come from the dataset tables. They are **not** the same as the 80 patients and 502 post-operative time points that Suter et al. 2023 analysed after excluding patients [S8]; section 6 describes that paper.

### 2.2 Real-data definitions (A15)

- Pair: two consecutive post-operative scans of one patient. Eligible if the enhancing volume is at least 500 mm3 at both scans in both tools.
- Class from volume ratio r = V1/V0 within a tool: "PD-like" if r >= 1.40, "PR-like" (PR-threshold) if r <= 0.35, otherwise "SD-like". These are consecutive-scan ratio categories, not RANO calls (section 1). In the A15 tables below, PD, SD and PR are shorthand for these -like categories.
- Flip: the two tools give different categories for the same pair. Patient-cluster bootstrap 95% CI (10,000 draws, seed 20261018 as declared).
- Reading rule declared in advance: "material" if the CI lower bound is at least 0.20. [A15.3]

### 2.3 Simulation model (A16, reused in A18, A19, A21)

- Units are weeks. Two arms 1:1, N per arm (50 or 100). True volume grows exponentially from baseline; baseline is the true nadir. Control median time to a true +40% rise is 34.8 weeks (8.0 months), **an assumption**. Analysis plan A16.2 says it lies in a published range (8.8 to 11.8 months in a 73-patient cohort) and cites PMC11241823. **That citation does not support the value.** I opened PMC11241823 in the 2026-10-21 audit: it is Belue et al., Diagnostics 2024;14(13):1374, 92 patients with glioblastoma. The text I read gives no median PFS and no 73-patient cohort (its survival plots are images; the last ~1,240 characters of the page were not read). Also, 8.0 months is below the 8.8-month lower end that the plan itself states. Treat 34.8 weeks as a plain assumption with no verified source. Log growth rate has an assumed spread of 0.6.
- Scans every D weeks (8 or 12) to week 104. Administrative censoring at the last scan.
- Measured log volume = true log volume + e. e is normal, independent across scans, non-differential by arm. The error SD of one scan is s/sqrt(2), where s is the SD of the *change* error between two scans.
- Single-scan rule: progression at the first scan whose volume is at least 1.40 times the running minimum of earlier measured volumes (a simplified nadir reference; see section 1 for what RANO 2.0 actually uses in each setting). Two-scan rule: the criterion must hold at that scan and at the next scan against the same reference.
- Response, the 500 mm3 floor and death are not modelled in A16 (death is added in A18 Part C). Growth is monotone.
- Estimands, all relative to the same rule with s = 0: E1 shift in control-arm Kaplan-Meier median PFS; E2 change in power at hazard-ratio setting 0.7; E4 change in Peto log hazard ratio; E5 share of events called earlier than the true crossing minus one scan interval; type I error at no treatment effect.
- Materiality, declared in advance: E1 at least 4.345 weeks in absolute value with a Monte Carlo CI excluding 0; power loss at least 0.05 with CI excluding 0.
- 2,000 simulated trials per cell. Seed 20261018. Common random numbers across cells.

### 2.4 Analysis history and limits on preregistration

The amendments record test declarations, and the ledger says each analysis was run once. **This audit did not independently verify run counts or whether scripts were run before their code was committed.** It checked commit order only, from `git log` author dates (2026-10-18, local time):

| Amendment | Amendment text committed | Code committed | Result committed |
|---|---|---|---|
| A15 | 51fd053, 16:53:33 | ad9c55c, 16:55:13 | 945e2ab, 16:57:32 |
| A16 | 1b7f4ec, 17:00:35 | 40c2db6, 17:05:44 | 93d030c, 17:06:30 |
| A17 (calculator) | c1a1d6a, 17:26:20 | 28a75ab, 22:28:41 | tests/ledger 2bv: no separate result commit checked |
| A18 | 17c95dc, 17:33:55 | 0e8435c, 17:36:16 | efa5551, 17:38:27 |
| A19 | d9b924b, 22:24:05 | df5cacd, 22:25:36 | 29bdb56, 22:26:27 |
| A20 (viewer) | 390fac2, 22:31:30 | 29b70b3, 22:32:49 | 6e860bb, 22:35:25 |
| A21 | 2d3e6a4, 22:39:59 | 48c5933, 22:41:01 | 131495f, 22:42:15 |

For A15, A16, A18, A19, A20 and A21, the recorded amendment commit precedes the code commit, which precedes the result commit. A17 has no separate result commit checked. **This shows commit order only.** It cannot show that scripts were not run locally before their code was committed, that each was run once, or that test choices were made without sight of earlier results. A16 was declared about 3 minutes after the A15 result commit, and A15 records that aggregate noise summaries (script 104) were seen before its declaration. The plan file was also edited after the runs by errata commits (4700f31 and 2147d75 on 2026-10-19, 48c8850 on 2026-10-20). The 2026-10-19/20 audit errata (end of `analysis_plan_v1.md`) withdraw several labels and are applied throughout this note.

---

## 3. Results

Each block gives: what it shows, what it supports, the weakest limitation, and the evidence.

### 3.1 Real-data tool disagreement (A15, A18 Part D)

This section reports **disagreement between two tools**. It does not estimate measurement noise. The two tools use different label definitions and are not a scan-rescan pair [analysis_plan_v1.md A15.1, errata].

#### 3.1.1 Class flip rate (A15 test 1)

| Quantity | Value | Source |
|---|---|---|
| Consecutive post-operative pairs | 430 (91 patients in tables) | `rano_tool_dependence.json`: `n_consecutive_postop_pairs`, `n_patients_in_lumiere_tables` |
| Excluded: tool volume missing, zero, or below 500 mm3 | 309 (71.9%) | `n_pairs_excluded`, `share_excluded_by_floor_or_missing` |
| Eligible pairs / patients | 121 / 43 | `n_pairs_eligible`, `n_patients_eligible` |
| Flip rate | 0.545 (CI 0.456-0.638) | `TEST1_flip_rate` |
| One tool PD-like, other not | 0.471 (CI 0.381-0.552) | `TEST1a_one_tool_PD_other_not` |
| Declared reading | "material tool dependence" | `TEST1_reading` |

Cross-tabulation, HD-GLIO-AUTO (rows) by DeepBraTumIA (columns), from `class_crosstab_hd_by_dbt`:

| HD-GLIO-AUTO \ DeepBraTumIA | PD-like | SD-like | PR-like |
|---|---|---|---|
| PD-like | 30 | 17 | 9 |
| SD-like | 25 | 25 | 5 |
| PR-like | 6 | 4 | 0 |

(The file labels these PD, SD and PR; the "-like" is added here because they are consecutive-scan ratio categories, not RANO calls.)

Derived by me from that table (not stored in the file): the tools agree on 55 of 121 pairs (30 + 25 + 0). **The tools never agree on the PR-like category: 14 pairs are PR-like in HD-GLIO-AUTO and 10 in DeepBraTumIA (marginal counts in `class_counts`), and the PR-like/PR-like cell is 0, so none of those 24 pair-calls coincide.** This says nothing about RANO partial response.

- **Supports:** on large post-operative lesions, the consecutive-scan class depends heavily on which tool produced the volume.
- **Weakest limitation:** the flip rate is not an error rate. It mixes label-definition differences, tool error and real change. Seventy-two percent of pairs were dropped, mostly because HD-GLIO-AUTO had no volume for a scan (159 and 134 missing cells; see 3.1.4), which may mean "no enhancing tumor found" rather than noise. Eligible pairs come from later scans and a higher share of expert progression. One centre.
- **Evidence:** `output/rano_noise/rano_tool_dependence.json`; ledger 2br.
- **CI difference, code-reading explanation not confirmed by rerun:** the point estimate is 0.545 in both outputs. A15's interval is 0.4559-0.6381; A18 Part D's 500-floor interval is 0.4556-0.6354. Both use the patient-cluster bootstrap with the same seed, number of draws, and percentile method. A15 computes this interval first in its seeded stream; A18 computes the 500-floor interval third in a separately initialized stream with the same seed, after the 0- and 250-floor bootstraps. Different call order would produce different resamples. This is consistent with the interval differences, but other causes were not excluded and no bootstrap was rerun to confirm the explanation.

#### 3.1.2 Agreement with expert rating (A15 test 1b)

| | HD-GLIO-AUTO | DeepBraTumIA |
|---|---|---|
| Sensitivity for expert PD (83 expert-PD pairs of 106) | 0.554 (0.463-0.681) | 0.506 (0.387-0.652) |
| Specificity (23 non-PD pairs, derived as 106 - 83) | 0.609 (0.444-0.800) | 0.739 (0.500-0.917) |

- **Supports:** neither tool's PD-like category reproduces the expert progression rating closely on these pairs.
- **Weakest limitation:** expert ratings are 2D RANO at selected follow-ups, judged against baseline or nadir, not against the previous scan. Specificity rests on 23 pairs. The comparison is not like for like.
- **Evidence:** `rano_tool_dependence.json`: `TEST1b_vs_expert_rating`; ledger 2br.

#### 3.1.3 Minimal detectable change versus the +40% line (A15 test 2)

Line: ln(1.40) = 0.3365 (`line_ln_1.40` = 0.33647). For each size band, two **proxy** change-SDs were computed. The JSON keys call them `s_change_upper_bound` and `s_change_lower_bound`. **Those names are withdrawn** by the 2026-10-20 errata; both are proxy scenarios, not bounds. Below I call them *cross-tool proxy* and *expert-correction proxy*.

| Band (median volume) | Pairs / patients | Cross-tool proxy s | MDC95 (ln) | Expert-correction proxy s | MDC95 (ln) | Reading against 0.3365 |
|---|---|---|---|---|---|---|
| 0 (1.6 cm3) | 40 / 22 | 1.286 | 2.52 | 1.881 | 3.69 | both above the line |
| 1 (3.7 cm3) | 40 / 19 | 1.133 | 2.22 | 0.013 | 0.026 | proxies on opposite sides |
| 2 (11.4 cm3) | 41 / 19 | 0.718 | 1.41 | 0.000 | 0.000 | proxies on opposite sides |

(Band 0's "lower" proxy, 1.881, is larger than its "upper" one, 1.286. That alone shows the old names were wrong.)

- **Supports:** for the smallest lesions in this set, the +40% line is far below the disagreement between the tools. For larger lesions the two proxies disagree about which side of the line the noise is on.
- **Weakest limitation:** the cross-tool proxy assumes independent tool errors and can contain definition differences; the expert-correction proxy anchors on an AI mask and is degenerate (near 0) in larger bands. The expert-correction proxy uses 10, 24 and 180 scans in the three bands. Derived by me: the cross-tool values (0.72-1.29) are roughly 2 to 16 times the published-anchor conversions (the full converted range 0.079-0.342, section 3.3; not every value in it was simulated), so they cannot be read as noise.
- **Evidence:** `rano_tool_dependence.json`: `TEST2_mdc_vs_progression_line`; ledger 2br and its 2026-10-20 erratum.
- The JSON also stores `"reading": "inside the noise"` for band 0 and `"undetermined"` for bands 1 and 2. The "inside the noise" wording is withdrawn by the errata.

#### 3.1.4 Sensitivity of the flip rate to the floor and to the nadir reference (A18 Part D, post hoc)

| Floor (mm3) | Pairs / patients | Flip rate | 95% CI |
|---|---|---|---|
| 0 (any positive volume) | 214 / 66 | 0.626 | 0.562-0.693 |
| 250 | 148 / 50 | 0.581 | 0.497-0.665 |
| 500 | 121 / 43 | 0.545 | 0.456-0.635 |
| 1000 | 87 / 33 | 0.529 | 0.426-0.645 |
| Each tool against its own nadir (floor 500) | 68 / 20 | 0.544 | 0.377-0.663 |

Exclusion flow of the 430 pairs: 216 with a missing tool volume, 93 with a positive volume below 500 mm3 (and none missing), 121 eligible. Missing cells by column: HD-GLIO-AUTO 159 (first scan) and 134 (second scan); DeepBraTumIA 16 and 17. Eligible pairs have a later median time since the first post-operative scan (56 vs 45 weeks) and more expert PD ratings (78.3% of 106 rated vs 57.6% of 262 rated).

- **Supports:** the flip rate stays between 53% and 63% across floors, and 54% against each tool's own nadir. The declared prediction PD1 (lower CI >= 0.20 everywhere) is supported.
- **Weakest limitation:** the nadir analysis keeps only 68 pairs from 20 patients. Missing HD-GLIO-AUTO volumes drive most exclusions. Post hoc. Still a tool-disagreement figure.
- **Evidence:** `output/rano_sensitivity/rano_sensitivity_summary.json`: `partD`, `PD1_supported`; ledger 2bt.

### 3.2 Simulations with assumed noise (A16, A18 Parts A-C)

**All results in this section follow from the model in section 2.3. They show what noise of a stated size could do. They are not evidence about any real trial.**

#### 3.2.1 A16: effect of noise at N = 100 per arm, D = 8 weeks

True control median to a +40% rise is 34.8 weeks; the no-noise simulated control median is 38.856 weeks (discretisation adds about 4 weeks).

| s (change-SD) | Rule | Control median PFS (wk) | E1 shift (wk) | Power at HR setting 0.7 | Type I error | E5 early-call share |
|---|---|---|---|---|---|---|
| 0 | single | 38.856 | 0 | 0.969 | 0.052 | 0 |
| 0.25 | single | 23.544 | -15.312 | 0.288 | 0.045 | 0.549 |
| 0.5 | single | 16.072 | -22.784 | 0.087 | 0.044 | 0.678 |
| 0.75 | single | 15.952 | -22.904 | 0.061 | 0.045 | 0.716 |
| 1.0 | single | 15.732 | -23.124 | 0.054 | 0.044 | 0.734 |
| 0 | two-scan | 38.856 | 0 | 0.967 | 0.053 | 0 |
| 0.25 | two-scan | 30.632 | -8.224 | 0.583 | 0.047 | 0.388 |
| 0.5 | two-scan | 23.972 | -14.884 | 0.177 | 0.054 | 0.529 |
| 0.75 | two-scan | 22.256 | -16.600 | 0.098 | 0.049 | 0.576 |
| 1.0 | two-scan | 20.620 | -18.236 | 0.073 | 0.049 | 0.599 |

(The volume-dependent scenario "vol" is listed in the same file; for single scan E1 = -23.08, power 0.051.) E4 (Peto log HR attenuation, single scan): +0.40 (s = 0.25), +0.55 (0.5), +0.59 (0.75), +0.61 (1.0).

- **Supports:** if per-scan noise is independent and Gaussian, the nadir rule gives earlier median PFS and less power; the declared predictions H-S1 and H-S2 are supported. Material E1 occurs in 39 of 40 (N, D, s > 0, rule) combinations; material power loss in 40 of 40. The one non-material E1 cell is N = 100, D = 12, two-scan, s = 0.25 (E1 = -3.69, below 4.345).
- **Weakest limitation:** the size of every effect is set by the assumed noise model (independent, normal, constant or size-banded SD) and by assumed monotone growth with an assumed spread. Noise at s = 0.5 or more is already as large as the +40% line (0.3365 in ln), so those cells show saturation rather than a realistic range. The A16 single-scan cell at s = 0.25 lies inside the converted range (0.079-0.342) and gives a -15.312-week shift; this is one modeled scenario, not a measured effect. There is no response phase and no stable disease.
- **Negative and null results, kept as declared:**
  - **H-S3 not supported.** Its first clause held (|E1| 16.6 weeks with two-scan vs 22.9 with single-scan at s = 0.75) but its second clause failed: with no noise the two-scan and single-scan control medians are both 38.856 weeks, because monotone growth means confirmation never changes the call. The ledger calls this a flaw in the prediction.
  - **H-S4 not supported as declared.** Four of 48 type I error cells lie outside the declared band [0.04, 0.06] (0.039, 0.0395, 0.039, 0.0605), each by less than the Monte Carlo SE of about 0.005. No sign of real inflation, but the declared rule failed.
  - **Discretisation check E0:** +4.32 to +4.36 weeks at D = 8 (expected about 4) and +5.75 (N = 50) and +4.97 (N = 100) at D = 12 (expected about 6). The N = 100, D = 12 value is about one week under the expected order.
- **Evidence:** `output/trial_sim/trial_sim_summary.json` (`H_S1`..`H_S4`, `E0_discretisation`, `E1_material_cells`, `E2_material_loss_cells`); `trial_sim_cells.csv` (table above, N = 100, D = 8 rows); ledger 2bs.

#### 3.2.2 A18 Part A: finer noise grid (N = 100, D = 8)

| s | E1 single (wk) | Power change single | E1 two-scan (wk) | Power change two-scan |
|---|---|---|---|---|
| 0.05 | -0.50 | -0.0085 | -0.23 | -0.006 |
| 0.10 | -4.69 | -0.0535 | -1.30 | -0.0265 |
| 0.15 | -8.29 | -0.2345 | -3.98 | -0.099 |
| 0.20 | -14.05 | -0.4945 | -6.41 | -0.222 |
| 0.25 | -15.31 | -0.6815 | -8.22 | -0.384 |

- **Supports:** the bias does not need large noise. Single-scan E1 first becomes material at s = 0.10 (about 10% change-SD). Two-scan E1 first becomes material at 0.20 and two-scan power loss at 0.15. PA1 and PA2 supported. The s = 0.25 cells reproduce A16 exactly (the file's reproduction check).
- **Weakest limitation:** same model as 3.2.1.
- **Evidence:** `rano_sensitivity_cells.csv` (rows `A`); `rano_sensitivity_summary.json` (`partA_first_material`, `reproduction_check_vs_A16`); ledger 2bt.

#### 3.2.3 A18 Part B: correlated noise (AR(1) between successive scans; marginal per-scan SD sigma)

E1 single-scan (weeks):

| sigma | rho 0 | rho 0.5 | rho 0.8 |
|---|---|---|---|
| 0.10 | -7.50 | -5.48 | -2.35 |
| 0.25 | -20.94 | -15.18 | -9.44 |
| 0.50 | -22.89 | -22.66 | -17.23 |

- **Supports:** serial correlation reduces the bias but does not remove it, at least in this AR(1) form. PB1 supported.
- **Weakest limitation:** AR(1) is a convenience. Real correlation structure (same reader, same scanner, tool version) is unknown. Whether independent noise is the worst case in general was **not** shown, only in these cells.
- **Evidence:** `rano_sensitivity_cells.csv` (rows `B`); `rano_sensitivity_summary.json` (`partB_single_scan`); ledger 2bt.

#### 3.2.4 A18 Part C: death as a competing event (median death 78 weeks, assumed)

E1 single-scan: -11.82 weeks at s = 0.25 and -16.05 at s = 0.5; two-scan -7.20 and -11.03. Control median with no noise 32.0 weeks. Same direction, smaller size than A16. No prediction was declared.

- **Evidence:** `rano_sensitivity_cells.csv` (rows `C`); ledger 2bt. **Weakest limitation:** death time is independent exponential and assumed.

### 3.3 The published repeatability anchor and its conversion (A19, A21)

#### 3.3.1 The anchor (Abu Khalaf 2021)

- Twenty patients with recurrent glioblastoma, one 1.5 T scanner, repeat scans 45.1 hours apart on average (SD 7.6, range 23.2-54.6), automated BraTumIA 1.2 and 2.0 [reference S2: full text fetched, summarized].
- Enhancing tumour repeatability coefficient in percent-change terms (RC%) 46% (CI 33-67) for BraTumIA 1.2 and 39% (28-57) for 2.0; RC in mL 6.9 (4.9-10) and 5.2 (3.7-7.5) [`analysis_plan_v1.md` A19.1, A21.2; the 46% [33-67] row was also seen in the fetched article table].
- Enhancing volume was larger on the later scan (mean 24.2 vs 26.3 mL for version 1.2 and 22.9 vs 24.9 mL for 2.0, per A19.1; ratio of means 26.3/24.2 = 1.087 and 24.9/22.9 = 1.087, i.e. 8.7% by my arithmetic). The article's Discussion separately reports median volume increases of 10.0% (version 1.2) and 10.3% (version 2.0) and median specific growth rates of 2.1% and 2.3% per day, as seen in a summarizing fetch on 2026-10-21. The mean-based 8.7% and the median-based 10% are different statistics, not a conflict. The authors suggest incomplete gadolinium washout; they say tumour growth cannot be fully excluded but seems unlikely to explain it.
- The patients were postoperative, in a bevacizumab and irinotecan study, with large enhancing volumes (about 25 mL). 19 of the 20 have public imaging in the RIDER Neuro MRI collection on TCIA.
- **The paper's text, as fetched, does not say whether RC% was computed on raw, percent or log data. The within-subject SD appears only in the Online Supplemental Data. I could not open the supplement (HTTP 403 from the journal; the article text gives no supplement URL; Europe PMC lists no supplementary files and marks the article as not open access, checked 2026-10-21).** The article text also describes no adjustment for the later-scan rise before RC was computed, so whether it was removed is unknown. The authors warn that the enhancing-tumour repeatability figures "should be interpreted with caution". The scale of the RC% remains unresolved, and no scale is inferred here.

#### 3.3.2 Conversion used in A19 and the three readings in A21

The conversion from RC% to the change-SD s is an **assumption** (the paper gives no formula) [analysis_plan_v1.md A19.1, A21.1]:

| Reading | Rule | BraTumIA 1.2: point (CI) | BraTumIA 2.0: point (CI) |
|---|---|---|---|
| R1 linear (used in A19) | s = RC% / 1.96 | 0.235 (0.168-0.342) | 0.199 (0.143-0.291) |
| R2 log back-transform | s = ln(1 + RC%/100) / 1.96 | 0.193 (0.145-0.262) | 0.168 (0.126-0.230) |
| R3 absolute | s = (RC in mL / mean volume in mL) / 1.96 | 0.139 (0.099-0.202) | 0.111 (0.079-0.160) |

**Full converted range (all readings, both versions, point estimates and confidence limits): 0.079 to 0.342.** Source: `conversion_audit_summary.json` (`conversion_table`, `readings`, `union_range`); `conversion_table.csv`.

**What was actually simulated is smaller than that range.** A19 simulated four values (0.143, 0.199, 0.235, 0.342). A21 simulated a selected seven (0.08, 0.111, 0.126, 0.139, 0.168, 0.193, 0.262), which span about 0.08 to 0.262. Together that is 11 distinct values (ledger 2bx erratum). The file key `union_range` names the full converted range, not the simulated set. Per the A21 erratum the ledger counts 17 distinct converted values (the table has 18 entries; by my reading 0.1680 and 0.1684 coincide when rounded to three decimals, which would explain 17, but I did not confirm how the ledger counted); six converted values the ledger lists as not simulated are 0.099, 0.145, 0.16, 0.202, 0.23 and 0.291, and the lower end 0.079 was represented only by the nearby 0.08. The top of the full range (0.342) was simulated in A19 but not in A21.

#### 3.3.3 A19: simulation at the four converted values (N = 100, D = 8)

| s | Single-scan E1 (wk) | Single-scan power change | Two-scan E1 (wk) | Two-scan power change |
|---|---|---|---|---|
| 0.143 | -7.64 | -0.20 | -3.55 (not material) | -0.077 |
| 0.199 | -14.00 | -0.48 | -6.38 | -0.216 |
| 0.235 | -14.99 | -0.63 | -7.68 | -0.327 |
| 0.342 | -20.44 | -0.82 | -12.28 | -0.603 |

At D = 12 and N = 100, single-scan E1 is -3.70 (not material), -6.65, -10.59 and -15.40. Across all 32 cells (N 50/100, D 8/12, both rules) E1 is material in 23 and power loss in 32. Type I error stays within [0.04, 0.06] in all eight N = 100, D = 8 cells.

- **Supports:** with this one study's converted noise, in this model, single-scan nadir calls come about 8 to 20 weeks earlier (7.6 to 20.4 at N = 100, D = 8) and power falls a lot. A longer scan interval shrinks the bias but power loss stays material.
- **Weakest limitation:** one study (20 patients, one scanner, one protocol, large recurrent tumours, a different tool from the one used for the A15 data). The interval is two days, so the figure may include real change and a systematic later-scan rise that the model does not include. The conversion is assumed.
- **Evidence:** `output/anchored_sim/anchored_sim_summary.json` (`PX1`..`PX4_supported`, `N100_D8`, `material_E1_cells`, `material_E2_cells`) and `anchored_sim_cells.csv`; ledger 2bu.

#### 3.3.4 A21: simulation at a selected seven of the converted values (about 0.08 to 0.262; N = 100, D = 8)

| s | Single-scan E1 (wk) | Single-scan power change | Two-scan E1 (wk) | Two-scan power change |
|---|---|---|---|---|
| 0.080 | -2.57 (not material) | -0.0285 (not material) | -0.60 | -0.0185 |
| 0.111 | -5.77 | -0.0715 | -1.78 | -0.037 |
| 0.126 | -6.71 | -0.1185 | -2.53 | -0.053 |
| 0.139 | -7.32 | -0.1795 | -3.32 | -0.0715 |
| 0.168 | -10.46 | -0.326 | -5.06 | -0.143 |
| 0.193 | -13.54 | -0.459 | -6.20 | -0.2005 |
| 0.262 | -15.72 | -0.703 | -8.74 | -0.409 |

- **Supports:** the direction (earlier PFS, lost power) holds for every simulated value from s = 0.111 up. At the point estimates of the readings the single-scan shift is 5.8 to 15.0 weeks and the power change 0.07 to 0.63 (the upper ends come from A19's R1 values). At the lowest reading (R3 for version 2.0, s = 0.111) the two-scan rule keeps the bias to -1.8 weeks.
- **What does not survive:** the size. For the same version, the single-scan shift differs by 7.7 weeks (version 1.2) and 8.2 weeks (version 2.0) between readings R1 and R3 (declared prediction PU4, supported). These correspond to factors of about 2.0 (14.988/7.32) and 2.4 (13.996/5.768), respectively.
- **Weakest limitation:** the seven simulated values are a selected grid of the converted values (about 0.08 to 0.262), not all of them and not the whole 0.079-0.342 range (A21 erratum). Results at the other converted values are not available. Not identifiable from the published numbers: the scale of the RC, whether the later-scan rise was removed first, the error shape, and any volume dependence. Normal independent error only.
- **Evidence:** `output/conversion_audit/conversion_audit_summary.json` (`N100_D8_single_scan`, `N100_D8_two_scan`, `PU1`..`PU4_supported`, `PU4_detail`); `conversion_audit_cells.csv` (not opened by me); ledger 2bx and its erratum.

### 3.4 Tools (no scientific result)

`src/rano_calc.py` turns a stated change-SD and a normal-error model into class probabilities; `output/viewer_rano/viewer_rano.html` displays volumes with uncertainty bands. Seven declared calculator tests are author-reported as passed (ledger 2bv: the transcript is not retained; **I did not rerun them**). The viewer's integrity checks are recorded in ledger 2bw. They check arithmetic and display. Neither is validated on patients or reviewed by a clinician [ledger 2bv, 2bw]. They are mentioned only so that nobody mistakes them for results.

---

## 4. What was not measured

1. **Scan-rescan noise from the project's own pipeline.** None. The project never ran a segmentation tool on a repeat-scan pair. LUMIERE volumes are as supplied.
2. **The scale of the one published repeatability number.** Unresolved; the supplement was not obtained.
3. **Serial correlation of noise across visits** beyond an assumed AR(1).
4. **Volume dependence of noise.** Only the cross-tool proxy (not a noise measure) and one banded scenario "vol" were used.
5. **Systematic shift between scans** (the later-scan rise seen by Abu Khalaf). Not modelled.
6. **Response, pseudoprogression, steroid effects, new lesions, clinical status** and the actual RANO 2.0 confirmation algorithm.
7. **Any real trial.** No patient-level survival data were used.
8. **Separation of biology from noise.** These project analyses cannot separate them; a repeat scan 2 days apart is the closest published case, but it can still include biological or systematic change.

A feasibility check of the public QIN-GBM paired scans was done and **stopped**. Access rules, compute limits and an unvalidated segmentation path blocked it (`research_notes/final_feasibility_check_qin_path_2026-10-20_claude.md`). That memo is a feasibility record, not a result.

---

## 5. Limitations

- **Noise is assumed or converted everywhere in the effect-size results.** The simulations tell us what noise of a stated size would do. They cannot say what noise is.
- **Independent normal noise** is the base model. Correlation was tried only as AR(1); heavy tails and asymmetry were not tested.
- **Monotone exponential growth** with an assumed spread and an assumed 8-month true median. No response phase, so the nadir-selection effect in treated tumours is not represented.
- **Two-scan confirmation is a simplification,** not the RANO 2.0 rule.
- **One centre, two supplied tools** for the real-data part. Seventy-two percent of consecutive pairs were excluded. The floor of 500 mm3 is the project's operational choice.
- **One published repeatability study** (20 patients, large recurrent tumours, a 2-day gap, a different tool).
- **Exploratory, many cells.** Materiality was judged by rules recorded in the amendments, not adjusted for multiplicity. Monte Carlo CIs reflect simulation error only.
- **The author's earlier descriptions had errors** that the audit errata correct (floor labelled as a RANO rule, "upper/lower bound" labels, "union" of 17 values). This draft follows the corrected text.

---

## 6. Discussion

**What the evidence supports, in one sentence.** If per-scan volume noise is independent and Gaussian with a change-SD of about 0.10 or more, then in a model of exponentially growing tumours a single-scan nadir-referenced progression call shifts median PFS earlier by about a month or more and costs substantial power; one published repeatability study, read three ways, converts to a change-SD between 0.079 and 0.342 depending on reading and confidence limit. The 11 distinct simulated values span about 0.08 to 0.342; they do not cover every converted value. They show the same direction. The modeled shift is 5.8 to 15.0 weeks across the lowest and highest point-estimate readings, a factor of about 2.6 across endpoints that combine different readings and tool versions; within-version R1-to-R3 factors are about 2.0 and 2.4.

**What it does not support.** It does not show that RANO 2.0 volume thresholds are within the noise of any specific tool. It does not show that real trials are biased. The tool-disagreement figures cannot be used as noise levels; they are 2 to 16 times larger than the published-anchor values and mix definitions with error.

**The basic false-progression calculation is routine.** The chance that noise alone produces a false call rises with the number of scans, because each scan is another opportunity and the running minimum drifts down. A hand argument written for the project's decision notes bounds the chance between half the range-statistic tail and the full tail. It is standard order-statistics reasoning, was not checked by code, and is not a theorem or a new result. The false-progression-by-visit figure was therefore **not** produced; it would restate what A16 and A18 already imply.

**Relation to prior work, kept modest.**
- Hong et al. 2012 report in their abstract that measurement variability can attenuate the PFS treatment effect and may inflate type II error (abstract only; not read in full) [S4].
- Beaumont et al. 2015 derive thresholds of -35% and +55% for lung CT from an assumed 15% SD for a two-measurement ratio, and state they did not examine nadir-reference bias (full text fetched, summarized) [S5].
- Yoon et al. 2019 simulate the effect of reader variability on RECIST response and progression rates on CT (full text fetched, summarized) [S6].
- Suter et al. 2023 evaluate automated longitudinal measurements against expert RANO assessments in 80 LUMIERE patients with 502 post-operative time points. They used expert-identified reference time points because severe under-segmentation could wrongly reassign the reference time point, and they note that failing to detect a large lesion would wrongly update the nadir. Their best trend agreement with the expert was 81.1% (HD-GLIO volumes), which they call insufficient for use without expert input, and they recommend the tools as a time-saving aid with manual correction [S8; summarizing fetch of the web page, 2026-10-21; relevant sections only]. That paper analyses time to progression against overall survival; it does not analyse PFS or model the effect of noise on PFS or power. It supports the relevance of segmentation errors. It does not establish this project's noise-to-PFS result, and the 80 patients and 502 time points are its analysed subset, not the whole LUMIERE dataset (91 patients, 638 study dates).
- The RANO 2.0 neuroradiology guide gives no quantitative noise figure and says volumetric approaches "would arguably result in a better inter-reader agreement" (fetched, summarized) [S1].
- **Novelty is not claimed.** Searches were limited and no systematic review was done. The honest wording is: not found in our searches, with the above papers as the nearest items.

**What would change the picture.** A measured, same-pipeline scan-rescan noise value with a defensible scale; evidence on serial correlation; and a model that includes response and the actual confirmation rules. None is available in this project. The QIN-GBM route that might have supplied the first of these was assessed and stopped.

**Suggested handling of this note.** Present it as an assumption-audit and sensitivity analysis. Lead with the limitation (no measured noise). Have a clinician or imaging physicist review the framing before anyone cites it.

---

## 7. References (access status as actually obtained)

Status key: **FT-summary** = full text page fetched this session, read through a summarizing tool (not raw text); **abstract** or **record** = abstract or bibliographic index record only (for example Europe PMC); **not verified** = not checked. No reference is labelled "full text read" unless every section was read directly; none meets that bar here.

- **S1.** Ellingson BM, Sanvito F, Cloughesy TF, et al. A Neuroradiologist's Guide to Operationalizing the Response Assessment in Neuro-Oncology (RANO) Criteria Version 2.0 for Gliomas in Adults. AJNR Am J Neuroradiol 2024;45(12):1846-1856. doi:10.3174/ajnr.A8396. PMC11630866. **FT-summary** (article text via a summarizing tool, twice; tables and figures as quoted by that tool, not inspected as images); bibliographic details match the Europe PMC record (2026-10-21). Note: the earlier frontier note in `research_notes/` attributes PMC11630866 to "Wen et al."; the fetched page names Ellingson et al. The primary RANO 2.0 article is S9.
- **S9.** Wen PY, van den Bent M, Youssef G, Cloughesy TF, Ellingson BM, Weller M, et al. RANO 2.0: Update to the Response Assessment in Neuro-Oncology Criteria for High- and Low-Grade Gliomas in Adults. J Clin Oncol 2023;41(33):5187-5199. doi:10.1200/JCO.23.01059. PMC10860967. **FT-summary** (article text read through a summarizing tool; the final ~4,400 characters of the page were checked on 2026-10-21 and contain only site navigation, not article text; the tables, which are images, and the Data Supplement were **not read**). Authors, volume and pages match the Europe PMC record (2026-10-21).
- **S2.** Abu Khalaf N, Desjardins A, Vredenburgh JJ, Barboriak DP. Repeatability of Automated Image Segmentation with BraTumIA in Patients with Recurrent Glioblastoma. AJNR Am J Neuroradiol 2021;42(6):1080-1086. doi:10.3174/ajnr.A7071. PMC8191681. **FT-summary** (main text only, re-fetched 2026-10-21; bibliographic details match Europe PMC). **Online Supplemental Data not read** (journal pages returned HTTP 403; Europe PMC lists no supplementary files for this article).
- **S3.** Suter Y, Knecht U, Valenzuela W, Notter M, Hewer E, Schucht P, Wiest R, Reyes M. The LUMIERE dataset: Longitudinal Glioblastoma MRI with expert RANO evaluation. Sci Data 2022;9:768. doi:10.1038/s41597-022-01881-7. PMC9755255. **Record and abstract level only:** authors, volume, article number and the counts 91 / 638 / 2,487 come from the Europe PMC record (2026-10-21). The publisher page redirected to a login and the full text was not read in this audit.
- **S4.** Hong S, Schmitt N, Stone A, Denne J. Attenuation of treatment effect due to measurement variability in assessment of progression-free survival. Pharm Stat 2012;11(5):394-402. doi:10.1002/pst.1524. **Abstract** (Europe PMC, summarized). Full text not obtained (Wiley HTTP 403).
- **S5.** Beaumont H, Souchet S, Labatte JM, Iannessi A, Tolcher AW. Changes of lung tumour volume on CT - prediction of the reliability of assessments. Cancer Imaging 2015;15:17. doi:10.1186/s40644-015-0052-2. PMC4628325. **FT-summary** (bibliographic details match the Europe PMC record, 2026-10-21).
- **S6.** Yoon JH, Yoon SH, Hahn S. Development of an algorithm for evaluating the impact of measurement variability on response categorization in oncology trials. BMC Med Res Methodol 2019;19:90. doi:10.1186/s12874-019-0727-7. PMC6498480. **FT-summary** (bibliographic details match the Europe PMC record, 2026-10-21).
- **S7.** Sridhara R, Mandrekar SJ, Dodd LE. Missing data and measurement variability in assessing progression-free survival endpoint in randomized clinical trials. Clin Cancer Res 2013;19(10):2613-2620. doi:10.1158/1078-0432.CCR-12-2938. **Abstract** (Europe PMC record; title, authors, volume, issue and pages 2613-2620 match that record, 2026-10-21). The abstract and article were not read for content in this audit.
- **S8.** Suter Y, Notter M, Meier R, Loosli T, Schucht P, Wiest R, Reyes M, Knecht U. Evaluating automated longitudinal tumor measurements for glioblastoma response assessment. Front Radiol 2023;3:1211859. doi:10.3389/fradi.2023.1211859. **FT-summary** (web page read through a summarizing tool on 2026-10-21; authors, patient count 80 and time-point count 502 as returned by that tool; relevant sections only, not a full reread).
- **Tool references** (HD-GLIO-AUTO, DeepBraTumIA/BraTumIA): HD-GLIO and DeepBraTumIA project pages were read as READMEs this session, but the primary tool papers were **not read** and are not cited here.
- **Not cited because unresolved:** a "Walovitch 2013" paper (not found in four searches); the Schmitt journal-club slides (file retrieved, text not readable with available tools; whether they cover noisy nadirs is **unconfirmed**).
- **Methods named but not referenced:** Bland-Altman repeatability coefficient, Peto log hazard ratio, patient-cluster bootstrap. Standard methods; no source checked.

---

## 8. Unresolved facts that need your review

1. **Abu Khalaf RC% scale and within-subject SD.** The supplement is unread. The conversion (R1/R2/R3) stays an assumption until it is resolved. Optional: open the supplement in a normal browser, or send the unsent author question from `research_notes/final_feasibility_check_qin_path_2026-10-20_claude.md`.
2. **The later-scan rise.** Partly resolved: the article reports mean volumes 24.2 vs 26.3 mL and 22.9 vs 24.9 mL (8.7% by ratio of means, my arithmetic) and, in its Discussion, median increases of 10.0% and 10.3% (different statistics). **Still unknown:** whether the rise was removed or modelled before RC% was computed. The article text describes no adjustment, and the supplement was not accessible.
3. **Commit order and "run once".** Partly resolved: `git log` shows amendment, then code, then result commits for A15, A16, A18, A19, A20 and A21; A17 has no separate result commit checked (section 2.4). **Still unverified:** that no script ran before its code was committed, that each was run once, and that test choices were made without sight of earlier results (A15 records that script-104 noise summaries were seen before its declaration). The plan file was edited by errata commits after the runs.
4. **The calculator's 7/7 test pass** is author-reported with no transcript. I did not rerun it.
5. **Source for the 34.8-week assumed true median.** Checked and **not supported**: PMC11241823 (Belue et al. 2024) gives no median PFS in the text I read and its cohort has 92 patients, not the 73 the plan states; 8.0 months also lies below the plan's own 8.8-month lower end. The value stays an assumption. Not checked: whether the plan meant a different paper, and the final ~1,240 characters and the plot images of that page.
6. **Scope of source review.** S3 bibliographic details and dataset counts (91 / 638 / 2,487) are from the Europe PMC record; the dataset paper's full text was not read. S8 (Suter 2023: authors, 80 patients, 502 post-operative time points, 91 reviewed) was read through a summarizing fetch of relevant sections only. Neither article was comprehensively reread for this note.
7. **Walovitch 2013** and **Schmitt slides' coverage of noisy nadirs:** unresolved; not used.
8. **RANO source distinction:** PMC11630866 is Ellingson et al. 2024 guide; PMC10860967 is Wen et al. 2023 primary RANO 2.0 article. This note distinguishes them. The primary article's tables and Data Supplement (including how the 40% and 65% volumetric values were derived, its Table S2) were not read.
9. **Whether independent noise is the worst case** for serial correlation: shown only in the A18 Part B cells, not in general.
10. **Files I did not open:** `rano_pairs.csv`, `rano_sensitivity_pairs.csv`, `conversion_audit_cells.csv`, `trial_sim_E0.csv`, `viewer_rano_info.json`, `paper/phd_level_plan_probabilistic_rano.md`, `paper/overview_for_review_rano.md`. Statements about them come from the ledger or the summary files.
11. **A15/A18 confidence intervals:** A15's interval is the first bootstrap call in a stream initialized with the declared seed. A18's 500-floor interval is the third bootstrap call in its separately initialized stream with the same seed, after the 0- and 250-floor intervals. This would produce different resamples; other causes were not excluded, and no rerun confirmed the explanation. The bounds are 0.45588-0.63810 in A15 and 0.45555-0.63542 in A18.

---

## Appendix A. Number check log

**Read from files this session and matched to the ledger:**
- `rano_tool_dependence.json`: pairs 430, excluded 309 (0.7186), eligible 121, patients 43, flip 0.5455 (0.4559-0.6381), one-tool-PD 0.4711 (0.3814-0.552), class counts, crosstab, TEST1b values, MDC band values.
- `trial_sim_summary.json` and `trial_sim_cells.csv`: all values in 3.2.1, H-S1..H-S4, E0, 39 and 40 material cells.
- `rano_sensitivity_summary.json` and `rano_sensitivity_cells.csv`: Parts A, B, C, D values, reproduction check.
- `anchored_sim_summary.json` and `anchored_sim_cells.csv`: the A19 values, D = 12 values, 23 and 32 material cells.
- `conversion_audit_summary.json` and `conversion_table.csv`: the conversion table, the N = 100, D = 8 results for seven values, PU4 detail.

**Read directly from output files in the 2026-10-21 audit (no code run):** `trial_sim_summary.json` (E1 material cells 39, E2 material loss cells 40, H-S4 outside list, E0 values) and every N = 50/100, D = 8/12 row of `trial_sim_cells.csv`. The 40 cells are 2 N x 2 D x 2 rules x 5 noise settings (s = 0.25, 0.5, 0.75, 1.0, vol). The one cell with |E1| below 4.345 weeks is N = 100, D = 12, two-scan, s = 0.25 (E1 = -3.69); the next smallest is -4.5 (N = 50, D = 12, two-scan, s = 0.25). The smallest signed power change is -0.319; its absolute loss, 0.319, exceeds the 0.05 materiality threshold. The A19 table was checked against `anchored_sim_summary.json` and the A21 tables against `conversion_audit_summary.json`.

**Taken from the ledger or amendment text only (not recomputed):** commit identifiers; the tests of the calculator; the A19.1 published RC numbers beyond the single table row seen in the fetched article; A16's assumed parameters; the reading rules.

**Derived by me (arithmetic on file values, not stored in a file):** 55 of 121 agreements; PR-like/PR-like cell = 0; 23 non-PD pairs (106 - 83); "2 to 16 times" ratio of cross-tool proxy to the converted range (0.718/0.342 and 1.286/0.079); the 8.7% later-scan rise from the A19.1 volumes (26.3/24.2); the cross-version endpoint factor 14.988/5.768 = 2.60; within-version factors 14.988/7.32 = 2.05 (version 1.2) and 13.996/5.768 = 2.43 (version 2.0).
