# Track C calibration targets from published sources (2026-11-03)

Status: Stage 1 working record of `paper/trackc_restart_plan_2026-11-03.md`. No simulation was run for this file. Claude, with five read-only Haiku 5.5 reference checkers (report 1: trial outcomes; 2: MGMT, growth, plasticity, real data; 3: gap-closing for trial numbers and the label; 4: response criteria and scan interval; 5: trial time origins and scan schedules). The reports themselves are not stored in the repository; this file is the record.

## How to read this table

- **P** = a checker opened a primary page (journal abstract record, PMC full text, or the drug label) and reported the number. The checkers read pages through a text summary, so numbers are *not* checked character by character against the PDF. Check against the PDF before quoting in a paper.
- **S** = secondary source only (review, news, repackaged summary). Do not use as a constraint.
- **X** = search snippet only, or not found. Do not use.
- The checkers could not open PubMed, NEJM, ASCO or Elsevier pages (cookie walls or 403). They used Europe PMC records, PMC pages and DailyMed.
- Errors in my own requests, found by the checkers: I cited RTOG 0525 as "NEJM 370:699". That page is the RTOG 0825 bevacizumab trial. RTOG 0525 is Gilbert MR et al., J Clin Oncol 2013;31(32):4085-4091. I also gave the mRANO paper the wrong journal.

## A. Standard-of-care schedule (what one "course" is)

| Item | Value | Level | Source |
|---|---|---|---|
| Maintenance temozolomide (TMZ) | Days 1-5 of each 28-day cycle; cycle 1 at 150 mg/m2/day; cycles 2-6 may rise to 200 mg/m2/day | P (label, via summary) | DailyMed Temodar label, section 2.2: https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=046a9011-3911-4d3f-a15f-fbb56d5aad56 |
| Number of maintenance cycles | 6; the label gives no separate maximum | P (label, via summary) | same |
| Concomitant phase | 75 mg/m2 once daily for 42-49 days with focal radiotherapy | P (label, via summary) | same |
| Planned adjuvant cycles in the Stupp trial | 6 (150-200 mg/m2, 5 of 28 days) | P | Stupp 2005 abstract, Europe PMC record of PMID 15758009 |

Model meaning: one course = 5 treated days; one cycle = 28 days; the standard arm has 6 courses. The old model used 12 courses.

## B. Outcomes of the standard arm (what the model must be able to reproduce)

| Item | Value | Level | Source |
|---|---|---|---|
| Stupp 2005 (EORTC 26981 / NCIC CE.3), n = 573 | Median overall survival 14.6 vs 12.1 months (radiotherapy + TMZ vs radiotherapy); HR for death 0.63 (95% CI 0.52-0.75); 2-year survival 26.5% vs 10.4% | P | N Engl J Med 2005;352:987-996, PMID 15758009, DOI 10.1056/NEJMoa043330 (Europe PMC abstract; OS also on the DailyMed label, section 14.1) |
| Stupp 2005 median progression-free survival | 6.9 vs 5.0 months in one summary page; 7.2 vs 5.0 months in a news snippet; 6-month PFS not found | S / X | OncoLink summary (S); NCI snippet (X). **Not primary-verified. Do not use as a constraint.** |
| Stupp 2009 5-year analysis | 2-year survival 27.2% (22.2-32.5) vs 10.9% (7.6-14.8); 5-year 9.8% (6.4-14.0) vs 1.9% (0.6-4.4); HR 0.6 (0.5-0.7) | P | Lancet Oncol 2009;10:459-466, PMID 19269895 (Europe PMC abstract) |
| RTOG 0525 (standard-dose arm), n = 833 randomised | Median PFS **5.5 months (95% CI 4.7-6.1)**; median OS 16.6 months (14.9-18.0). Randomisation took place **after** the concomitant radiotherapy + TMZ phase and before adjuvant TMZ; PFS and OS run from random assignment. Progression: at least 25% increase in the product of the two largest diameters, or a new lesion (Macdonald-type rule). Median adjuvant cycles given: **3** (standard arm), 4 (dose-dense arm); at least 6 cycles: 37% (155 patients) standard, 43% (181) dose-dense. | P (checker report 5, Methods and Results read through a page summary) | Gilbert MR et al., J Clin Oncol 2013;31(32):4085-4091, PMID 24101040; https://pmc.ncbi.nlm.nih.gov/articles/PMC3816958 |
| GEINO 14-01, stop-at-6 arm, n = 79 | Median PFS 7.7 months (5.7-9.8); 6-month PFS 55.7% (45.8-67.8); median OS 23.3 months (17.9-28.7). **Not used as a standard-arm anchor:** patients entered only after 6 adjuvant cycles without progression, and PFS runs "from the time of inclusion", so the clock starts about 6 months into adjuvant treatment and is conditional on being progression-free then (checker report 5). | P | Neuro-Oncology 2020;22(12):1851-1861, DOI 10.1093/neuonc/noaa107; https://pmc.ncbi.nlm.nih.gov/articles/PMC7746946 |

**Anchored value for the model's standard arm (proposed):** median progression-free time of **5.5 months (95% CI 4.7-6.1)**, the RTOG 0525 standard arm, which is 143 to 186 model days (30.44 days per month). RTOG 0525's clock starts at random assignment after chemoradiation, near the first adjuvant course, where the model's day 0 sits. Weakest points: one trial, one arm; its progression rule and scan schedule differ from the model's (section F); patients who progressed during chemoradiation are not in RTOG 0525 (833 randomised of 1,173 registered), so it describes a selected group; only 37% of the standard arm received 6 or more cycles (median 3), so "6 courses" is the planned schedule, not what most patients received.

## C. Negative controls (the model should show little benefit from these changes)

| Comparison | Published result | Level | Source |
|---|---|---|---|
| Dose-dense vs standard TMZ (RTOG 0525) | OS: HR 1.03 (0.88-1.20), P = .63, median 14.9 (13.7-16.5) vs 16.6 (14.9-18.0) months. PFS: HR 0.87 (0.75-1.00), P = .06, median 6.7 (6.2-7.7) vs 5.5 (4.7-6.1) months. | P | PMC3816958 |
| 12 vs 6 cycles, randomised (GEINO 14-01) | PFS: HR 0.99 (0.71-1.38), P = .95. OS: HR 1.30 (0.90-1.88), P = .16. Phase II, not powered for medians. Randomisation came after cycle 6 (entry required 6 cycles without progression), so this is a **landmark** comparison: the model version must start at the end of cycle 6 and include only scenarios still progression-free then. Progression by RANO; MRI every 3 months. | P | PMC7746946 |
| More than 6 vs 6 cycles, post hoc (Blumenthal et al.; secondary analysis of EORTC and NRG/RTOG trials) | Landmark analysis from the start of cycle 6 + 28 days, only for patients progression-free then. PFS: adjusted HR 0.80 (0.65-0.98), P = .03. OS: HR 0.92 (0.71-1.19), P = .52. 333 patients at 6 cycles, 291 at more than 6. The authors call it post hoc and non-randomised. | P | Neuro-Oncology 2017;19(8):1119-1126, DOI 10.1093/neuonc/nox025; PMC5570239. Weakest point (my inference): who received extra cycles was a decision, not a randomisation, so confounding by indication is likely. The PFS benefit disappears for OS. |
| 12 vs 6 cycles, small randomised trial in high-grade glioma | Median PFS 18 (14.8-21.1) vs 16 (11.0-20.9) months, log-rank P = .88. | P, but not GBM-only and the PFS is much longer than in the GBM trials | Front Oncol 2024;14:1357789; PMC11106464. **Not used as a GBM target.** |

## D. Growth of untreated tumours

| Item | Value | Level | Source |
|---|---|---|---|
| Stensjoen 2015, n = 106 | Median specific growth rate 1.4% per day (range -1.9 to 13.2% per day); median equivalent volume doubling time 49.6 days | P | Neuro-Oncology 2015;17(10):1402-1411; PMC4578579 |
| Roder 2024, n = 70 | Median volume doubling time 31 days (range -42 to 3013; the authors flag negative values as artefacts) | P | Neuro-Oncol Adv 2024;6(1):vdae053; PMC11046984 |

The model draws r (growth rate) log-uniformly from 0.001 to 0.05 per day. Net tumour growth at burden B is r(1 - B) with B0 between 0.05 and 0.5. Whether the model's typical regrowth speed is slower than these measured pre-treatment rates is a Stage 2 test C item (not yet computed).

## E. MGMT and published aggregates (side branch; not used by the simulation)

| Item | Value | Level | Source |
|---|---|---|---|
| Hegi 2005, 206 assessable tumours, 45% methylated | Methylated, radiotherapy + TMZ: median OS 21.7 months (17.4-30.4); methylated, radiotherapy alone: 15.3 months (13.0-20.9), P = 0.007; HR for methylation 0.45 (0.32-0.61) | P (abstract) | N Engl J Med 2005;352:997-1003, PMID 15758010. Unmethylated medians and 2-year survival were **not found** in any opened page. |
| RTOG 0525, both arms pooled | Methylated vs unmethylated: OS 21.2 (17.9-24.8) vs 14.0 (12.9-14.7) months; PFS 8.7 (6.6-11.2) vs 5.7 (5.1-6.1) months. By arm, standard dose: methylated OS 21.4 (17.6-29.0), unmethylated 14.6 (13.2-16.5). | P | PMC3816958. The checker saw the OS hazard ratio printed as 1.74 in one place and 0.58 in another (reciprocal directions of one comparison). |
| CeTeG / NOA-09 (MGMT-methylated only), n = 141 randomised | Median OS 31.4 (27.7-47.1) months with TMZ vs 48.1 (32.6-not reached) with lomustine + TMZ; HR 0.60 (0.35-1.03). The authors call it a small trial. | P (abstract) | Herrlinger U et al., Lancet 2019;393:678-688 |
| Stupp 2009 by MGMT status | Not found in any opened page | X | |

## F. Response criteria and scan schedule

| Item | Value | Level | Source |
|---|---|---|---|
| RANO 2010 (Wen et al.) full text | **Not open.** Only the abstract record was reached. | X | Europe PMC record of PMID 20231676 |
| mRANO (Ellingson et al.), PMC5398984, DOI 10.1007/s13311-016-0507-6 | Progression: sum of products of perpendicular diameters rises by 25% or more relative to the **smallest** measurement (not relative to baseline). A volumetric threshold of 40% is given as the equivalent. Measurable lesion: both diameters at least 10 mm. Confirmation: two sequential scans at least 4 weeks apart. | P, but summary-level only; mRANO, not RANO 2010 | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5398984/fullTextXML |
| RANO 2.0 (Wen 2023) | Record found, full text not fetched (HTTP 500) | X | PMC10860967 |
| MRI schedule inside the trials | RTOG 0525: MRI about 4 weeks after radiotherapy and before cycle 4 (also before cycles 7 and 10 if given), then every 3 months until progression. GEINO 14-01: MRI every 3 months until progression; progression by RANO. | P (checker report 5) | PMC3816958; PMC7746946 |
| MRI interval, EANO 2021 guideline | "An initial interval between scans of 2-6 months is common practice" after completion of therapy; a scan 3-4 weeks after radiotherapy gives a new baseline. No statement for the adjuvant-TMZ phase was found (the last 28,000 characters of the text were not read). | P (summary-level) | Weller M et al., Nat Rev Clin Oncol 2021;18:170-186; PMC7904519 |
| Volume change per TMZ cycle | **Not found** in any open paper | X | |

Consequence: the model's progression rule (true burden reaches 1.2 times the **baseline**) is not the RANO rule (about 1.4 times the **smallest earlier volume**). A tumour that first shrinks to 10% of baseline would need to regrow to 14% of baseline to count as RANO progression, but to 120% of baseline in the model. The model's time to progression is therefore longer than RANO progression-free time for the same tumour. The baseline rule follows the adaptive-therapy convention (progression when burden exceeds the pre-treatment level). The spec reports both.

## G. Real patient data and plasticity

| Item | Value | Level | Source |
|---|---|---|---|
| Ho et al. 2020 | "Leukopenia is a biomarker for effective temozolomide dosing and predicts overall survival of patients with glioblastoma", Mol Clin Oncol 2020;13(6):80, DOI 10.3892/mco.2020.2150. n = 141, retrospective, single centre, 2012-2017. Per patient: blood counts at baseline and up to six adjuvant cycles, leukopenia grade, completion of treatment (88 of 141 completed). Median OS 19 months; leukopenia 21 vs 14 months (P = 0.008). Licence CC BY-NC-ND 4.0. Per-patient values may sit in supplementary Table SI on the article page; the checker also saw "on reasonable request". **Not downloaded.** | P | https://spandidos-publications.com/10.3892/mco.2020.2150 |
| Liau 2017 (Cell Stem Cell 20:233) and Rabe 2020 (Cell Death Dis 11:19) | Reversible drug-tolerant states. **No switching or reversion rate (per day or per week) in the abstracts.** Full texts not read. | P (abstract only) | Europe PMC records of the two DOIs |
| Corigliano et al. bioRxiv 2024 | Model of drug-induced persisters in colorectal cell-line data; best schedules had non-zero recovery periods and doses below current clinical doses. Not peer reviewed. | P (record) | DOI 10.1101/2024.11.29.626082 |

Consequence: mu (sensitive-to-resistant switching on drug) and nu (reversion off drug) remain unmeasured. Any claim that depends on them is conditional on an assumed region.

## H. What is still unverified (do not use as a constraint)

1. Stupp 2005 median PFS and 6-month PFS from a primary page (NEJM was blocked).
2. Hegi 2005 unmethylated medians and 2-year survival; Stupp 2009 MGMT-stratified numbers.
3. RANO 2010 full text (nadir rule, 12-week post-radiotherapy rule, steroid and new-lesion rules); RANO 2.0 full text. mRANO figures are summary-level.
4. (Settled by checker report 5 for the two trials above: scans every 3 months. Still unverified: the Macdonald rule's reference point, smallest earlier size or baseline, in RTOG 0525.)
5. (Settled by checker report 5: RTOG 0525 median adjuvant cycles 3 and 4. Read through a page summary; check the PDF before quoting.)
6. Any open paper that reports tumour volume change per TMZ cycle (needed to anchor the kill strength m; none found).
7. Second and third estimates of untreated growth; the full texts of Liau 2017 and Rabe 2020 for any time scale.
