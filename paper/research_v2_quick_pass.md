# Research prompt v2: quick pass (INTERIM)

Scope: 8 web searches, search snippets only. NO full paper was read. Every row is "abstract/snippet level". This does NOT cover all 13 questions. Treat as leads, not findings.

## A. Enhancing recurrence

| Method / study | Link | Key number | Test design | Access | Full paper? | T4? |
|---|---|---|---|---|---|---|
| Tran et al. 2025, AI target from metabolic + diffusion MRI (npj Digit Med) | https://doi.org/10.1038/s41746-025-01861-2 | Specificity 0.89 vs 0.79 (p=0.004); sensitivity 0.92 vs 0.95 (p=0.10) vs uniform 2 cm | Not verified from snippet; UCSF group. Possible single-centre | Paper open; data not verified | No | Not verified |
| ISMRM 2025 abstract, Kukreja (CNN/ViT, pre-surgery MRI) | https://archive.ismrm.org/2025/2561.html | Claims better inclusion/exclusion than standard; no number seen | Abstract | n/a | No | Not verified |
| Post-op + white-matter pathlength CTV (ASTRO 2025 abstract) | https://am25.astro.org/sessions/pqa-03-21637/planning-assessment-of-pre-rt-white-matter-mapping-enhanced-ctv-for-glioblastoma-radiation-th-107567 | 125 train / 12 test | Tiny test set | n/a | No | Not verified |

Why hard (Q1), partial:
- Pattern studies: most first recurrences are central/in-field. Smaller GTV had more distant recurrence (Freiburg, n=97). https://doi.org/10.1186/1748-717x-8-97
- NOT found in this pass: pseudoprogression, registration error, residual tumor effect. Not verified.
- Q2 (equal-volume gain), Q3 (post-op MRI gain), Q4, Q5, Q6: NOT researched. Not verified.

Weakest point: Tran 2025 compares to a 2 cm margin on sensitivity/specificity, not equal-volume coverage. Not the same metric as yours.

## B. New locked test set

| Dataset | Link | Key fact | Access | Verdict |
|---|---|---|---|---|
| RHUH-GBM | https://arxiv.org/pdf/2305.00005 | n=40, pre-op, early post-op, recurrence | TCIA limited licence | Already in your training set. Not eligible. |
| UCSD-PTGBM | https://escholarship.org/content/qt66w5w496/qt66w5w496.pdf | 178 patients, 243 timepoints, post-treatment | CC BY 4.0 (per TCIA notes) | Pre-op scans NOT confirmed. Likely not eligible. |
| BraTS 2025 Task 1 (pre + post treatment) | https://zenodo.org/records/15094823 | Train/val on Synapse; hidden test closed | Synapse registration | Recurrence/enhancing-label pairing NOT verified. Check Synapse. |
| GLIS-RT, IvyGAP, CPTAC-GBM, TCGA-GBM, ReMIND | none checked | Not searched | n/a | Not verified |

## C. Growth and digital twin

| Study | Link | Key number | Design | Full paper? |
|---|---|---|---|---|
| BMC Cancer 2025, two-species reaction-diffusion data assimilation | https://bmccancer.biomedcentral.com/articles/10.1186/s12885-025-14557-3 | n=21 calibrated on weekly MRI | Held-out forecast and baseline NOT confirmed | No |
| Counseling model, 132 intervals, 46 recurrent GBM patients | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12515658/ | 86% of intervals: at least one scenario within 20% of volume | Best-of-scenarios, so lenient. No persistence baseline seen | No |

- Q8: no paper confirmed that beats a persistence baseline on held-out scans. Supports your null result. Not verified for the full literature.
- Q9: not researched.
- Q10: NASEM 2023, https://www.nationalacademies.org/read/26894/chapter/16 (DOI 10.17226/26894). Secondary sources say a twin needs bidirectional data flow and VVUQ (verification, validation, uncertainty quantification). Exact wording NOT retrieved. Fetch chapter before quoting.

## D. Treatment effect

| Dataset / study | Link | Key fact |
|---|---|---|
| Burdenko B-GBM-PD | https://wiki.cancerimagingarchive.net/x/YgywC (per search) | 180 pts, RTDOSE/RTPLAN subset. You rejected it already. |
| RHUH-GBM | arXiv link above | RT dose, systemic therapy, extent of resection. In your training set. |
| Mayo hypoxia + radiation model | https://mayoclinic.elsevierpure.com/en/publications/a-patient-specific-computational-model-of-hypoxia-modulated-radia | Radiosensitivity from 2 MRIs; proof of principle only |
| Leder et al. Cell 2014 | https://www.jci.org/references/scholar/35152/B49 (listing only) | Mouse, not patients |

- Q11: TMZ dates and steroids not found in any dataset snippet.
- Q13: no patient-validated virtual schedule comparison found in this pass. Not proven absent.

## Best next steps
1. Fetch BraTS 2025 Synapse page and check for pre-op + recurrence pairs with an enhancing label.
2. Fetch Tran 2025 full text; check its test design and whether it reports equal-volume coverage.
3. Fetch NASEM chapter 16 and quote the definition for the digital-twin wording.

## Could NOT verify
Everything not linked above. All Section A Q2-Q6. Pseudoprogression/registration/residual evidence. Whether any listed paper tests only on its own centre. All numbers are from search snippets.

---
# Update: 3 fetches done (INTERIM, still partial)

## Tran 2025 (full text read, first 100k of 114k chars)
https://pmc.ncbi.nlm.nih.gov/articles/PMC12331921/
- 101 patients, 2002-2012, split 67/16/18 -> test about 18 patients. Single centre (UCSF). Same-centre hold-out. No external test.
- Input is POST-operative, pre-RT (within 1 week before RT). Not pre-op. Not your task.
- Target mixes enhancing (CEL) and non-enhancing (T2L) progression. No enhancing-only result.
- No volume-matched comparison. Metrics are sensitivity, specificity, PCC.
- vs RTOG 2 cm: sens 0.92 vs 0.95 (p=0.09), spec 0.89 vs 0.79 (p=0.004). vs EORTC 1.5 cm: sens 0.92 vs 0.74.
- Text mismatch: abstract p=0.10, results p=0.09. Data private. Code on request only.
- Verdict: weak evidence for your question. Cannot reproduce. Small n.

## NASEM chapter 16 (read)
https://www.nationalacademies.org/read/26894/chapter/16
- WARNING: this chapter is workshop proceedings (speaker remarks), NOT the consensus report text. Do not cite it as the NASEM definition.
- Speaker definition (Baron, p. 165): "a dynamic virtual copy of a physical asset, process, system, or environment."
- Speaker point (Farrar, p. 165): limits "define the domain in which one has confidence."
- The terms "VVUQ" and "fit-for-purpose" are not in this page. Get the consensus report (DOI 10.17226/26894) main chapters for the VVUQ and bidirectional-flow wording.

## BraTS 2025 Task 1 (partial)
https://zenodo.org/records/15094823
- Listed as "Pre- and Post-Treatment Adult Glioma". CC BY 4.0 record.
- 2025 case count, labels, and Synapse terms NOT found. 
- Precursor BraTS 2024 post-treatment set had labels ET, NETC, SNFH and resection cavity RC; about 2200 training samples (https://arxiv.org/pdf/2411.17617). That is post-treatment only.
- Open question: does 2025 give pre-op and later recurrence for the SAME patient? Not verified. Check the Synapse project directly (needs your login).
