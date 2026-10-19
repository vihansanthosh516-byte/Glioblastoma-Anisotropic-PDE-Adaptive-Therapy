# Literature search: digital-twin parts (2026-10-17)

Status: FIRST PASS. Four web searches only. Every row below came from a search snippet, not from reading the full paper. Open each link and check numbers before you cite them. All third-party figures, study designs, dataset counts, fields, access and license descriptions remain provisional; exact primary-source tables/sections have not been recorded here.

## Part 1: recurrence target
| Item | Link | Key number (from snippet) | Test design | Access | Fits compute |
|---|---|---|---|---|---|
| Multiparametric DL vs uniform 2 cm CTV (diffusion + metabolic MRI) | https://www.nature.com/articles/s41746-025-01861-2 | Specificity 0.89 vs 0.79 (p=0.004); sensitivity 0.92 vs 0.95 (p=0.10) | Single centre, retrospective, small n. Not equal volume. Not enhancing-only. | Paper open. Data not checked. | Unknown (needs metabolic MRI) |
| Earlier abstract, same idea | https://doi.org/10.1093/neuroncol/noad179.0215 (PMC10640264) | 72 patients (54/18). 2 cm CTV sens 0.833 vs 0.795 | Conference abstract. Preliminary. | Open | n/a |
| Transformer-DenseNet 1-year recurrence, HGG | https://pmc.ncbi.nlm.nih.gov/articles/PMC13156239/ | not read | Predicts recurrence, not CTV replacement | Not checked | Not checked |
| Multicentre protocol: DL to predict recurrence spread | https://pmc.ncbi.nlm.nih.gov/articles/PMC13395678/ | none (protocol) | Protocol only | Not checked | n/a |

The snippets do not establish a directly comparable equal-volume enhancing-recurrence test. The local H-2 result was mean enhancing-recurrence coverage gain +0.0088 versus the 15 mm margin (98 patients; interval includes zero; ledger 2ax; `output/predictgbm/final/test_results.json:PRIMARY_final_minus_standard_enh`). This is distinct from the geometric growth baseline in ledger 2r. The short search found no directly comparable result; it cannot rank these methods.
Item 2 (new datasets) and item 3 (why enhancing recurrence is hard): NOT searched yet.

## Part 2: growth forecast
| Item | Link | Note |
|---|---|---|
| Bayesian glioma growth, tissue heterogeneity | https://arxiv.org/pdf/2209.12089 | Calibration method. No persistence baseline seen in snippet. |
| Bayesian personalised RT design, multimodal | https://pmc.ncbi.nlm.nih.gov/articles/PMC7170051/ | Single time point. FET-PET. |
| Proliferation-saturation index, adaptive Bayes | https://www.tandfonline.com/doi/full/10.1080/09553002.2019.1589013 | R2 0.84, c-index 0.89. Lung cancer data, NOT glioma. |
| LUMIERE: 91 GBM, 638 study dates | https://www.nature.com/articles/s41597-022-01881-7 | Already used by us. |
| CFB-GBM v2.0: 264 GBM, GTV at t0-t2, RANO | https://arxiv.org/html/2608.17884v1 | Possible new data. Dose/TMZ fields NOT confirmed. |
| NASEM 2023 digital twin report | https://nap.nationalacademies.org/read/26894/chapter/3 | Twin = model + physical counterpart + two-way data flow, with verification, validation, uncertainty quantification. Exact wording NOT read. |

Not found: a paper where a patient-specific GBM growth model beat a no-change baseline on held-out real scans. This does NOT prove none exists. Search was too short. Locally, the PDE forecast did not significantly beat persistence, whereas the population-rate geometric baseline had positive GBM-only mean gain +0.0149 (ledger 2w); those are different models.
Honest term for our work: "patient-specific simulation", not "digital twin" (no two-way data flow, no validated forecast).

## Part 3: treatment effect
| Item | Link | Note |
|---|---|---|
| SAILOR (Oslo): RT dose map, 3-19 follow-ups, 60 Gy + TMZ | not found by search | 25 patients. Public status NOT verified. Find the source. |
| UCSF-ALPTDG: treatment history incl. chemo | https://doi.org/10.1148/ryai.230182 | Dose not confirmed. We already use it. |
| Rectified-flow post-RT MRI prediction | https://arxiv.org/pdf/2603.08385 | Uses pre-RT priors. Not read. |
| 3D-GlioPREDICT | https://arxiv.org/pdf/2606.05113 | Not read. |

Not searched: MGMT/IDH/ADC/rCBV biomarkers, validation of virtual-treatment schedules.

## Next steps
1. Open the CFB-GBM v2.0 paper and the SAILOR source. Check dose, TMZ dates and access.
2. Search Part 1 items 2 and 3, and Part 3 items 2-4, in a second pass.
3. Read the Nature 2025 paper. Check the lesion-size-weighted loss idea against our soft top-k loss.

## Could not verify
- Any held-out growth win over persistence in real GBM.
- SAILOR public access.
- CFB-GBM dose and TMZ fields.
- NASEM exact wording on validation.
- Whether the PSI result holds in glioma.

## Update 2026-10-17: locked-test candidates (Part 1 item 2). Snippets only, NOT verified.
- CFB-GBM v2.0 (https://arxiv.org/html/2608.17884v1, TCIA, CC BY 4.0): REJECTED for H-4. t0 is ~1 week before radiochemo (likely post-op); one merged GTV label, no separate enhancing; many GTVs are nnU-Net output.
- Burdenko-GBM-Progression (https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=133073470): 180 pts, 1-8 follow-ups, GTV from RT plans. Pre-op baseline and enhancing recurrence labels NOT confirmed. Check first.
- Duke longitudinal + landmarks (https://zenodo.org/records/20548528): 61 pts, pre-op and post-recurrence pairs, auto-segmentation, made for registration. Labels weak.
- RHUH-GBM, LUMIERE, UCSF, TUM, MU: already used. Not new.
- OpenNeuro ds007045 (337 pts): pre-op only, no recurrence scan.
- UCSD-PTGBM (CC BY 4.0): post-op glioma. Pre-op pairing not confirmed.
- BraTS 2025: pairing and access NOT confirmed.
