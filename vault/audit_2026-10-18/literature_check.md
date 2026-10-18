# Literature check: what published work says about the open problems

Scope: full texts read for three papers (open-access PDFs fetched and searched). A few other items surfaced in search snippets only (for example a preoperative-postoperative recurrence-risk model and TumorFlow); they were NOT read in full and are not relied on. Verify each citation yourself before using it.

## 1. PREDICT-GBM platform paper (npj Digital Medicine 2026, 9:686; 243 patients, multicentre)
doi 10.1038/s41746-026-03194-0. This is the benchmark your H-2 and H-3 study is built on.
- The authors describe the gains over the guideline plan as modest. Mean enhancing-recurrence coverage: U-Net 79.37 +/- 2.08%, GliODIL 78.91 +/- 2.08% (BH-adjusted p 2.9e-5 and 1.0e-3). All plans are iso-volumetric to the guideline plan.
- Ablations (the authors call them minor gains): retraining on only the three PREDICT-GBM datasets gave 79.26%; tissue maps and tumour segmentation only, no ensemble, gave 78.86%. Extra training cases and extra input channels together are worth well under one percentage point. Removing the distance-based post-processing cost 0.67 points. A white-matter-weighted prior reduced coverage (-0.72 points at weight 0.2, -14.7 points at 0.5).
- Sensitivity to labels: on RHUH (n=40, a cohort in your H-2 test set), with the expert-corrected recurrence segmentations in place of the BraTS-derived ones, the gap narrowed substantially: SOC 79.24 +/- 5.05%, U-Net 80.47 +/- 4.95%, GliODIL 79.23 +/- 5.01%. The authors conclude the model advantage is sensitive to annotation quality.
- Failure cases they report: distant recurrences, very small initial tumours (iso-volumetric redistribution has little room), and multifocal tumours (the guideline plan outperformed model plans). They also note frequent ties under saturated coverage.

## 2. Metabolic and diffusion MRI target definition (npj Digital Medicine 2025; 101 patients)
doi 10.1038/s41746-025-01861-2. Reports that lesion-size-weighted loss and metabolic (MRSI) inputs improved the model. The comparison is sensitivity and specificity against a 2 cm T2 expansion, not iso-volumetric enhancing-recurrence coverage, and MRSI coverage was a limitation in the retrospective data. Not directly transferable to your data (you have no MRSI).

## 3. Glioma growth models after resection (IEEE TMI 2023)
doi 10.1109/tmi.2023.3298637. Fourteen untreated low-grade glioma patients. DTI-informed anisotropic diffusion showed a small but significant improvement over homogeneous isotropic diffusion. Different tumour grade, no treatment, n=14, so it neither supports nor contradicts a null on treated GBM.

## What this means for your results
1. **Your enhancing-recurrence nulls match the size of effects in the literature.** Published iso-volumetric gains are around a percentage point or two, shrink under expert-corrected labels, and the 15 mm plan already saturates in many patients (your 32 of 49 UCSF). A +0.009 locked-test result with a CI crossing 0 is not evidence of a bug.
2. **A DTI or white-matter prior has not been shown to help on treated GBM.** A published ablation found a white-matter-weighted prior lowered coverage. Your exact-solver result (-0.0004) and the 14-patient LGG result do not conflict.
3. **Documented but small levers:** more training data, extra channels and ensembling (each worth under a point in the benchmark), and a size-weighted loss (promising, but from a different design). Any of these must be tried on development data first, under a declared rule, because the H-2 and H-3 test sets are now dev (A12.1).
4. **A label-quality analysis is well motivated.** The benchmark authors show the advantage depends on annotation quality, and you have already measured segmentation disagreement (script 104) and UPENN rater noise (script 109). A declared, exploratory re-score on RHUH with expert-corrected labels is a legitimate extension. It is post hoc because RHUH was used in the H-2 test.
5. **Growth forecasting:** I found no published model that beats persistence on real treated-GBM longitudinal MRI. That fits your results and the README's "unverified" list; it is not proof that none exists.

## Limits of this check
- Three papers read in full; this is not a systematic review.
- Metrics, cohorts and designs differ from yours; the numbers are context, not benchmarks to match.
- Post-op MRI, perfusion and PET inputs, and small-data methods (topics in your research prompt) were not covered.
