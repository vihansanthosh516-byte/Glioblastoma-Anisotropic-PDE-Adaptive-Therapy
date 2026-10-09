# Gap analysis: what is missing to reach the end goal (2026-10-09)

End goal: (H-1) measure how forecastable glioblastoma growth is; (H-2) test whether a growth model can place the
radiotherapy margin better than the guideline 15 mm margin, at the same treated volume; (H-3) forecast next-scan volume
with calibrated intervals and an abstain rule. Numbers for our results: `vault/PAPER_FINDINGS_LEDGER.md`.
Literature below was found by web search on 2026-10-09; each entry is [S] (search summary) until the full text is read.

## 1. Where the field is (fresh literature)

| Finding | Source | What it means for us |
|---|---|---|
| Guideline margin is 15 mm around enhancing tumour / cavity (was 20 mm in 2016); no consensus on a T2/FLAIR margin (0-15 mm) | Niyazi et al. 2023, ESTRO-EANO, Radiother Oncol 184:109663 | M0 in H-2 = this plan; confirm PREDICT-GBM `create_standard_plan` uses the same rules |
| PREDICT-GBM: 243 patients, multi-centre; at equal treated volume, U-Net 79.37 +/- 2.08% and GliODIL 78.91 +/- 2.08% enhancing-recurrence coverage; only U-Net and GliODIL significantly beat the standard plan pooled; LMI and nnU-Net were below it; distant recurrences missed by all; gains "modest" | Zimmer et al., arXiv 2509.13360 (v2 Mar 2026); npj Digit Med 2026 | The bar to beat is about 79%. A plain population Fisher-Kolmogorov (like LMI) can fall BELOW the standard plan. Our M1/M2 need per-patient fitting to compete |
| GliODIL: physics-informed discrete loss, infers full tumour-cell map from MRI (+ FET-PET); best recurrence prediction among tested models, 152 test patients | Balcerak et al. 2025, Nat Commun, doi 10.1038/s41467-025-60366-4 | The strongest biophysical competitor; it fits per patient. Its predictions are in the PREDICT-GBM release, so we can compare without running it |
| Prospective pilot: ML-guided personalised dose escalation in new GBM; survival comparison preliminary | Nat Commun 2026, s41467-026-72545-y | Field is moving to clinic; context for the poster, not a method to copy |
| Atlas fibre orientation warped to the patient gives a significant but small Dice gain, mostly in butterfly gliomas crossing the corpus callosum | arXiv 2507.17707 (2025) | Matches our null overall (ledger 2ak). Suggests a pre-registered subgroup: tumours touching the corpus callosum |
| No head-to-head study of patient DTI vs atlas tensor for recurrence | search 2026-10-09 | Open question we could answer with UPENN-GBM DTI (317 patients) if recurrence or follow-up masks exist |
| No shared benchmark for next-scan glioma prediction; treatment-aware diffusion model predicts future masks with uncertainty | IEEE TMI 2025 (PubMed 40031286); arXiv 2509.10824 | H-3 has no standard competitor; persistence and the geometric rule stay the baselines |
| Conformal intervals for longitudinal lesion size; time-varying nonconformity score for irregular follow-up; exchangeability breaks under site shift (weighted conformal) | arXiv 2609.21197 (2026); Tassopoulou et al. (OpenReview 2025); arXiv 2407.19938 | H-3 method: split conformal on MU, time-varying score (our intervals are irregular), report coverage on LUMIERE as a shifted site |

## 2. What we are missing

### H-2 (the possible strong positive)
1. **PREDICT-GBM data not downloaded** (Hugging Face LZimmer/PREDICT-GBM, about 25.5 GB). Needs the author's OK. Amendment 6 is already committed, so the download is allowed by our own rules.
2. **Per-patient fitting.** Our models use one population cell per fold. The PREDICT-GBM result says population models (LMI) can lose to the 15 mm margin. Fix: fit lambda = sqrt(D/rho) per patient from the pre-op scan (enhancing edge vs FLAIR edge; Konukoglu 2010 eikonal idea; Lipkova 2019 Bayesian) using the exact solver. Must be added to Amendment 6 as a model before any test-set run.
3. **Adapter to PREDICT-GBM inputs** (their tissue maps, SRI24/MNI space, our atlas tensor registered to that space) and use of `topk_plan` for equal volume.
4. **Compare against published predictions** (U-Net, GliODIL, LMI files in the release): report our model next to them on the same patients, not only M2 vs M0.
5. **Distant recurrence**: nobody captures it. Report the share of distant recurrence separately; do not hide it.

### H-3
6. **Conformal intervals and abstain rule not built.** Split conformal on MU development folds, time-varying nonconformity score, then frozen test on LUMIERE; report coverage at 80/90%.
7. **LUMIERE volumes are noisy** (tools disagree about follow-up change by more than the change, ledger 2ab). State that the LUMIERE test bounds what any forecast can show.

### Direction question (follow-up to ledger 2ak)
8. **Patient DTI test.** UPENN-GBM has DTI for 317 patients (CC BY 4.0). Check first whether it has follow-up or recurrence scans; if not, patient DTI can only be tested where both DTI and a later mask exist.
9. **Corpus-callosum subgroup**, pre-registered before looking.

### H-1
10. **True scan-rescan noise**: not available openly (RIDER Neuro MRI and QIN-GBM are dbGaP controlled access). Measured substitutes: tool disagreement (2ab, 2ag), expert correction (2ai). Simulated rescans (scanner perturbation, then re-segmentation) are still open.

### Paper and poster
11. Remove every "direction helps" claim (ledger 2ak shows it was a solver artefact).
12. Write the solver-artefact finding as its own result, with the Selling method and the GPU rerun.

## 3. Order of work (recommended)
1. Fix old claims (item 11). Small, urgent.
2. Amendment 7: per-patient fitting model and the corpus-callosum subgroup, committed BEFORE the PREDICT-GBM download.
3. PREDICT-GBM download (needs OK), adapter, development-set runs on TUM only.
4. Freeze (`configs/h2_frozen.yaml` + tag `h2-frozen`), then test on LUMIERE + RHUH.
5. H-3 conformal (CPU, can run in parallel).
6. UPENN-GBM DTI check (item 8).

## Sources
- Niyazi et al. 2023 ESTRO-EANO: https://pubmed.ncbi.nlm.nih.gov/37059335/
- PREDICT-GBM: https://arxiv.org/abs/2509.13360 ; https://www.nature.com/articles/s41746-026-03194-0
- GliODIL: https://www.nature.com/articles/s41467-025-60366-4
- ML-guided dose escalation pilot: https://www.nature.com/articles/s41467-026-72545-y
- Butterfly glioma fibre-tract model: https://arxiv.org/pdf/2507.17707
- Treatment-aware diffusion model: https://pubmed.ncbi.nlm.nih.gov/40031286/
- Multi-task diffusion glioma progression: https://arxiv.org/html/2509.10824v1
- Conformal lesion-size forecasting: https://arxiv.org/pdf/2609.21197
- Robust conformal volume estimation: https://arxiv.org/pdf/2407.19938
- Konukoglu et al. 2010 IEEE TMI: https://nmr.mgh.harvard.edu/node/4758
