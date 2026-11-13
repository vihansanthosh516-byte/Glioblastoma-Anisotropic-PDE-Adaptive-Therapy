# Integration protocol: BAMD Gate 1b on region-matched biopsies — locked 2026-11-12

Same question as BAMD (`paper/integration_protocol_bamd_2026-11-12.md`): does pre-operative MRI carry tumor biology? Gate 1 on TCGA failed (`paper/integration_bamd_gate1_result_2026-11-12.md`); it is spent and is not re-run or changed. This is a **new test on new data**: MRI values measured at the exact biopsy site, matched to RNA from that biopsy. Written before any value in the data was analysed (only column names, row counts and the CE/NE label counts were seen).

## Data

Hu LS et al. 2023, Nat Commun 14:6066 (doi 10.1038/s41467-023-41559-1). Public tables from github.com/HuLiLab/Multi-Regional-GBM-Imaging-and-Genetics, folder "MEM Modeling Pathway Signatures", downloaded 2026-11-12 to `data/external/hu2023/`:
- `MultiReg_imaging_pathclass.csv`: 95 MRI-localized biopsies, 25 patients. MRI at biopsy: T2, RCBV, STDRCBV, MD, FA, SPGRC (T1 post-contrast), EPI; CE/NE label. Biology: Garofano pathway scores GPM (glycolytic/plurimetabolic), MTC (mitochondrial), NEU (neuronal), PPR (proliferative/progenitor).
- The advanced-imaging table (107 biopsies) is not used (its perfusion parameters are not in UPENN).
Prior use: Hu et al. reported per-variable mixed-model associations. They did not report leave-one-patient-out prediction or transfer to an external cohort (not verified from full text).

## Gate 1b (outcome-free; no survival anywhere)

- Inputs X: the 7 MRI values + CE indicator (8 features). Targets Y: the 4 pathway scores. Rows with any missing X or Y dropped (count reported).
- Model: ridge, alpha chosen by inner leave-one-patient-out CV from {0.1, 1, 10, 100}, features standardized on training patients.
- Evaluation: leave-one-patient-out (LOPO); out-of-patient predictions pooled; per target Spearman rho.
- Null: 2,000 permutations. Each permutes the biology rows **within** the CE group and **within** the NE group (keeps the CE/NE balance, breaks the MRI-biology link at the same site), then re-runs the full LOPO.
- **Pass** for a target: rho >= 0.25 and Holm-adjusted p < 0.05 over the 4 targets.
- Secondary (descriptive): the same with the CE indicator removed (does local MRI add beyond CE vs NE?).

## If any target passes (stage 2; separate locked note before running)

Transfer needs the same MRI types in UPENN-GBM (T2, T1 post-contrast, DSC rCBV, DTI MD/FA are distributed by TCIA for UPENN-GBM; not yet downloaded). Stage 2 would compute tumor-region median values in UPENN, apply the frozen biopsy-level map, and test survival exactly as H1 in the BAMD protocol. Scaling differences between Hu's units and UPENN's maps are a known risk; stage 2 will define within-patient normalization before any UPENN value is read.

## If nothing passes

Report it. UPENN stays unopened.

## Weaknesses

25 patients; biology is a 4-class pathway framework (Garofano 2021), not Neftel states; values are as published (normalization not re-checked); this re-analyses published data, so the novel part is prediction in held-out patients and the planned transfer.
