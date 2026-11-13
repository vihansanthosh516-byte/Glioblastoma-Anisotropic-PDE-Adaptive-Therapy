# Integration prior-art check: MRI + tumor biology for GBM survival (2026-11-12)

Method: one independent AI literature search (about 20 queries; PubMed, arXiv and journal pages; several full texts blocked) plus two targeted searches by me for MRI-predicted Neftel cell states. Numbers below come from abstracts or pages and are not all verified. No human expert reviewed this.

## Closest prior work

| Work | What it did | Overlap with BAMD |
|---|---|---|
| Matsumoto et al. 2025, Neuro-Oncology 27(12):3200, doi 10.1093/neuonc/noaf175 | Multi-region DNA methylation deconvolution; stem-like malignant cell proportion; "radiomethylomic" MRI models predicting cell state and survival (UPenn group) | Closest. MRI predicts a cell-state quantity linked to survival. Differences: methylation, one stem-like proportion, same-centre data. Full text not read. |
| Sun et al. 2021, Radiology, doi 10.1148/radiol.2021203281 | MRI signature mapped to RNA pathways ("RadGene"), survival HR in external sets | MRI to RNA programs, then survival. Not Neftel states; no residual split. |
| Chen et al. 2026, IJMS 27:9635, doi 10.3390/ijms27219635 | Radiomics predicts a NANS-high macrophage signature (n = 40, TCIA) | MRI predicts an RNA signature; no survival use of the prediction. |
| Wijethilake et al. 2020 | Radiomics predicts Verhaak subtypes (TCGA) | Subtype, not Neftel states; survival use not verified. |
| Braman et al. 2021 (MICCAI, arXiv 2107.00648); Robinet et al. 2024 (DRIM) | Orthogonal / disentangled multimodal fusion | Related to the residual idea, but in learned embeddings, not a frozen biology-anchored split with a pre-registered test. |
| Kazerooni et al. 2022, Sci Rep, doi 10.1038/s41598-022-12699-z | UPenn IDH-wt n = 516; clinical 0.65, + radiomics 0.70, + genomics panel 0.75 C-index (replication) | Best reference for incremental value tables on UPenn-type data; no RNA. |
| Beig et al. 2018, Sci Rep | TCIA n = 115, radiomics tied to a hypoxia RNA score | MRI to RNA score; small validation. |

Flagged: Zhang X et al. 2025 (Front Biosci) and Nicolasjilwan et al. 2014 appear retracted or withdrawn; not used as evidence.

## Typical C-index ranges reported (abstract level, mostly internal validation)

Clinical only about 0.62-0.67; radiomics only about 0.66-0.74; combined about 0.70-0.79.

## Novelty verdict

- MRI predicts Neftel four-state scores learned from bulk RNA, transferred to an RNA-free external cohort for survival: **no prior work found**; the nearest (Matsumoto 2025) is a different molecular source and one state.
- Splitting MRI into a biology-anchored part and a residual part and testing each for survival: **no direct prior work found** (orthogonal fusion is related).
- Pre-registered incremental-value test: **not found** in this field.

Unverified: full texts of Matsumoto 2025, Sun 2021, Choi 2020, Wijethilake 2020; preprints behind blocked pages. A targeted multiple-kernel-learning search returned no MRI+RNA GBM survival hit.
