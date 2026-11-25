# Track A TMZ data lead: Mayo GBM PDX National Resource (2026-11-24)

Status: lead only. Nothing downloaded. Not checked against prior project use beyond a repo text search (no hit for Vaubel, Sarkaria or Mayo PDX).

Source: Vaubel et al., Clin Cancer Res 2020;26(5):1094-1104, PMC7056576.

Verified from the paper text:
- 96 GBM PDX models; RNA-seq for 68 ("RNA-sequencing (n=68)").
- In vivo TMZ response: orthotopic tumours, mice randomized 5-10 per group, placebo vs TMZ (66 or 50 mg/kg daily x 5), survival endpoint (Kaplan-Meier, log-rank).
- MGMT promoter methylation measured; it correlated with TMZ survival benefit.
- "All models and genomic data are openly available to investigators."

Unverified (from a model summary of the page, must be checked before use):
- About 37 PDX with in vivo TMZ response; per-model data in supplementary Table S7.
- RNA-seq SRA PRJNA548556, WES PRJNA543854, methylation GEO GSE132650; cBioPortal and Mayo portal access.

Why it fits: measured in vivo TMZ response on independent patient tumours, pre-treatment expression, MGMT as the required baseline.

Weak points:
- Small: about 37 treated models if the summary is right; overlap with the 68 RNA-seq models unknown.
- MGMT is a strong known predictor here; a new method must beat MGMT plus expression baselines.
- Prior art on the same resource: Zhao et al., Neuro-Oncology 2019 (PMC6736132), xenograft-based TMZ/RT gene signatures; Nagel et al., Cancer Res 2017 (12 PDX, DNA repair capacity).
- No sister cultures, so the within-only estimator from dev5/dev6 does not apply; the method question changes.

Next: read Table S7 and the SRA run table; count models with both RNA-seq and TMZ survival; check Zhao 2019 for which models it used.

## Supplement check (2026-11-24)

The user downloaded the 8 supplement files from PMC7056576 into `data/external/trackA_mayo_pdx_supp/` (git-ignored). Table S7 (supplement-8.docx) lists 37 PDX lines with standard-therapy survival; 30 have a TMZ survival ratio (median days TMZ / placebo, orthotopic). RNA-seq molecular group is called for 68 lines (Table S6A, supplement-7.xlsx). Lines with both a TMZ ratio and RNA-seq: **29** (MGMT methylated 13, unmethylated 16).

Exposure disclosure: while counting, the median TMZ ratio by MGMT group was printed (methylated 3.53, unmethylated 1.50). No expression data were opened. This cohort is therefore not fully outcome-blind; MGMT-stratified response is known, and any protocol must treat MGMT as a known strong baseline.

Assessment: n = 29 is small. With MGMT already separating response, a new expression method must add to MGMT on 29 models. Power for that is low. Treat as a development or limited test cohort, not a definitive test.
