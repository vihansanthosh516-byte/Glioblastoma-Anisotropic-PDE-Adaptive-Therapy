# Track A TMZ: frozen repair-axis score, human test on GSE7696 (LOCKED, 2026-11-24)

Locked before the survival analysis. Code `src/trackA_tmz_gse7696_test.py` (hash in result). One run: `py -3.14 src/trackA_tmz_gse7696_test.py`. The script refuses to overwrite `output/trackA_tmz_gse7696/result.json`.

## Question

Does the same frozen repair-axis score that passed on Mayo PDX (`trackA_tmz_mayo_repair_axis_result_2026-11-24.md`) predict which patients **benefit from adding TMZ** to radiotherapy, beyond MGMT methylation and age?

## Data

GSE7696 (Murat et al., J Clin Oncol 2008), Affymetrix U133 Plus 2.0, log2 values as deposited. Matrix SHA-256 6b4cb5d3...caa25; GPL570 annotation SHA-256 d7cd4435...f8bde8f394.

Eligibility (fixed in code): disease status "GBM" (primary), treatment "TMZ/radiotherapy" or "radiotherapy", MGMT M or U, first sample per patient in GEO order. Result: **68 patients, TMZ+RT 42, RT 26** (M/U: TMZ+RT 21/21, RT 15/11). Counts were taken from labels only.

Prior exposure: an earlier metadata listing showed a few survival times among sample characteristics; no survival analysis was done. GSE7696 has not been used in this project before (repo search).

## Score (frozen, unchanged from Mayo)

score = z(-MGMT) + z(mean z(MSH2, MSH6, MLH1, PMS2)), z across the 68 patients. Gene expression = mean of all GPL570 probes whose symbol list contains the gene (MGMT 204880_at; MSH2 209421_at; MSH6 1569034_a_at, 202911_at, 211449_at, 211450_s_at; MLH1 202520_s_at; PMS2 209805_at, 221206_at).

## Primary test

Cox model: overall survival ~ arm (TMZ+RT = 1) + score + arm x score + MGMT status + age. One-sided Wald test that the arm x score coefficient is negative (higher score, larger TMZ benefit).

**PASS only if all hold:** coefficient < 0; one-sided p < 0.05; its z is more negative than the same model with MGMT expression alone (ablation); its z beats at least 95% of 1,000 random swaps (four MMR genes replaced by four random single-symbol genes, MGMT kept). Otherwise FAIL.

## Code checks (generated data only)

Planted interaction (n = 300): PASS. No interaction: FAIL. 30 null datasets shaped like this cohort (68 patients, 42/26 split): 0 false passes.

## Weak points

- **Low power.** Interaction tests need far more patients than main effects; with 68 patients and 26 in the control arm, only a very large interaction can pass. A fail means "not shown".
- Randomization status of every sample is not verified from GEO; the cohort mixes EORTC/NCIC trial patients and others per the source paper. Treatment may not be fully randomized in this subset.
- PMS2 probes also match PMS2CL (pseudogene).
- Different platform (microarray) and species context (human vs PDX) from the Mayo test, which is the point of this shift test.
