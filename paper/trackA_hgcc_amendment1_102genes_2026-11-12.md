# Track A HGCC bridge, amendment 1: 102-of-119 exact-mapped score (2026-11-12)

**Frozen before any TMZ response magnitude was read.** The commit adding this file is the freeze. The gate check (`paper/trackA_hgcc_gate_check_2026-11-12.md`, commit bf0e18d) used IDs, gene names, and TMZ presence only. The user chose this amendment on 2026-11-12.

Base protocol: `paper/trackA_hgcc_tmz_response_bridge_2026-11-12.md`. Everything there stands except the changes below.

## Changes

1. **Gene set.** Use the 102 of the 119 frozen panel genes whose symbols appear exactly once in `HGCC_RNA_expression.txt`. Omitted (17, no exact symbol): AMY2B, CFAP44, EXTL3-AS1, KANTR, LENG8-AS1, LINC00907, LINCR-0001, MRNIP, NBPF11, PRKCQ-AS1, PSMA3-AS1, PTOV1-AS1, PTOV1-AS2, RSRP1, SLC25A25-AS1, TWF2, ZNF547. No aliases, substitutions, or imputation. Original Ntafoulis directions, equal weights. The direction counts (resistance/sensitivity) are whatever the 102 contain and are recorded in the result. Controls use the same counts.
2. **Score.** Within each expression profile, percentile-rank (average ties, rank / n) all genes whose symbol is non-empty and occurs exactly once. Score = mean percentile of resistance-associated panel genes minus mean percentile of sensitivity-associated panel genes.
3. **Join.** Exact `U####` IDs only (83 origins). Drug-file columns with suffixes (`b`, `replicated`, `U3065_<n>`) are not used, because the source does not document that they are the same assay and dose range. One origin = one `U####` ID. Expression profiles sharing an ID (11 IDs have two) are scored separately and the scores averaged.
4. **Outcome.** The `Temozolomide` row of `HGCC_drug_response_AUC.txt`, exact-ID columns only. Higher AUC = more viability under drug.
5. **Primary test.** Spearman rho between origin score and origin TMZ AUC. One-sided permutation test in the predefined positive direction: 100,000 permutations of AUC across origins, NumPy `Generator(PCG64DXSM(20261112))`, p = (1 + count(rho_perm >= rho_obs)) / (1 + 100,000). Percentile 95% bootstrap interval over origins, 10,000 resamples, `PCG64DXSM(20261113)`.
6. **Decision (added now, before outcome access, because the base protocol named no alpha).** "Supports an in-vitro association" if rho > 0 and permutation p < 0.05. "Opposite direction" if rho < 0 and the bootstrap interval excludes 0. Otherwise "no association shown (null or inconclusive)". No other pass condition.
7. **Secondary, descriptive only.** Spearman rho of MGMT expression (exact symbol `MGMT`, origin mean) with TMZ AUC. 1,000 direction-preserving expression-matched control panels built exactly as base protocol step 5 (gene mean and SD across origin-averaged profiles, standardized; for each panel gene in a seeded random order, pick uniformly among its 20 nearest unused non-panel genes; seed `PCG64DXSM(20261114)`), built before the outcome is loaded; report the fixed score's percentile among control rhos.
8. **Outputs.** `output/trackA_hgcc_tmz_102/result.json` with aggregate statistics only, file hashes, and the amendment hash. No per-line values are written or committed (data terms unclear).

## Claim boundary

A positive result would show only that this fixed, reduced expression score ranks measured in-vitro TMZ viability in one previously published cell-line panel. It is a secondary reanalysis of data Ju et al. (2025) already modelled. It is not patient drug resistance, not a new method, not linked to any MRI patient, and it cannot change A5's failure or reopen Track C.
