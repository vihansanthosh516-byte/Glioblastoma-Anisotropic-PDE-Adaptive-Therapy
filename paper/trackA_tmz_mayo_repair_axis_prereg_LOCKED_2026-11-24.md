# Track A TMZ: repair-axis score on Mayo PDX (LOCKED pre-registration, 2026-11-24)

Locked before any Mayo expression value is read. Code: `src/trackA_tmz_mayo_transfer.py` (hash recorded in the result). One run: `py -3.14 src/trackA_tmz_mayo_transfer.py`. The script refuses to overwrite `output/trackA_tmz_mayo_transfer/result.json`. No retuning after the run.

Supersedes the draft `trackA_tmz_mayo_transfer_prereg_DRAFT_2026-11-24.md`, which was not locked because its new step was inert.

## Change of measure (recorded, user-approved)

Track A drug question stays: predict measured TMZ response from untreated tumour expression with a new method that beats existing methods. The measure changes from in-vitro viability to **in-vivo survival benefit**: TMZ survival ratio (median days TMZ / placebo) in orthotopic PDX (Vaubel 2020, Table S7). The user approved this change in chat on 2026-11-24 ("continue with whatever you think is right for the drug part"; "let's do the Mayo Clinic one").

## Data

- Test: 27 Mayo PDX with a TMZ ratio (Table S7) and cBioPortal RNA-seq (`data_mrna_seq_rpkm.txt`, SHA-256 48a69b409f5dbdf031ec93918bd5c2a955bbf8ad87334c1c10c5aaa465a283cc, ODbL). Expression = log2(RPKM + 1). MGMT status = Table S7 M/U.
- HGCC (83 origins) is used only to train comparators 5-6.
- Known exposure before lock: median TMZ ratio by MGMT group (M 3.53, U 1.50). No Mayo expression values seen.

## Primary candidate C2: repair-axis score

Fixed formula, no training on any outcome:

score = z(-MGMT expression) + z(mean of z(MSH2), z(MSH6), z(MLH1), z(PMS2))

z = standardize across the 27 models. Rationale: TMZ lesions are removed by MGMT (resistance), and TMZ cytotoxicity needs intact mismatch repair (MMR loss gives resistance). Expression-only design check on HGCC (no outcome): score vs MGMT-only correlation 0.59, so the MMR term changes predictions materially.

## Comparators (same 27 models)

1. MGMT status (M/U).
2. MGMT expression alone (**ablation**: C2 without the MMR term).
3. RNA-seq molecular group, best single indicator (favours the comparator).
4. Ntafoulis-derived 102-gene score.
5. HGCC-trained ridge on MGMT-residual TMZ AUC (2,000 variable genes, rank-transformed).
6. HGCC-trained ridge on raw TMZ AUC.
7. Zhao et al. 2019 signature: not run (gene list and weights not obtained).

Random control: 1,000 scores with the same formula, the four MMR genes swapped for four random genes (MGMT kept).

## Primary test

Partial Spearman between score and log TMZ ratio, adjusting for MGMT status. One-sided permutation p, 10,000 permutations within MGMT groups.

**PASS only if all hold:** C2 partial rho > 0; p < 0.05; C2 partial rho above every comparator 1-6; C2 above at least 95% of random swaps. Otherwise FAIL, reported in full.

## Code checks before lock (generated data only)

- Planted repair-axis signal: PASS. No signal: FAIL. MGMT-only signal: FAIL (ablation wins).
- 40 null datasets: 0 false passes.

## Weak points

- n = 27. About 80% power only for a partial rho near 0.45. A fail means "not shown".
- MMR expression is a proxy; MMR loss in TMZ resistance is mostly by mutation (Table S5A lists MMR-mutant PDX), which this score does not use.
- PDX RNA-seq may include mouse stroma reads.
- Novelty is limited. MGMT and MMR roles in TMZ response are established biology (e.g., Nagel et al. 2017 DNA repair capacity in 12 PDX). The new part is a fixed, untrained repair-axis expression score tested for in-vivo benefit beyond MGMT status against trained transfer models. A full prior-art check is still needed before any novelty claim.
- Measure changed from in-vitro to in-vivo (recorded above).
