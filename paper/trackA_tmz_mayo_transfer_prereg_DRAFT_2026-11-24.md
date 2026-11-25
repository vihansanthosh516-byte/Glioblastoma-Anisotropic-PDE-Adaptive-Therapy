# Track A TMZ: in-vitro to in-vivo transfer beyond MGMT (DRAFT pre-registration, 2026-11-24)

Status: DRAFT. Not locked. Mayo PDX expression has not been downloaded or opened. Lock needs user sign-off, then commit and push before any Mayo expression is read.

## Question

Does an expression model trained on measured in-vitro TMZ response (HGCC, 83 lines) predict in-vivo TMZ survival benefit in independent orthotopic PDX (Mayo, 29 models) **beyond MGMT promoter methylation**?

## Data

- Training (development only): HGCC expression and TMZ AUC, local in `data/hgcc/`. HGCC is spent as a test for the 102-gene score; here it is used only for training.
- Test: Mayo GBM PDX (Vaubel 2020). Outcome = TMZ survival ratio (median days TMZ / placebo), Table S7. Expression = cBioPortal `gbm_mayo_pdx_sarkaria_2019_rna_seq_mrna` (RPKM). MGMT = Table S7 PDX MGMT (M/U). Eligible: the 29 models with both a TMZ ratio and RNA-seq.
- Known exposure: median TMZ ratio by MGMT group was seen (M 3.53, U 1.50). Expression has not been seen.

## Candidate (one primary)

C1, MGMT-residual rank transfer:
1. Genes in both platforms. Within each sample, convert expression to ranks (removes platform scale).
2. In HGCC, regress TMZ AUC on MGMT expression; keep the residual as the target.
3. Ridge regression on the 2,000 most variable shared genes (variance computed in HGCC only); penalty chosen by 5-fold CV in HGCC, grouped by patient origin.
4. Apply the frozen model to Mayo ranks to get one score per model. Direction set so a higher score means more TMZ benefit.

## Comparators (all scored on the same 29 models)

1. MGMT status alone (M/U).
2. MGMT expression (Mayo RNA-seq).
3. RNA-seq molecular group (classical / mesenchymal / proneural).
4. Ntafoulis-derived 102-gene score (earlier Track A panel).
5. Same pipeline without the MGMT residual step (ridge on raw AUC).
6. 1,000 random gene sets of 2,000 genes through the same pipeline.
7. Zhao et al. 2019 xenograft TMZ signature, if its gene list and weights are available before lock; otherwise listed as not run.

## Primary test

Partial Spearman correlation between the score and log TMZ ratio, adjusting for MGMT status. One-sided permutation p (10,000 permutations, shuffle within MGMT group).

**Pass** only if all hold:
- C1 partial rho > 0 with p < 0.05;
- C1 partial rho is higher than every comparator 2-6 (and 7 if run);
- C1 beats at least 95% of the random gene sets.

Anything else is a fail, reported in full.

## Ablation

Comparator 5 is the ablation of the new step (MGMT residual target).

## Power and weak points

- n = 29. About 80% power only for a partial rho near 0.45 at one-sided 5%. A small real effect can be missed. A fail means "not shown", not "no effect".
- In vitro AUC and in vivo survival are different endpoints; failure may reflect the assay gap.
- RNA-seq for PDX was taken from untreated tumours at stated passages; mouse stroma reads may remain.
- Novelty is modest. The new part is training on MGMT-residual in-vitro response and testing transfer to in-vivo survival under an MGMT-adjusted rule. Closest prior art: Zhao et al. 2019 (signatures built on these xenografts), Ju et al. 2025 (HGCC TMZ modelling), PERCEPTION and CaDRReS-Sc. A full prior-art check is needed before any novelty claim.
- One run only. No retuning after the Mayo expression is opened.
