# Integration protocol: method development on TCGA, one locked confirmation on fresh data (2026-11-12)

Same question: does pre-operative MRI carry tumor biology? Gate 1 (TCGA, hand-made features, ridge) and Gate 1b (Hu 2023 biopsies) failed. The user's view: the method, not the data, is the problem. This protocol lets the method change **without** fishing.

## Code audit before this protocol (development data)

- ID matching verified: in the 54 matched TCGA patients with a BraTS MGMT label, MGMT methylation predicts lower MGMT RNA (AUC 0.741, p 0.0025; shuffled 0.583). So the MRI-RNA pairing is correct.
- RNA state structure is as expected (OPC-NPC rho 0.84; MES vs NPC -0.38).
- No bug found. The failure is a method or signal problem.

## Split

- **Development set (spent, used freely): TCGA anchor, n = 91.** Any features, targets and models may be tried here. Nothing from it is a result.
- **Confirmation set (fresh, untouched until the lock): CPTAC-GBM (23 with BraTS 2021 MRI + masks and open GDC STAR RNA-seq) and IvyGAP (26 with BraTS 2021 MRI + masks and Allen regional RNA-seq).** Only IDs and counts have been seen. RNA values and MRI features of these 49 patients are not looked at until the method is locked.

## Development rules

1. Try methods from the literature report (`paper/integration_methods_literature_2026-11-12.md`, to be written): richer features (texture, location, habitat, deep features), intensity normalization, non-linear models, and dichotomized or continuous targets.
2. Use nested CV on TCGA only. Record every method tried and its CV score in `output/integration_methoddev/log.json` (all tries, not only the best).
3. Pick **one** target (one biology score) and **one** method: the best nested-CV result. Write them, with all code hashes, into a dated lock note before any confirmation data is read.

## Confirmation rule (fixed now)

- Fit the locked method on all 91 TCGA patients. Freeze it.
- Compute the target score in each confirmation cohort with the same gene list. Scores are z-scored within each cohort (platform differences). CPTAC: STAR counts to log2(TPM+1). IvyGAP: mean log2(FPKM+1) over the patient's cellular tumor (CT) samples; if a patient has no CT sample, all samples.
- Predict from MRI in each cohort. Spearman rho per cohort; combined rho by Fisher z weighted by n - 3. One-sided p by 10,000 within-cohort permutations of the combined statistic.
- **Pass:** combined rho >= 0.25 and p < 0.05. One test only; no second target, no second method.

## Weaknesses (known now)

Confirmation n = 49 gives limited power (a true rho of 0.35 has about 70% power at one-sided 0.05). Different RNA platforms in each cohort. IvyGAP RNA is regional and microdissected. Development on 91 patients may overfit even with nested CV.
