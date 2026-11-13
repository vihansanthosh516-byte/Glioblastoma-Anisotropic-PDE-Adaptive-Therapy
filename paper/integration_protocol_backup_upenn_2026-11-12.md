# Integration backup test: does the TCGA-trained MRI biology score predict survival in UPENN? (2026-11-12)

Written before any UPENN survival analysis. Second, independent pre-registered test. It does not replace the CPTAC + IvyGAP confirmation (`paper/integration_protocol_methoddev_2026-11-12.md`); both are reported, whatever they show.

## What has been seen

- UPENN MRI features (images and masks only; no outcome).
- The column names of `data/external/upenn_meta/clin.csv` and its first 2 rows (seen in a format check on 2026-11-12). No survival model, curve or statistic has been computed on UPENN.

## Question

Integration chain: tumor RNA biology (TCGA) -> MRI -> patient outcome (UPENN, no RNA). If the MRI biology score carries real biology, it should carry prognostic information in a new hospital.

## Method (fixed now)

1. Score: the single target and single method in the lock note (`paper/integration_methoddev_lock_2026-11-12.md`), fit on all TCGA development patients, frozen, then applied to UPENN new features (`output/integration_methoddev/features_UPENN-GBM.csv`). No refit on UPENN.
2. Cohort and coding: as `paper/integration_protocol_bamd_note1_2026-11-12.md` (IDH1 Mutated excluded; event = Deceased; time = `Survival_from_surgery_days_UPDATED`, else `Survival_Censor`).
3. Primary test: Cox model, time ~ score_z + age + log(whole-tumor volume). One-sided Wald test on score_z. Direction fixed by biology: if the locked target is MES (or a MES contrast), higher score = higher hazard; for any other target the direction is the sign of that target's TCGA RNA correlation with MES, stated in the lock note.
4. **Pass:** one-sided p < 0.05 in the fixed direction. One test only.
5. Secondary (descriptive, no pass rule): Harrell C-index of age + volume with and without the score, 10 x 5-fold CV.

## Weaknesses (known now)

MES subtype is only weakly prognostic in the literature. Volume adjustment may absorb part of the signal. A pass shows the MRI score adds prognostic information; it does not prove the score measures MES cells.
