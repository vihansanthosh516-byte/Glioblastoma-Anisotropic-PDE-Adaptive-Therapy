# Draft manuscript: Scorer choice can reverse patient-level recurrence directions in paired glioblastoma single-nucleus data

Status: draft, 2026-11-16. Based only on the locked result in `paper/recurrence_reliability_result.md`. Not peer reviewed. Every number traces to `output/recurrence_reliability/run_20261116_locked_a1/result.json` or `patient_changes.tsv`.

## Abstract

**Background.** Paired primary–recurrent single-nucleus RNA-seq is increasingly used to describe how an individual glioblastoma changes at recurrence. Patient-level directions are usually reported with nucleus-level sampling uncertainty, but the sensitivity of those directions to equally standard technical choices is rarely measured.
**Methods.** We reanalysed 21 IDH-wildtype patients with an untreated primary and an author-described first-recurrence specimen from GSE174554, restricted to author CNV-labelled malignant nuclei. One fixed program (Neftel MES-like) was scored with two standard procedures: the Tirosh/Neftel control-gene score and the UCell rank score. For each patient we computed the recurrent-minus-primary change and called its direction with a nucleus-bootstrap indeterminate zone. The primary endpoint, opposite-direction discordance between procedures, its tolerance (5%), and all secondary analyses were fixed in a committed protocol before scoring.
**Results.** Directions were opposite in 3 of 21 patients (14.3%; exact 95% CI 3.0–36.3%), so the prespecified reading was inconclusive. Any disagreement occurred in 7 of 21 (33.3%). Split-half checks showed the indeterminate zone was calibrated for nucleus sampling (5.5% and 6.3% false determinate calls), yet all three reversals were determinate under both procedures. Whole-tissue scoring reversed 2–4 patients relative to malignant-only scoring, and patient versus nucleus weighting changed the sign of the cohort summary (−0.097 vs +0.004).
**Conclusions.** In this cohort, scorer choice produced patient-level reversals that nucleus-bootstrap intervals did not anticipate. The sample is too small to bound the reversal rate below a useful tolerance. Patient-level recurrence directions should be reported with scorer sensitivity.

## 1. Introduction

Recurrent glioblastoma shows changes in malignant cell states, including mesenchymal programs (Wang et al. 2022; Spitzer et al. 2025). Recent work reports patient-level recurrence trajectories with patient-level uncertainty (Alexopoulos and Karampelas 2026). Separately, scoring-method benchmarks show that signature scores differ across methods at the cell level (Noureen et al. 2022; Ciernik et al. 2026) and in donor-level pseudobulk designs (PathwayBench 2026), and that composition can dominate mesenchymal scores (Kreis et al. 2024). We did not find a study that measures how often a single patient's paired recurrence direction for one fixed program depends on the scoring procedure. We address that measurement question. We do not propose new recurrence biology.

## 2. Methods

Full protocol: `paper/recurrence_reliability_protocol.md`; amendment 1 (made before scoring): `paper/recurrence_reliability_amendment1_2026-11-16.md`.

- **Data and cohort.** GEO GSE174554 processed per-library matrices and the authors' tumor/normal labels (CONICSmat). Specimen pairing from Supplementary Table 1 cross-checked to GEO by exact specimen ID and stage. Eligibility and batch rules: `paper/recurrence_reliability_eligibility.md`; technical QC: `paper/recurrence_reliability_qc_plan.md`. Flow: 27 candidate pairs → 23 → 21.
- **Program.** Neftel et al. (2019) MESlike1 and MESlike2; per nucleus MES-like = maximum of the two scores.
- **Procedures.** R: log2(CPM/10+1), 30 aggregate-expression bins over the cohort's malignant nuclei, 100 controls per gene, module minus control mean. A: UCell, max rank 1500.
- **Patient change and direction.** delta = recurrent mean − primary mean; SE from 500 nucleus-bootstrap draws per specimen; direction determinate when |delta| > 1.96·SE.
- **Endpoint and inference.** Opposite-direction discordance over all 21 patients; exact Clopper-Pearson 95% CI; reading rule fixed against a 5% tolerance. Secondary analyses S0–S6 descriptive.

## 3. Results

See `paper/recurrence_reliability_result.md` for the full table and per-patient values. Key figure (to be drawn from `patient_changes.tsv`): one row per patient showing R and A changes with their SE bars, coloured by agreement category.

## 4. Discussion

The main observation is a mismatch between two sources of uncertainty. The nucleus bootstrap is well calibrated for sampling nuclei (S6) but blind to the choice between two standard scorers, which reversed three patients' directions. Composition (S1), construct (S4), and weighting (S5) choices produced comparable changes. For the one specimen with a sizeable independent second capture (SF12243), the direction was stable; a smaller second capture (SF10433, 132 nuclei) moved one call to indeterminate.

Limitations: 21 patients; two scorers; one cohort; operational CNV labels; author-described first recurrence; symbol-only features. The interval (3–36%) is compatible with both acceptable and poor reliability. A confirmatory bound below 5% would need roughly 59–72 patients with no reversals.

Next steps: a prespecified multi-scorer extension and an external test in CARE after a metadata and exposure audit.

## Data and code

Code: `src/recurrence_reliability/`; tests: `tests/test_recurrence_reliability.py` (19 synthetic tests). Config: `configs/recurrence_reliability.json`. Raw GEO files are not redistributed; download URLs and SHA-256 are in `output/recurrence_reliability/metadata_manifest.json` and `qc_report.json`.

## References

Wang L. et al. Nat Cancer 2022, doi 10.1038/s43018-022-00475-x. Spitzer A. et al. Nat Genet 2025, doi 10.1038/s41588-025-02168-4. Nomura M. et al. Nat Genet 2025, doi 10.1038/s41588-025-02167-5. Alexopoulos G., Karampelas I. OMICS 2026, doi 10.1177/15578100261499487. Neftel C. et al. Cell 2019, doi 10.1016/j.cell.2019.06.024. Tirosh I. et al. Science 2016, doi 10.1126/science.aad0501. Andreatta M., Carmona S.J. Comput Struct Biotechnol J 2021 (UCell), doi 10.1016/j.csbj.2021.06.043. Verhaak R.G.W. et al. Cancer Cell 2010, doi 10.1016/j.ccr.2009.12.020. Noureen N. et al. eLife 2022, doi 10.7554/eLife.71994. Ciernik L. et al. Genome Res 2026, doi 10.1101/gr.280926.125. Kreis J. et al. Cancer Res Commun 2024, doi 10.1158/2767-9764.CRC-23-0383. PathwayBench 2026, doi 10.21203/rs.3.rs-10271634/v1.
