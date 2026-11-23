# Track A crossover data feasibility

## Petti provides a technical development lead

The [Petti methods](https://www.biorxiv.org/content/10.1101/2025.04.04.647096v3.full) state that the third Xenium cohort comprises six sections taken from the same blocks as Visium specimens. This independently supports a same-block cross-platform comparison. It does not resolve the specimen-number discrepancy, establish ROI registration, or add independent patients. Serial sections contain section variation as well as assay variation. The Xenium and Visium platforms also differ in resolution, so a cell-level spatial graph score cannot be subtracted across them without validating a common representation. [Specimen audit](trackA_petti_acquisition_audit.md)

Treat this as a possible technical development panel. Establish an authoritative block-to-section-to-assay crosswalk, same-condition repeats or separate evidence quantifying section effects, a common feature definition, and coordinate registration quality. If only expression summaries are comparable, restrict the assay comparison to those summaries. Do not use a spatial feature unless the challenge covers that feature at its working resolution. Do not count Visium and Xenium from a shared block as independent test patients.

## Independent leads

- [Greenwald et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11088502/) describe near-adjacent CODEX and Visium sections and an explicit registration analysis. This is useful prior art for registration checks. Different molecular modalities and incomplete independent region pairing preclude assuming this is a clean assay crossover or an eligible external cohort.
- [GBM-space](https://www.nature.com/articles/s41586-026-11234-8) includes Visium and Xenium plus other modalities. The paper points to [the portal](https://www.gbmspace.org/) and controlled raw EGA data. A paper-level sample count does not establish independent core/edge patients, matched technical acquisitions, or nonoverlap with previous Track A work. Screen the portal's annotation definitions and join only public metadata until eligibility is verified.
- [Spatial platform reproducibility benchmarking](https://www.nature.com/articles/s41587-025-02811-9) provides a multi-site serial-section design. It can inform a calibration experiment, while its general tissue benchmark does not establish transport to GBM core and infiltrating edge.

No new expression or outcomes were opened. A real test of the candidate requires both acquisition alignment evidence and untouched patient pairs. The [latest simulation](trackA_probe_mismatch_result.md) supplied the alignment flag by design. None of these sources yet establishes that the alignment flag is observable and valid for the actual Track A spatial feature.
