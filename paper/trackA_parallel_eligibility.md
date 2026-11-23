# Track A parallel eligibility findings

## Regional expression lead

[Larroquette et al.](https://doi.org/10.1016/j.canlet.2026.218552) describe matched MRI-navigated core and peripheral bulk RNA profiles. The screen did not identify a public sample-to-patient and region crosswalk or accession. This could be a whole-tissue expression comparison after confirming data access, independence and endpoint compatibility. Bulk profiling cannot reconstruct cell neighborhoods. The authors already explicitly addressed normal-cell admixture, making their filtering strategy relevant prior art.

[Alfardus et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8342598/) describe matched core, rim and margin microRNA. The article's accession placeholder gives no actionable download. This is a different feature modality from the locked mRNA/spatial candidates. No cohort files were opened.

Request the Larroquette coded patient-to-specimen map, independently assigned region and pathology fields, public accession or access procedure, and reuse terms. Ask about any overlap with previously released cohorts. A new expression-only plan must be frozen separately before opening the results.

## Measured TMZ response lead

[Ju et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11700395/) used measured TMZ AUC for treatment-naive cultures in [GSE232173](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232173). The [primary preprocessing script](https://raw.githubusercontent.com/ErasmusMC-Bioinformatics/two-step-TL/main/GSE232173preprocess.R) names a separate response workbook and joins it to expression by culture ID. The public code inventory did not establish access to that workbook. The [GEO series metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232173&targ=self&form=text&view=full) explicitly identifies this expression release as the Ntafoulis study used for prior signature derivation. Exclude this series as independent validation of that signature. The exact overlap of individual culture IDs still matters for any development use.

Request the coded culture-to-origin map, baseline sample IDs, exact TMZ viability outcome and assay provenance, MGMT fields, and a lineage overlap check against Ntafoulis and HGCC. Do not treat a trained weight file, response figure or sample count as a source of joined patient outcomes. Use the existing transfer method as prior art and a comparator where implementable.

## Spatial design additions

Include a segmentation-free comparator such as [Points2Regions](https://doi.org/10.1002/cyto.a.24884) and an orthogonal development annotation audit. Two segmentation methods can agree while sharing spillover. Refer to the [external design](trackA_external_design.md) for exact design power. No external test was run, no messages were sent, and no spent cohort reopened.
