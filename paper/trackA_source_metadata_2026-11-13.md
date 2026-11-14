# Track A metadata audit and access requests

## Harwood CosMx

Use the latest [Zenodo record](https://zenodo.org/records/22013127). It adds Cell_annotation.rds to the archives listed by the first release. The current and earlier archive checksums match for the sample ZIPs. The annotation file remains unopened. It may contain expression as well as labels. The published methods describe two distinct routes to a region label. Fields of view were selected using adjacent pathology, while the differential-expression analysis clustered regions by the inferred malignant-to-resident-cell proportions. For independent validation use only a traceable pathology label. The cohort comes from the same periphery collection described in the earlier Harwood work. Treat the earlier study as a potential source of overlapping donors and verify the relationship by coded IDs before assigning it any role. Sources: [first release](https://zenodo.org/records/16957940), [current release API](https://zenodo.org/api/records/22013127), [published methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12771808/).

Proposed metadata request for the corresponding author:

- Export coded donor ID, slide or archive ID, FOV ID, independently assigned pathology region, and any ambiguity or exclusion flags. Identify the region assignment procedure without using expression or inferred cell proportions.
- Confirm which sequenced cases are primary IDH-wildtype, whether the archive represents the entire paper, and whether coded donors overlap the earlier Notch and synaptic study.
- Provide either a separate cell annotation table or its schema, including join keys, malignant call, state, and validation evidence. Provide aggregate cells per donor, region, and state if available. The requested metadata exclude expression values.

## Gill bulk RNA

The [GEO SOFT metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59612/soft/GSE59612_family.soft.gz) identifies each biopsy's GEO accession, submitted title, and MRI region. The only sample characteristic is tissue type. The identifiers are specimen IDs. Their prefixes cannot establish donor pairing. The paper points to its SI Table S1. PMC returned a verification page for the attachment and PNAS denied source access in this pass. Public expression matrices were not read. The next useful input is the authors' coded participant-by-biopsy map and molecular diagnosis table, or an accessible copy of the published supplement if it contains those fields. Sources: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59612), [source paper and SI link](https://pmc.ncbi.nlm.nih.gov/articles/PMC4151734/#s17).

## Nam direct assay

The pinned [repository README](https://github.com/WangLabHKUST/TMZscript/blob/48a1a13fd946a1ea50013dc6b835ca45d509b1a0/README.md) identifies the featureCounts table as tumor–culture expression pairs used in the published figure and the labeled workbook as training input. The repository names TCGA as its published test dataset, which does not supply direct cultured TMZ response. Checked workbook sheet names only. Ask for the coded pretreatment RNA to measured continuous viability AUC crosswalk, assay definitions and processed data terms. The repository alone does not establish an independent assay validation cohort. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10010007/).

Wu's [specimen audit](trackA_wu_crosswalk_2026-11-13.md) remains applicable. No new candidate expression or response values accessed. No external message sent.
