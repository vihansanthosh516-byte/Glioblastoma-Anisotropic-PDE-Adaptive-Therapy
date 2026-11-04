# Track A Research Plan — Science Fair MD Final

**Scope:** Track A only

**Current status:** A5 is complete and failed its frozen success rule. External confirmation is unavailable. Broad dataset scouting is paused. Track A New is a separate proposed biology study; it has no approved cohort or frozen malignant-cell signature.

Track A contains two distinct research plans. This document is the overview. It keeps their questions, evidence, and claims separate.

## Track A plan structure

1. **A5 — completed frozen study:** [A5 Research Plan](Track%20A5%20Research%20Plan%20Science%20Fair%20MD%20Final.md). Its original research question, methods, result, limits, and materials are recorded there. The A5 protocol, result, and figures remain unchanged.
2. **Track A New — proposed biology study:** [Track A New Research Plan](Track%20A%20New%20Research%20Plan%20Science%20Fair%20MD%20Final.md). It tests a narrower malignant-cell question. It can proceed only after a scope amendment, a qualifying independent validation cohort, and a frozen analysis plan.

## Work completed across Track A

Early exploratory classifier work used the UCSC Cell Browser multiomic-GBM atlas. A historical random-cell split reported 78.7% C-GAT accuracy, while a patient-ID lookup reached 82.8%. A later patient-held-out C-GAT run reached 63.8%, near the 64.6% cells-per-patient baseline. An exploratory expression region classifier reached AUC 0.904 versus 0.713 for cell-type-only prediction, but region was confounded with sequencing library and malignant periphery cells were sparse. These results motivated donor-level testing and composition controls; they are exploratory.

Three bulk regional gene-program tests failed their prespecified adjusted success rules: SenMayo senescence/SASP raw p = 0.0035, hypoxia-adjusted p = 0.37; Hallmark OXPHOS raw p = 0.0024, adjusted p = 0.67; TAM activation raw p = 0.0039, adjusted p = 0.44. These results do not prove the covariates explain the regional effects or that the effects are absent. See the [Track A output record](../archive_old_gbm/vault/TRACK_A_RESULTS.md), [spatial handoff](trackA_spatial_handoff.md), and [composition finding](trackA_composition_finding.md).

A5 then tested its frozen donor-held-out prediction question on Yu et al. 2020 / GEO GSE117891. Its candidate mean AURC was 0.1479 versus 0.1764 for logistic regression, but won in only four of six donors, below the frozen 5/6 rule. Candidate Brier score was worse, and the no-gate ablation had slightly lower AURC. Full methods and caveats are in the A5 plan and [execution audit](trackA_execution_audit_20261031.md).

## Current Track A New feasibility

No screened external validation cohort passes every gate. GSE286419 is the closest label-ready lead, with four retained Core/Margin pairs. GSE318564 has three explicit IDH-wildtype Core/Edge pairs and an unresolved malignant-annotation join. SNUH-GliomaAtlas has six pairs but overlaps prior project use. HRA009021 is a named potential discovery lead: its paper reports 26 IDH-wildtype patients and a 10-person, 30-tissue single-cell subset; an existing metadata-only audit links the donor codes to region labels but not GSA sample accessions to regions. GSA marks the record Open, and NGDC policy permits research use under restrictions; no accession-specific restriction or license is displayed. Exact sample-to-region assignment, malignant-cell label validation, and patient-level independence remain unresolved, so it cannot yet qualify as discovery or validation. Harwood Zenodo 22013127 / concept record 16957940 remains unconfirmed and ineligible; its donor/sample mapping, per-case IDH status, FOV-region map, annotation fields and validation, and patient-level overlap remain unresolved. Its dataset-record license is separate from the article license. The ZIP access deviation is recorded in the [Track A dataset plan](trackA_alternate_dataset_plan_20261101.md).

Keep the six-donor minimum and the independence requirement. Do not present Track A New as external validation until the validation cohort has at least six new, eligible IDH-wildtype donors with both regions mapped to each donor and every other gate passes. Broad scouting remains paused; reopen metadata screening only for a specific named accession that could plausibly pass the gates.

## Shared materials and data record

- **Computer used:** Dell 16 Plus DB16250, Windows 11, 32 GB RAM, CPU-only, 8-thread setting.
- **Python and libraries used across Track A scripts:** Python 3.14.4; NumPy, pandas, scikit-learn, PyTorch, SciPy, PyTorch Geometric, Scanpy, h5py, Matplotlib, Seaborn, openpyxl, and pypdf. A5 package versions appear in its plan.
- **Datasets used in completed Track A work:** Yu et al. 2020 / GEO GSE117891; GBmap Core / CELLxGENE metadata-only annotation snapshot; UCSC Cell Browser multiomic-GBM; Ivy GAP / GSE107559; Angers / E-MTAB-5092; E-MTAB-8743. GBmap expression was not used as A5 features. The Cell Browser atlas is previously exposed and cannot serve as untouched validation.
- Metadata-only candidate screens were not scored datasets. No future validation dataset is approved or selected. The full candidate ledger is in the [Track A dataset plan](trackA_alternate_dataset_plan_20261101.md).
