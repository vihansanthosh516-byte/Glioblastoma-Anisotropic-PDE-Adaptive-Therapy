# New biological tracks for a stronger GBM paper

Research date: 2026-11-16. These are proposals, not frozen protocols, confirmed novel methods, or permission to open reserved data. Sources are papers and dataset descriptions; no expression matrices or new patient outcomes were opened. START_HERE.md remains authoritative. Existing Track A regional-label and TMZ-score failures are not repaired by changing endpoints.

## Which new biological questions are worth adding?

### Takeaway
The best near-term biological addition is a patient-paired recurrence robustness study, or an orthogonal RNA-to-protein spatial validation study. Both need a precise question beyond already-published mesenchymal shifts and hypoxic niches. A drug-perturbation track has a stronger experimental anchor than the failed expression-to-TMZ proxy, but very few independent donors.

### Cited findings
- Wang et al. profiled longitudinal GBM with RNA, chromatin, DNA and spatial assays; mesenchymal recurrence shifts and AP1 regulation are already reported. Their RNA results describe 36 matched primary/first-recurrence patients, whereas the indexed GEO GSE174554 description says 40 primary and 40 matched recurrent specimens. These counts describe different release/analysis scopes and must not be silently equated. [Original paper](https://www.nature.com/articles/s43018-022-00475-x), [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554).
- Greenwald et al. already established multilayered spatial GBM organization and its relation to hypoxia. Spatial RNA is available at GSE237183; processed spatial RNA and CODEX imaging/antibody validation are deposited at Zenodo DOI 10.5281/zenodo.8105466. [Paper and data-availability statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC11088502/), [processed resource](https://doi.org/10.5281/zenodo.8105466).
- Zhao et al. used acute slice cultures from six GBM resections for conserved drug responses, plus an additional resection for a patient-specific screen. Etoposide and panobinostat cell-type effects were already examined; GSE148842 is the accession. Transcriptional perturbation is an experimental treatment readout, but not a clinical response measurement. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8114529/), [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148842).
- Onubogu et al. studied nine matched primary/recurrent GBM cases using multiplexed imaging, and performed single-cell spatial profiling in four tumors. Perivascular collagen/macrophage niches at recurrence are already published. [Original study](https://insight.jci.org/articles/view/179853).

### Inferences
1. **Patient-paired recurrence robustness track (priority 1 for biological coherence).** Hypothesis: the estimated malignant-cell state shift from primary to recurrence is stable under independent state definitions and remains distinguishable from changes in nonmalignant composition. Use GSE174554 as a paired discovery/reanalysis cohort. Freeze two published, non-overlapping state definitions and compare patient-pseudobulk changes with cell-composition-only and simple gene-set baselines. Unit of inference is patient; libraries and cells are nested. Primary endpoint: patient-level paired state-score change, with an interval and a predeclared robustness/replication criterion. Reserve an actually independent second cohort only after origin overlap and paired identifiers are checked. Failure: insufficient paired metadata, state assignment changes the sign materially, or no independent replication. Claim would concern stability of recurrence-state measurement, not discovery of the known MES shift or proof that TMZ caused it. Estimated effort: 3–6 weeks after metadata access; estimate, not measured runtime. Connection to C: motivates a reversible-state model, but does not estimate state-transition rates from two destructive tissue snapshots.
2. **Orthogonal spatial niche validation track (priority 2; potentially strongest evidence addition).** Hypothesis: an RNA-defined niche relation survives measurement with independent protein markers and a spatial null that preserves tissue density. Candidate resource: Greenwald GSE237183/Zenodo 8105466. Freeze one niche, distance scale, protein endpoint, and matched section/sample crosswalk before inspection. Compare spatial graph model with cell abundance only, distance to necrosis alone, and density-preserving label permutations. Primary endpoint: donor-level excess protein-confirmed adjacency or distance effect, not millions of cell-level tests. Hold out complete donors and ideally an independent institution. Failure: no verifiable RNA/CODEX section mapping, effect disappears under spatial null, or cohort-level effect fails. RNA/protein measurements may be adjacent sections, not the same cell; this must be established before analysis. Novelty could be a rigorous cross-assay transportability test, not another generic hypoxia/macrophage story. Estimated effort: 4–8 weeks; existing spatial literature makes novelty narrow.
3. **Experimental perturbation transportability track (priority 3; small but distinct).** Hypothesis: a frozen baseline malignant-state score predicts the direction of a named drug-induced transcriptional response across held-out donors, beyond proliferation alone. GSE148842 supplies treated and control acute slices. Restrict to one drug with sufficient donor coverage after a metadata census; do not assume all seven resections have every drug. Use donor-paired pseudobulk, vehicle controls, proliferation-matched gene sets, and cell-fraction sensitivity analysis. Primary endpoint: held-out donor treatment-versus-control pathway change prediction error or sign agreement. Failure: insufficient matched treated/control donors, no gain over proliferation-only baseline, or result reverses after accounting for captured cell composition. This is not a viability AUC or TMZ clinical-effect model. Repeating published etoposide effects is replication only. The 2025 scFOCAL paper already integrates drug and single-cell expression, limiting broad novelty. [scFOCAL/data description](https://www.nature.com/articles/s41467-025-67783-5). Estimated effort: 2–5 weeks; no wet lab cost, but n≈6 independent donors makes this exploratory unless independently replicated.
4. **Digital pathology plus molecular concordance track (priority 4; genuinely new modality, heavier work).** Hypothesis: interpretable H&E features predict one independently measured molecular phenotype across cohorts, beyond tumor purity and tissue composition. Candidate resources: public CPTAC-GBM pathology images and its multiomics study, with TCGA-GBM as a possible separately locked external cohort only after exposure and overlap audits. Use simple stain-normalized morphology or frozen image embeddings with a linear model before any large neural model. Group all slides/tiles by patient, control tissue source/stain, and report patient-level calibration and cross-cohort performance. Failure: no sample-level slide–RNA join, inadequate independent patients, or no improvement over purity/composition baseline. Broad image-to-expression prediction is established; novelty needs a specific robust endpoint or falsification result. Original CPTAC study covers 99 treatment-naive tumors, not 99 guaranteed matched imaging/RNA cases. [CPTAC primary paper](https://pubmed.ncbi.nlm.nih.gov/33577785/), [TCIA collection](https://www.cancerimagingarchive.net/collection/cptac-gbm/), [IDC collection](https://portal.imaging.datacommons.cancer.gov/collections/cptac_gbm/). Estimated effort: 6–12 weeks, substantial storage and GPU work; sample mapping must precede downloading slides.

### Gaps
- No candidate has a verified eligible cohort in this research session. Publicly listed files do not establish donor independence, identical section origin, rights to redistribute, or successful download.
- GSE174554 direct GEO retrieval failed/was CAPTCHA-blocked; original paper was retrieved. The 36-versus-40 paired count remains unresolved and must be retained explicitly.
- No external paired recurrence cohort was confirmed independent of Wang/GLASS lineage; do not label GSE139533 as a guaranteed external replication set. It is bulk RNA and a recent paper already reanalyzes paired recurrence. [Recent reanalysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC13628663/).
- CPTAC collection full-text verification was incomplete because retrieval failed; collection-level availability was supported by indexed official descriptions. No matched-case number was established.
- Existing project patient exposure is broad. CPTAC and TCGA already participated in integration analyses. They cannot be called untouched overall; a different prespecified endpoint still needs explicit reuse disclosure.

## What would make these additions scientifically strong?

### Takeaway
Choose one new biological endpoint with independent measurement, rather than collecting several additional expression signatures. The publication gain comes from testing whether a proposed link survives a demanding comparison.

### Cited findings
- Published spatial work already links tissue organization, hypoxia and cell states; simply reproducing those associations is not a new biological mechanism. [Greenwald et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11088502/).
- The recurrence atlas already combines longitudinal tissue and spatial assays; generic therapy-induced state switching is prior art. [Wang et al.](https://www.nature.com/articles/s43018-022-00475-x).
- The acute-slice study itself shows cell-type-specific effects and discusses limitations of drug-response models; post-treatment transcription and cell recovery do not establish patient efficacy. [Zhao et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8114529/).

### Inferences
- Keep a single primary endpoint per added track; declare secondary analyses and multiplicity in advance.
- Include composition-only, proliferation-only, or distance-only baselines appropriate to the biology. These baselines directly test likely confounding.
- Separate technical cross-assay concordance from biological causality and clinical utility. Spatial proximity is not verified signaling; paired recurrence is not randomized treatment assignment.
- Use at least two independent cohorts for a strong transportability claim, but assess independence by patient origin, institution and publication lineage, not different accession strings.
- Recommended paper shape: one main question about when GBM measurements remain reliable across patients/modalities, with A and B failure results retained and C explicitly simulation-only. A new biological track can supply an orthogonal validation section if it passes gates.

### Gaps
- The search was focused, not systematic or exhaustive. No claim of first-ever method or guaranteed publishability is justified.
- Original clinical covariates, complete assay matching, public-file licenses, storage requirements and final patient counts need metadata-only checks before any protocol lock.

## Retrieval and scope audit

### Takeaway
Evidence was gathered through Firecrawl searches and full-text retrieval, with a final web fallback. A stalled Firecrawl feedback/retrieval batch was aborted after several minutes; no expression or outcome files were opened.

### Cited findings
- Full article text was retrieved for the recurrence atlas, acute-slice study, Greenwald spatial atlas and Onubogu recurrent niche study using the links above. TCIA/GEO direct page failures remain disclosed rather than treated as verified cohorts.

### Inferences
- The useful next step is a small metadata census of the selected candidate, not immediate model training or expansion to all four candidates.

### Gaps
- Some Firecrawl search-feedback operations may not have completed before interruption. This does not change the cited scientific evidence.
