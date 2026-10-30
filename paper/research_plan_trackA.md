# Research Plan — Track A: Core-vs-Edge Gene Programs in Labeled Glioblastoma Tissue

## 1. Scope

- **Data sources:**
  - Ivy Glioblastoma Atlas Project (Ivy GAP), Allen Institute, RNA-seq gene expression matrix release `gene_expression_matrix_2014-11-25`, with sample metadata.
  - Angers multi-zone glioblastoma cohort, ArrayExpress accession E-MTAB-5092 (Agilent arrays, raw files `raw.1` and `raw.2`, with the SDRF sample sheet).
  - MSigDB gene sets: SAUL_SEN_MAYO, HALLMARK_OXIDATIVE_PHOSPHORYLATION, HALLMARK_HYPOXIA, LEIN_NEURON_MARKERS.
- **Question:** Do fixed, externally defined gene programs differ between the tumor edge and the tumor core of the same patient, beyond what known tissue gradients explain?
- **Method:** A within-patient, gene-set score contrast. Each score is corrected by a matched random-gene-set null, then tested with an exact patient-level sign-flip test. The test then adjusts for confound scores.
- **Comparators:** (a) 1,000 random gene sets matched on expression level; (b) the hypoxia program; (c) a neuron-marker program for normal-brain content.
- **Unit of evaluation:** the patient. Each patient gives one paired edge-minus-core value.
- **Evidence for the hypothesis:** a significant patient-level shift that agrees in both cohorts and stays significant after confound adjustment. **Evidence against:** no significant shift, opposite signs between cohorts, or a shift that disappears after adjustment.
- **What the method cannot show:** It cannot show which cell type carries the signal, because bulk tissue mixes tumor and normal cells. It cannot show cause, treatment response, or clinical benefit. It is separate from Tracks B and C. No patient-level link exists between the datasets.

## 2. Rationale

**Background.** Glioblastoma (GBM) spreads into the brain beyond the visible tumor. Surgery removes the core but leaves infiltrating cells at the edge, where most recurrences begin. The Ivy GAP atlas (Puchalski et al., 2018) defined anatomic tumor regions by histology. Aubry et al. (2015) sampled tumor zones and the tumor–brain interface under MRI guidance. Both studies show that gene expression differs by region.

**Specific gap.** Many reported core-vs-edge differences use regions or signatures defined from the same expression data. They also seldom adjust for two strong tissue gradients: low oxygen (hypoxia) in the core, and normal-brain content at the edge. So it stays unclear whether a given program shows tumor biology or only tissue mix. Senescence programs have been linked to GBM progression in mouse and bulk patient data (Salam et al., 2023), but they have not been tested on labeled patient edge vs core with these adjustments. The same is true for oxidative-phosphorylation programs.

**Why it matters.** A careful, confound-adjusted test shows which core-vs-edge claims can be trusted. That helps future work target the residual tumor at the edge. Any benefit to patients would be indirect and in the future. This computational study will not establish clinical benefit.

Glioblastoma research faces the problem that core and edge samples differ in cell mix and oxygen, not only in tumor biology. Existing approaches compare regions but often derive signatures from the same data and do not adjust for these gradients. This project will examine whether fixed external programs differ beyond those gradients through three separate computational studies: Track A (labeled core-vs-edge gene programs), Track B (MRI annotation change between adjacent scans) and Track C (treatment-policy simulation). The work may help clarify which regional expression differences reflect tumor-intrinsic biology, but it will not establish clinical benefit.

## 3. Research Question and Hypothesis

- **Track A question:** In patient-matched, non-expression-labeled GBM edge and core samples, do the SenMayo senescence program and the Hallmark oxidative-phosphorylation program differ beyond expression-matched random gene sets, the hypoxia program, and normal-brain neuron content?
- **Track A hypothesis:** I hypothesize that at least one of these fixed programs will show a patient-level edge-vs-core shift. The shift will be significant by an exact sign-flip test, have the same direction in both cohorts, and stay significant after confound adjustment.
- **What would show it is wrong:** an unadjusted p-value above the alpha level, opposite signs between cohorts, or an adjusted p-value above the alpha level.

## 4. Materials

- **Datasets:**
  - Ivy GAP RNA-seq matrix (`gene_expression_matrix_2014-11-25.zip`, 37,770,060 bytes) and sample/tumor metadata. Source: http://glioblastoma.alleninstitute.org. Terms: Allen Institute Terms of Use, with citation of Puchalski et al. (2018).
  - E-MTAB-5092 raw arrays (`E-MTAB-5092.raw.1.zip`, 193,256,387 bytes; `E-MTAB-5092.raw.2.zip`, 123,796,105 bytes) and `E-MTAB-5092.sdrf.txt`. Source: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5092. Terms: EMBL-EBI Terms of Use, with citation of Aubry et al. (2015).
  - MSigDB gene sets listed in Section 1. Source: https://www.gsea-msigdb.org. Terms: MSigDB license terms.
- **Computer:** Windows 11 Home laptop, Intel Core Ultra 9 288V (8 cores), 31.6 GB RAM. CPU only. No GPU and no cloud service.
- **Software:** Python 3.14.4, NumPy 2.4.4, pandas 2.3.3, Python standard library (`zipfile`, `hashlib`, `itertools`). Git for version control and sealing.
- **Storage:** Local disk only. Raw data stay in a Git-ignored folder and are not shared or uploaded.

## 5. Procedure — Track A

1. **Data selection.**
   - Ivy GAP: include the donors that have both an LE-reference-histology (edge) and a CT-reference-histology (core) sample. Exclude all samples selected by gene expression (CT-control, CTpan, CTmvp, CThbv, CTpnz and similar).
   - Angers: include the GBM patients that have both an interface sample (edge) and a tumor-zone sample (core). The peripheral brain zone is used only in a secondary analysis. Exclude normal-brain controls and necrotic-zone samples.
   - The patient count will be confirmed from the metadata before any analysis. The study needs 10 or more patients.
2. **Gene-list selection and sealing.**
   - Select each fixed external list before any expression data is read.
   - Primary list: SenMayo (125 genes).
   - Confound lists: Hallmark hypoxia (200 genes) and LEIN neuron markers (69 genes).
   - Second list: Hallmark oxidative phosphorylation (200 genes). It will be tested only if the first list fails its success rule, at a stricter alpha of 0.025. No third list will be tested on this cohort.
   - Run a literature search for prior work on each list before sealing.
   - Write each protocol, hash it and the gene lists with SHA-256, and commit to Git before scoring. The analysis code will check the hashes and stop if any file has changed.
3. **Data preparation.**
   - Ivy GAP: log2(FPKM + 1). Map genes by symbol and average duplicate symbols.
   - Angers: use the Agilent processed signal of non-control probes, average replicate spots, transform with log2(x + 1), quantile-normalize across all arrays, and average probes by gene name.
   - Z-score each gene within each cohort. Drop genes with zero variance or missing from a platform, and report the count of genes kept.
   - Average replicate samples within each patient zone.
4. **Data split.** There is no training step. Each patient is one independent unit with one paired value. Samples from the same patient are never treated as independent.
5. **Method.**
   - Score = mean z-score of the list's genes.
   - Per patient, the edge-minus-core score difference is corrected by subtracting the median difference over 1,000 random gene sets. The random sets match the list's size and its distribution over expression-level deciles (seed 20261029).
   - The correction removes shifts shared by all genes of similar expression. Each patient's value is then divided by its cohort's root-mean-square value, so the two platforms are on one scale.
6. **Comparators.** The comparators are the matched random gene sets, the hypoxia score and the neuron-marker score. Each is computed exactly like the main score.
7. **Tuning.** None. All settings, lists, seeds and thresholds are fixed in the sealed protocol.
8. **Evaluation.**
   - Run an exact two-sided sign-flip test over all 2^n sign patterns. The statistic is the mean scaled patient value.
   - Run the adjusted test on the residuals after a no-intercept least-squares fit on the confound scores: hypoxia for list 1, and hypoxia plus neurons for list 2.
9. **Failure handling.**
   - A patient missing either zone is excluded and counted.
   - A missing gene is dropped and counted.
   - If the hashes do not match, the run stops.
   - If the result is positive before adjustment but not after, it is reported as explained by the confound, not as new biology.
10. **Reproducibility.** Record the dataset file names and sizes, the package versions, the seed, the sealed protocol commits, the analysis scripts and the JSON result files.

## 6. Risk Assessment

- **Privacy:** The data are de-identified, public expression data. No images or identifiers are used. No attempt will be made to identify or contact any patient.
- **Data use:** Follow the Allen Institute, EMBL-EBI and MSigDB terms, and cite each source.
- **Security:** Data stay on a local, password-protected computer. Raw files are not published or committed to Git.
- **Cloud:** None will be used.
- **Misinterpretation:** The results are tissue-level associations in bulk samples. They are not medical advice and not treatment evidence.
- **Materials:** No human participants, biological agents or chemicals are used. All work is computational. The fair's review rules for human-derived data will be followed.

## 7. Data Analysis

- **Primary metric:** the mean scaled per-patient edge-minus-core value. Positive means higher at the edge. Negative means higher in the core.
- **Unit:** the patient. Replicates are averaged within each zone.
- **Baseline:** expression-matched random gene sets (built into each patient's value), then the confound scores.
- **Uncertainty:** an exact sign-flip p-value. Per-cohort tests are reported as secondary results.
- **Decision rule (all required):**
  - unadjusted p below alpha (0.05 for list 1; 0.025 for list 2);
  - the same sign in both cohorts;
  - adjusted p below alpha, with the same sign.
- **Exclusions:** The number of excluded patients and dropped genes per cohort is reported.
- **Conclusions:**
  - Positive: evidence of a program shift beyond the named gradients.
  - Positive only before adjustment: the shift is explained by the confound.
  - Negative: no evidence of a shift.
  - A non-significant adjusted result means "not shown", not "absent".

## 8. Bibliography

Aubry, M., de Tayrac, M., Etcheverry, A., Clavreul, A., Saikali, S., Menei, P., & Mosser, J. (2015). From the core to beyond the margin: A genomic picture of glioblastoma intratumor heterogeneity. *Oncotarget, 6*(14), 12094–12109. https://doi.org/10.18632/oncotarget.3297

Lein, E. S., Hawrylycz, M. J., Ao, N., et al. (2007). Genome-wide atlas of gene expression in the adult mouse brain. *Nature, 445*(7124), 168–176. https://doi.org/10.1038/nature05453

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Systems, 1*(6), 417–425. https://doi.org/10.1016/j.cels.2015.12.004

Puchalski, R. B., Shah, N., Miller, J., et al. (2018). An anatomic transcriptional atlas of human glioblastoma. *Science, 360*(6389), 660–663. https://doi.org/10.1126/science.aaf2666

Salam, R., Saliou, A., Bielle, F., et al. (2023). Cellular senescence in malignant cells promotes tumor progression in mouse and patient glioblastoma. *Nature Communications, 14*, 441. https://doi.org/10.1038/s41467-023-36124-9

Saul, D., Kosinsky, R. L., Atkinson, E. J., et al. (2022). A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. *Nature Communications, 13*, 4827. https://doi.org/10.1038/s41467-022-32552-1

Subramanian, A., Tamayo, P., Mootha, V. K., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS, 102*(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102

Allen Institute for Brain Science. (2014). *Ivy Glioblastoma Atlas Project* [Data set]. http://glioblastoma.alleninstitute.org

EMBL-EBI. (2016). *E-MTAB-5092: Multi-zone glioblastoma transcriptome* [Data set]. ArrayExpress. https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5092
