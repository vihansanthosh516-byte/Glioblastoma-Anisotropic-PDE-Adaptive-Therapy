# Research Plan — Track A: Core-vs-Edge Gene Programs in Labeled Glioblastoma Tissue

## Rationale

**Background.** Glioblastoma (GBM) spreads into the brain beyond the visible tumor. Surgery removes the core but leaves infiltrating cells at the edge, where most recurrences begin. The Ivy GAP atlas (Puchalski et al., 2018) defined anatomic tumor regions by histology. Aubry et al. (2015) sampled tumor zones and the tumor–brain interface under MRI guidance. Both studies show that gene expression differs by region.

**Specific gap.** Many reported core-vs-edge differences use regions or signatures defined from the same expression data. They also seldom adjust for two strong tissue gradients: low oxygen (hypoxia) in the core, and normal-brain content at the edge. So it stays unclear whether a given program shows tumor biology or only tissue mix. Senescence programs have been linked to GBM progression in mouse and bulk patient data (Salam et al., 2023), but they have not been tested on labeled patient edge vs core with these adjustments. The same is true for oxidative-phosphorylation programs.

**Why it matters.** A careful, confound-adjusted test shows which core-vs-edge claims can be trusted. That helps future work target the residual tumor at the edge. Any benefit to patients would be indirect and in the future. This computational study will not establish clinical benefit.

Glioblastoma research faces the problem that core and edge samples differ in cell mix and oxygen, not only in tumor biology. Existing approaches compare regions but often derive signatures from the same data and do not adjust for these gradients. This project will examine whether fixed external programs differ beyond those gradients through three separate computational studies: Track A (labeled core-vs-edge gene programs), Track B (MRI annotation change between adjacent scans) and Track C (treatment-policy simulation). The work may help clarify which regional expression differences reflect tumor-intrinsic biology, but it will not establish clinical benefit.

## Research Question/Hypothesis

- **Track A question:** In patient-matched, non-expression-labeled GBM edge and core samples, do the SenMayo senescence program and the Hallmark oxidative-phosphorylation program differ beyond expression-matched random gene sets, the hypoxia program, and normal-brain neuron content?
- **Track A hypothesis:** I hypothesize that at least one of these fixed programs will show a patient-level edge-vs-core shift. The shift will be significant by an exact sign-flip test, have the same direction in both cohorts, and stay significant after confound adjustment.
- **What would show it is wrong:** an unadjusted p-value above the alpha level, opposite signs between cohorts, or an adjusted p-value above the alpha level.
- **Relation to my reasons for choosing the project:** Surgery leaves the tumor edge behind, so I want to know which edge-vs-core differences are real tumor biology and which only reflect oxygen level or normal brain. This question tests that directly, with fixed gene lists and labels that do not come from gene expression.

## Materials

- **Datasets:**
  - Ivy GAP RNA-seq matrix (`gene_expression_matrix_2014-11-25.zip`, 37,770,060 bytes) and sample/tumor metadata. Source: http://glioblastoma.alleninstitute.org. Terms: Allen Institute Terms of Use, with citation of Puchalski et al. (2018).
  - E-MTAB-5092 raw arrays (`E-MTAB-5092.raw.1.zip`, 193,256,387 bytes; `E-MTAB-5092.raw.2.zip`, 123,796,105 bytes) and `E-MTAB-5092.sdrf.txt`. Source: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5092. Terms: EMBL-EBI Terms of Use, with citation of Aubry et al. (2015).
  - MSigDB gene sets listed in Section 1. Source: https://www.gsea-msigdb.org. Terms: MSigDB license terms.
- **Computer:** Windows 11 Home laptop, Intel Core Ultra 9 288V (8 cores), 31.6 GB RAM. CPU only. No GPU and no cloud service.
- **Software:** Python 3.14.4, NumPy 2.4.4, pandas 2.3.3, Python standard library (`zipfile`, `hashlib`, `itertools`). Git for version control and sealing.
- **Storage:** Local disk only. Raw data stay in a Git-ignored folder and are not shared or uploaded.
- **Chemicals, biological agents, bacteria:** none. This is a computational project only.

## Procedures

All steps are done by me on my own computer. The data are public datasets made by others (sources in Step 3); I do not collect any new data.

### Study design

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



Each step is written so that another person can run the whole track from an empty folder. Only Python, NumPy and pandas are needed.

### Step 1. Set up
1. Make a project folder with the subfolders `data/external/ivygap/matrix`, `data/external/angers`, `paper/trackA_senmayo`, `paper/trackA_oxphos`, `src` and `output`.
2. Make `data/` Git-ignored. Start a Git repository so that every sealed file gets a dated commit.

### Step 2. Choose and seal the gene lists (before any expression data is downloaded)
1. Run a prior-art search for each candidate list. Search for papers that test the list in labeled GBM edge vs core. Record the queries, the closest papers and what is different.
2. Copy each list, one symbol per line, from its MSigDB gene-set page:
   - `paper/trackA_senmayo/senmayo_genes.txt` — SAUL_SEN_MAYO, 125 symbols.
   - `paper/trackA_senmayo/hallmark_hypoxia_genes.txt` — HALLMARK_HYPOXIA, 200 symbols.
   - `paper/trackA_oxphos/oxphos_genes.txt` — HALLMARK_OXIDATIVE_PHOSPHORYLATION, 200 symbols.
   - `paper/trackA_oxphos/neuron_genes.txt` — LEIN_NEURON_MARKERS, 69 symbols.
   - `paper/trackA_oxphos/hallmark_hypoxia_genes.txt` — the same hypoxia file.
3. Write `PROTOCOL.md` in each folder. It holds the question, the cohort, the method, the success rule and the known risks, exactly as in Steps 4–9 below.
4. In each folder, write `SHA256SUMS` with the SHA-256 hash of every file. Compute each hash with line endings normalized to LF.
5. Commit and push the list-1 folder (`trackA_senmayo`). The list-2 folder (`trackA_oxphos`) is sealed and committed only if list 1 fails its success rule, and before any list-2 score is computed.

### Step 3. Download the data
| File | Source | Size (bytes) |
|---|---|---|
| `gene_expression_matrix_2014-11-25.zip` | http://glioblastoma.alleninstitute.org/api/v2/well_known_file_download/305873915 | 37,770,060 |
| `rna_seq_samples_details.csv`, `tumor_details.csv` | Ivy GAP download page, http://glioblastoma.alleninstitute.org | small |
| `E-MTAB-5092.raw.1.zip` | https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-5092/ | 193,256,387 |
| `E-MTAB-5092.raw.2.zip` | same folder | 123,796,105 |
| `E-MTAB-5092.sdrf.txt` | https://www.ebi.ac.uk/biostudies/files/E-MTAB-5092/E-MTAB-5092.sdrf.txt | 16,025 |

1. Put the Ivy GAP files in `data/external/ivygap/` and unzip the matrix into `matrix/`. This gives `fpkm_table.csv`, `rows-genes.csv` and `columns-samples.csv`.
2. Put the Angers files in `data/external/angers/`. Leave the zips packed.
3. Check that each file size matches the table. If a size is wrong, download the file again.

### Step 4. Select the patients
1. **Ivy GAP.** In `rna_seq_samples_details.csv`, keep the rows where `structure_acronym` is `LE-reference-histology` (edge) or `CT-reference-histology` (core). Keep the donors (`donor_id`) that have both. The expected donors are 10865, 10926, 12108, 12109, 12110, 12112, 14734 and 292023109 (8). The `sample_id` column matches the column names in `fpkm_table.csv`. Do not use any other structure, such as CT-control-*, CTpan-*, CTmvp-*, CThbv-* or CTpnz-*.
2. **Angers.** In the SDRF, the patient is column 5 (`Characteristics[individual]`) and the zone is column 6 (`Characteristics[sampling site]`). Keep the patients with both `interface` (edge) and `tumor zone` (core). The expected patients are FT02, FT04, FT05, FT06, FT08, FT09 and FT13 (7). `Array Data File` gives each array's file name inside `raw.1` or `raw.2`. Use `peripheral brain zone` only for the secondary test. Exclude the normal controls (N1–N4) and the `necrotic zone` samples.
3. Confirm that the total is 15 patients, which is 10 or more. Record any patient that is excluded and why.

### Step 5. Prepare the expression data
1. **Ivy GAP.** Take the selected columns of `fpkm_table.csv`. Apply log2(FPKM + 1). Map `gene_id` to `gene_symbol` with `rows-genes.csv`. Average rows that share a symbol.
2. **Angers.** For each array file:
   - Read the table that starts at the line beginning with `FEATURES`.
   - Keep the rows with `ControlType` = 0.
   - Take `gProcessedSignal` and average the replicate spots of each `ProbeName`.
3. Put all arrays in the SDRF side by side. Drop probes that are missing on any array. Clip negative values to 0 and apply log2(x + 1). Quantile-normalize across all arrays. Map each probe to `GeneName` and average probes per gene. Then keep only the selected samples.
4. **In each cohort separately:** drop genes with zero variance, then z-score each gene across that cohort's selected samples.
5. If a patient has more than one sample in a zone, average the z-scores within that zone.

### Step 6. Score each list per patient
1. For each patient, compute the per-gene difference d(g) = edge z-score minus core z-score.
2. **List score difference D** = the mean of d(g) over the list genes found in that cohort. Record how many list genes were found.
3. **Matched null.** Rank all genes by mean log expression in the cohort and cut them into 10 equal bins (deciles). Draw 1,000 random gene sets with a random generator seeded with 20261029. Each set has the same number of genes per decile as the list and excludes the list genes. Compute D for each set.
4. **Per-patient value:** E = D minus the median of the 1,000 null D values.
5. Compute E the same way for the hypoxia list (H) and, for list 2, for the neuron list (N).

### Step 7. Run the main test
1. Within each cohort, divide each patient's E by the root-mean-square of E in that cohort. Pool the 15 scaled values.
2. Statistic = the mean of the pooled values.
3. **Exact two-sided sign-flip test:** compute the statistic for all 2^15 = 32,768 patterns of +1/−1 signs. p = the fraction of patterns whose absolute mean is at least the observed absolute mean.
4. Record the direction (edge > core if the mean is positive) and the mean E of each cohort.

### Step 8. Run the confound-adjusted test
1. **List 1:** fit E = b·H across the 15 patients by least squares with no intercept. Take the residuals.
2. **List 2:** fit E = b1·H + b2·N the same way. Take the residuals.
3. Scale the residuals within each cohort (Step 7.1) and run the same exact sign-flip test.

### Step 9. Apply the success rule
1. Success needs all three conditions:
   - the main p is below alpha (0.05 for list 1, 0.025 for list 2);
   - the Ivy GAP mean E and the Angers mean E have the same sign;
   - the adjusted p is below alpha, with the same sign as the main test.
2. If the main test passes but the adjusted test fails, record "explained by the confound".
3. If list 1 does not succeed, seal list 2 (Step 2.5) and repeat Steps 6–9 with it, using the same data and code. After list 2, test no more lists on this cohort.

### Step 10. Run the secondary tests (cannot change the main conclusion)
1. Run the test for each cohort alone (Step 7 on that cohort's values only).
2. Run the Angers test with the peripheral brain zone in place of the interface, for patients with both zones (list 1 only).

### Step 11. Failure handling
- **Hash mismatch** at the start of a run: stop. Do not change the sealed files.
- **Missing zone** for a patient: exclude the patient and count the exclusion.
- **Gene missing** from a platform: drop the gene and report the count kept per cohort.
- **Download size mismatch:** download the file again before any analysis.
- **Inconclusive result** (cohorts disagree, or p near alpha): report it as not supported. Do not re-run with new settings.

### Step 12. Record for reproducibility
1. Write one JSON result file per list: `output/trackA_senmayo/result.json` and `output/trackA_oxphos/result.json`. Each file holds the patient IDs, the genes found, every per-patient E, H and N value, all p-values, the coefficients, the seed and the protocol commit.
2. Commit the analysis scripts and the result files. Record the Python, NumPy and pandas versions.
3. Before reporting, check that each result file exists, is not empty and lists 15 patients.

## Risk Assessment

- **Privacy:** The data are de-identified, public expression data. No images or identifiers are used. No attempt will be made to identify or contact any patient.
- **Data use:** Follow the Allen Institute, EMBL-EBI and MSigDB terms, and cite each source.
- **Security:** Data stay on a local, password-protected computer. Raw files are not published or committed to Git.
- **Cloud:** None will be used.
- **Misinterpretation:** The results are tissue-level associations in bulk samples. They are not medical advice and not treatment evidence.
- **Materials:** No human participants, biological agents or chemicals are used. All work is computational. The fair's review rules for human-derived data will be followed.

## Data Analysis

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
- **Formulas:**
  - Gene z-score: z = (x − mean) / SD, within each cohort.
  - Per-patient value: E = D − median(D_null), where D is the mean of (edge z − core z) over the list genes.
  - Scaled value: E / sqrt(mean(E²)) within each cohort.
  - Sign-flip p-value: p = (number of the 2^n sign patterns with |mean| ≥ observed |mean|) / 2^n.
  - Adjustment: residuals of E = b·H (list 1) or E = b1·H + b2·N (list 2), least squares with no intercept.
- **Tables:** one table per list, with the patient count, the genes found per cohort, the main p, the adjusted p, the direction and the per-cohort means.
- **Graphs:** a paired dot plot of E per patient, colored by cohort; a scatter of E against the hypoxia (and neuron) score with the fitted line; a histogram of the 32,768 sign-flip statistics with the observed value marked.

## Bibliography

Aubry, M., de Tayrac, M., Etcheverry, A., Clavreul, A., Saikali, S., Menei, P., & Mosser, J. (2015). From the core to beyond the margin: A genomic picture of glioblastoma intratumor heterogeneity. *Oncotarget, 6*(14), 12094–12109. https://doi.org/10.18632/oncotarget.3297

Lein, E. S., Hawrylycz, M. J., Ao, N., et al. (2007). Genome-wide atlas of gene expression in the adult mouse brain. *Nature, 445*(7124), 168–176. https://doi.org/10.1038/nature05453

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Systems, 1*(6), 417–425. https://doi.org/10.1016/j.cels.2015.12.004

Puchalski, R. B., Shah, N., Miller, J., et al. (2018). An anatomic transcriptional atlas of human glioblastoma. *Science, 360*(6389), 660–663. https://doi.org/10.1126/science.aaf2666

Salam, R., Saliou, A., Bielle, F., et al. (2023). Cellular senescence in malignant cells promotes tumor progression in mouse and patient glioblastoma. *Nature Communications, 14*, 441. https://doi.org/10.1038/s41467-023-36124-9

Saul, D., Kosinsky, R. L., Atkinson, E. J., et al. (2022). A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. *Nature Communications, 13*, 4827. https://doi.org/10.1038/s41467-022-32552-1

Subramanian, A., Tamayo, P., Mootha, V. K., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS, 102*(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102

Allen Institute for Brain Science. (2014). *Ivy Glioblastoma Atlas Project* [Data set]. http://glioblastoma.alleninstitute.org

EMBL-EBI. (2016). *E-MTAB-5092: Multi-zone glioblastoma transcriptome* [Data set]. ArrayExpress. https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5092

## Additional Project Information
- **Human participants:** None. The project is a secondary analysis of publicly available, de-identified expression data (Ivy GAP, ArrayExpress E-MTAB-5092). I will not interact with any person, and I have no access to identifiers. Such projects are usually exempt from IRB review under ISEF rules. I will confirm this with my fair's Scientific Review Committee and attach any form it requires.
- **Vertebrate animals, potentially hazardous biological agents, hazardous materials/activities/devices:** None.
- **Continuation project:** [Answer yes or no truthfully. If yes, upload the earlier abstract and research plan.]
