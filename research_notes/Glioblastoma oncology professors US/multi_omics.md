# Multi-Omics Research in Glioblastoma: US Researchers and Computational Approaches

## Who Leads Labs in Glioblastoma Multi-Omics or Single-Cell Analysis?

### Takeaway
Multiple US institutions have established major research programs in glioblastoma multi-omics and single-cell analysis, with UCSF, University of Michigan, and the Institute for Systems Biology leading computational and experimental integration efforts.

### Cited Findings

**UCSF Department of Neurological Surgery:**
- Researchers conducted single-cell RNA-sequencing (scRNA-seq) and matched exome sequencing of human glioblastomas from 12 patients with over 37,000 cells to identify recurrent hierarchies of glioblastoma stem cells (GSCs) and their progeny — [UCSF single-cell atlas](https://www.nature.com/articles/s43018-022-00475-x)
- A related study profiled 86 primary-recurrent patient-matched paired glioblastoma specimens with single-nucleus RNA, single-cell open-chromatin, DNA and spatial transcriptomic/proteomic assays, finding that recurrent glioblastomas are characterized by a shift to a mesenchymal phenotype — [UCSF evolution under therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/)

**Institute for Systems Biology (Seattle):**
- Christopher L Plaisier and colleagues developed the SYGNAL (SYstems Genetics Network AnaLysis) pipeline for dissecting glioblastoma-related causal transcriptional regulatory networks from multi-omics and clinical patient data — [SYGNAL pipeline paper](https://www.cell.com/cell-systems/fulltext/S2405-4712(16)30189-2)
- The Baliga Lab at Institute for Systems Biology works on systems biology approaches to glioblastoma — [Baliga Lab glioblastoma project](https://baliga.systemsbiology.net/projects/glioblastoma/)

**University of Michigan:**
- Maria G Castro, R. C. Schneider Collegiate Professor of Neurosurgery and Professor of Cell and Developmental Biology, leads research on gene therapy-mediated immunotherapy for high-grade gliomas — [U Michigan gene therapy](https://www.michiganmedicine.org/health-lab/turning-tables-glioblastoma)
- Pedro Lowenstein collaborates with Castro on immunotherapy and gene therapy approaches — [Gene therapy combination immunotherapy](https://labblog.uofmhealth.org/lab-report/immunotherapy-gene-therapy-combination-shows-promise-against-glioblastoma)

**Glioblastoma Foundation Genomic Testing & Research Laboratory:**
- Newly established laboratory (opening fall 2024) providing comprehensive genomic testing and integrating somatic mutations and gene expression data for personalized glioblastoma treatment — [Glioblastoma Foundation lab](https://glioblastomafoundation.org/glioblastoma-foundation-announces-the-establishment-of-genomic-testing-research-laboratory/)

### Inferences
- UCSF and University of Michigan represent the two major US centers with active multi-modality single-cell profiling (RNA, chromatin, protein, spatial) in glioblastoma
- Institute for Systems Biology provides the strongest computational framework (SYGNAL) for causal inference in GBM gene regulatory networks
- Recent establishment of dedicated genomic testing laboratories indicates institutional recognition of translational potential

### Gaps
- Specific authors of the UCSF single-cell atlas studies not explicitly named in search results
- Limited information on ongoing multi-omics projects at major comprehensive cancer centers (e.g., MD Anderson, Memorial Sloan Kettering) in search results
- Current funding and team sizes for identified labs not detailed in search results

---

## Which Researchers Work on Gene Regulatory Networks or Causal Inference in Cancer?

### Takeaway
Researchers across US institutions are applying causal inference methods, particularly systems genetics and Mendelian randomization, to identify regulatory networks and drug targets in glioblastoma, with Christopher Plaisier's SYGNAL pipeline being a foundational causal approach.

### Cited Findings

**SYGNAL Pipeline and Mechanistic Causal Networks:**
- Christopher L Plaisier, Sofie O'Brien, Brady Bernard, Sheila Reynolds, and Zac Simon at Institute for Systems Biology developed SYGNAL, which integrates correlative, causal, and mechanistic inference approaches to infer the causal flow from mutations to regulators (transcription factors and microRNAs) to perturbed gene expression patterns — [SYGNAL paper](https://www.cell.com/cell-systems/fulltext/S2405-4712(16)30189-2)
- The SYGNAL network predicted 112 somatically mutated genes or pathways that act through 74 TFs and 37 miRNAs (67 not previously associated with GBM) to dysregulate 237 distinct co-regulated gene modules — [SYGNAL network analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC5001912/)

**Cell-Type-Specific Causal Inference in Multi-Omics:**
- Recent single-cell multi-omic integration studies identified fourteen cell-type-specific causal effects in glioblastomagenesis, including three high-confidence genes (EGFR in astrocytes, CDKN2A in OPCs, and JAK1 in excitatory neurons) — [2025 causal multi-omics](https://www.medrxiv.org/content/10.1101/2025.06.28.25330486v2.full)
- 11 high-confidence and 47 putatively causal genes were prioritized through pharmacogenomic analysis, with most identified as druggable — [Druggable causal genes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11662595/)

**Mendelian Randomization for Causal Drug Target Discovery:**
- Researchers combined eQTL and pQTL analyses using Mendelian randomization to identify causal glioblastoma-related genes, identifying 2,528 differentially expressed genes and novel causal plasma proteins including AKR1C4 — [Brain multi-omic MR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11780873/)
- Mendelian randomization studies identified GPX7 and CXCL10 as potential causal genes for glioblastoma, with thirteen genes showing evidence for causal effects across multiple brain regions — [Mendelian randomization GBM targets](https://link.springer.com/article/10.1186/s12885-025-13979-3)

**Dynamic Gene Regulatory Network Inference:**
- Recent evidence indicates that glioblastoma cell states are shaped by dynamic rewiring of gene regulatory networks in response to microenvironmental cues, with methodologies combining graphical lasso co-expression estimation, Jacobian matrix causal inference, and centrality-based gene importance analysis — [Dynamic network rewiring](https://link.springer.com/article/10.1186/s13040-024-00411-y)

### Inferences
- SYGNAL represents the most comprehensive causal inference framework validated in GBM with confirmed CRISPR-Cas9 perturbation studies
- Mendelian randomization approaches are becoming standard for moving beyond correlative associations to causal inference in multi-omics GBM studies
- Cell-type-specific causal effects require integration of snRNA-seq, snATAC-seq, and bulk genomics data
- Drug target prioritization increasingly relies on causal inference rather than differential expression alone

### Gaps
- Authors of individual Mendelian randomization GBM studies not fully identified in search results
- Validation status of most identified causal genes in wet-lab experiments limited to published CRISPR data
- Specific computational methods compared (e.g., graphical lasso vs. other network inference algorithms) not comprehensively reviewed

---

## Who Integrates Multi-Omic Data with Computational or Experimental Methods?

### Takeaway
US researchers increasingly combine computational multi-omics integration with experimental validation through machine learning, systems biology frameworks, and functional genomics (CRISPR, organoid models), with significant contributions from bioinformatics groups and dedicated computational oncology labs.

### Cited Findings

**Comprehensive Multi-Omics Integration Frameworks:**
- A 2025 comprehensive multi-omics framework combining machine learning with multi-omic data for glioma subtyping was led by corresponding authors Mian Ma and Jiandong Wu, with Yi Ding and Zhaiyue Xu as first authors — [Comprehensive multi-omics ML framework](https://www.nature.com/articles/s41598-025-09742-0)
- Integration of single-cell RNA-seq, snATAC-seq, and bulk RNA-seq datasets has identified core gene modules, candidate therapeutic drugs, and key transcription factors specific to glioblastoma subtypes, particularly CEBPG as a key regulatory factor in mesenchymal GBM — [Mesenchymal GBM transcription factors](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11662595/)

**Functional Precision Medicine with Transcriptomics and Organoid Modeling:**
- Comparative transcriptomics approaches combined with tumor organoid modeling have been developed to identify bespoke treatment strategies for glioblastoma through functional precision medicine pipelines — [Precision medicine pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8699481/)
- Genome-scale metabolic modeling combined with systems biology tools has been applied to examine global transcriptomics data from GBM patients from The Cancer Genome Atlas (TCGA) — [Metabolic modeling GBM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7181968/)

**Machine Learning and High-Throughput Screening Integration:**
- Machine learning predictors trained on high-throughput, image-based screening data for over 12,000 compounds across multiple patient-derived glioblastoma stem cell lines identified HSP90 inhibitor XL-888 with nanomolar potency across all tested cell lines — [ML-driven drug discovery](https://pmc.ncbi.nlm.nih.gov/articles/PMC13181845/)
- The CANDO (Computational Analysis of Novel Drug Opportunities) platform uses multiscale computational analysis, computing interaction scores between drug/compound libraries and proteins to predict new glioma therapies — [CANDO platform](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12139990/)

**Computational Medicinal Chemistry and BBB Optimization:**
- Computational modeling of physicochemical properties combined with synthesis and cellular testing has identified pyridine variants of benzoyl-phenoxy-acetamide with improved blood-brain barrier (BBB) penetration and water solubility for glioblastoma — [Pyridine variants BBB penetration](https://www.nature.com/articles/s41598-023-39236-w)

**Network Pharmacology and Molecular Dynamics Simulation:**
- Comprehensive computational pipelines combining network pharmacology, machine learning, virtual screening, and molecular dynamics simulations have been developed to identify multi-targeting natural compounds against glioblastoma — [Network pharmacology ML GBM](https://www.sciencedirect.com/science/article/pii/S1570180826001065)

**Spatial Transcriptomics and Computational Modeling:**
- Agent-based spatial computational models have been developed to predict combination chemotherapy, oncolytic virus, and immune checkpoint inhibitor effects against glioblastoma by simulating individual cells and local microenvironmental interactions — [Spatial computational modeling](https://www.nature.com/articles/s41540-024-00419-4)

**Bioinformatics Examination and Network Analysis:**
- Bioinformatics examinations of glioblastoma gene expression data have identified panels of therapeutic biomarkers through protein-protein interaction analysis, with 19 hub genes confirmed as potential therapeutic targets including CENPA as a prognostic target — [Bioinformatics biomarker panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC11996113/)

### Inferences
- Multi-omic integration increasingly requires both computational (ML, network inference) and experimental (CRISPR, organoid, BBB penetration testing) validation
- TCGA data combined with patient-derived samples provides the most powerful validation framework
- BBB penetration and CNS pharmacokinetics are becoming integral computational design criteria for GBM drug discovery
- Spatial methods (spatial transcriptomics, agent-based modeling) are addressing the critical gap of tumor microenvironment heterogeneity

### Gaps
- Publication status and current status of most identified candidate drugs from machine learning screens not specified
- Full team composition and institutions for comprehensive 2025 multi-omics framework study not detailed in search results
- Limited information on accessibility and availability of computational tools (CANDO, SYGNAL, others) to broader research community

---

## What Are Their Institutions and Recent Work?

### Takeaway
Leading glioblastoma multi-omics research is concentrated in 8-10 major US institutions, with recent work (2024-2025) focusing on cell-type-specific causal inference, dynamic regulatory networks, and precision medicine pipelines combining computational and experimental validation.

### Cited Findings

**Institute for Systems Biology (Seattle, WA):**
- Christopher L Plaisier (SYGNAL pipeline developer) — [ISB Seattle](https://www.cell.com/cell-systems/fulltext/S2405-4712(16)30189-2)
- Nitin Baliga and Baliga Lab (systems biology approaches) — [Baliga Lab](https://baliga.systemsbiology.net/projects/glioblastoma/)

**University of California, San Francisco (UCSF):**
- Department of Neurological Surgery - single-cell atlas studies with 37,000+ cells profiled — [UCSF Neuro](https://www.nature.com/articles/s43018-022-00475-x)
- Recent work includes spatial transcriptomic/proteomic profiling of primary-recurrent paired specimens with multi-modality data collection — [UCSF spatial omics](https://pmc.ncbi.nlm.nih.gov/articles/PMC13031279/)

**University of Michigan Health:**
- Maria G Castro, R. C. Schneider Collegiate Professor of Neurosurgery, Department of Neurosurgery and Department of Cell and Developmental Biology
- Pedro Lowenstein, collaborator
- Recent work: Phase 1 first-in-human trial of adenoviral vectors expressing HSV1-TK and Flt3L with $4.5M NIH grant for understanding tumor biology and developing new treatments — [U Michigan gene therapy trial](https://www.michiganmedicine.org/health-lab/path-forward-glioblastoma-treatment)

**Michigan State University:**
- 2024: Identified potential glioblastoma treatment using drug-like compound Ogremorphin (OGM) — [MSU 2024 treatment](https://humanmedicine.msu.edu/news/2024-MSU-researchers-find-early-promising-glioblastoma-treatment%20.html)

**University of Toronto:**
- 2024: Led team uncovering new targets for treating glioblastoma through genetic vulnerability screening — [U Toronto genetic targets](https://www.sciencedaily.com/releases/2024/11/241104142208.htm)

**Washington University School of Medicine:**
- 2024: Research on making glioblastoma cells visible to immune cells for immunotherapy — [Wash U immunotherapy](https://www.sciencedirect.com/releases/2024/11/241107160742.htm)

**Collaborative Multi-Institutional Effort:**
- 2024-2025: System biology approach to identify biomarkers involving Prince Sattam Bin Abdulaziz University (Saudi Arabia), Foundation University (Pakistan), Georgia Institute of Technology, and Quaid-i-Azam University (Pakistan) — [Multi-institutional biomarkers](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11150670/)

**University of Alabama:**
- Among leading institutions in ion channel-related glioblastoma research 2005-2024 — [Ion channel research](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12174454/)

**Glioblastoma Foundation:**
- Established Genomic Testing & Research Laboratory (opening fall 2024) for comprehensive genomic testing and personalized treatment identification — [GBM Foundation Lab](https://glioblastomafoundation.org/glioblastoma-foundation-announces-the-establishment-of-genomic-testing-research-laboratory/)

**Sanger Institute Collaboration (International):**
- GBM-space: spatial genomic atlas of glioblastoma project — [Sanger GBM-space](https://www.sanger.ac.uk/collaboration/gbm-space-spatial-genomic-atlas-of-glioblastoma/)

### Inferences
- US glioblastoma multi-omics research is geographically concentrated in the Upper Midwest (Michigan) and West Coast (Washington, California), with emerging centers in Texas and Georgia
- Recent 2024-2025 work shows shift from bulk genomics to single-cell/spatial multi-modality and causal inference approaches
- Major funding support (NIH grants >$4M) is concentrated in translational research groups (Castro lab) and established cancer centers
- International collaboration (Sanger, European institutions) is increasingly common for large-scale cohort studies

### Gaps
- Complete list of all US institutions with glioblastoma multi-omics research not available
- Specific funding amounts and NIH R01 details for most researchers not disclosed in search results
- Postdoc and graduate student team composition for major labs not documented in search results
- Publication pipeline and preprint status of ongoing research often unclear from institutional announcements alone

