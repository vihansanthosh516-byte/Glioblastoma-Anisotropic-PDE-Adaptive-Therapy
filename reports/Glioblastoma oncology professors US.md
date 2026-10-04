# Mapping computational glioblastoma research across America

US research in computational glioblastoma modeling, adaptive therapy, and multi-omics has consolidated around **three geographic hubs and 40+ active researchers** distributed across distinct but overlapping domains. **Kristin Swanson at Mayo Clinic pioneers mathematical neuro-oncology**, with **Tom Yankeelov's Oden Institute at UT Austin emerging as the largest digital twin hub** and **University of Minnesota's Therapy Modeling & Design Center leading adaptive therapy validation**. The field has matured from foundational reaction-diffusion PDEs (1990s–2000s) toward patient-specific inverse parameter estimation, anisotropic diffusion tensors informed by DTI imaging, and AI-integrated adaptive treatment frameworks funded by ARPA-H ($142 million ADAPT consortium announced 2025). Reaching out to 8–10 "Tier 1" researchers whose work directly overlaps yours—inverse parameter estimation, PDE tumor models, anisotropic diffusion, and adaptive control—will position your glioblastoma research within the established translational pipeline.

## Tier 1 — Perfect fit for your research

These 10 researchers directly address inverse parameter estimation, anisotropic reaction-diffusion PDEs, digital twins, or adaptive therapy—the core technical areas of your work.

**Kristin R. Swanson (Mayo Clinic & Arizona State University)**
- **Institution:** Mayo Clinic, Phoenix; Professor, Arizona State University Department of Mathematical and Statistical Sciences  
- **Department:** Mathematical Neuro-Oncology Lab (Mayo); Neurosurgery Vice Chair of Research  
- **Research Focus:** Pioneer in mathematical tumor forecasting; anisotropic Fisher-Kolmogorov models with DTI-informed diffusion tensors; glioblastoma growth prediction from MRI using reaction-diffusion PDEs; clinical trial validation of model-guided treatment planning  
- **Recent Work:** Developed spatially-dependent anisotropic diffusion model ∂c/∂t = ∇·(D(x)∇c) + ρc accounting for preferential tumor migration along white matter tracts; published foundational work on subthreshold tumor burden estimation ([arXiv 2210.01537](https://arxiv.org/pdf/2210.01537)) showing anisotropic models predict **higher proportion of undetectable disease** than isotropic approaches  
- **Contact:** Phone: 480-342-1914; Address: Support Services Building, Second Floor, 5777 E. Mayo Blvd., SSB 2-700, Phoenix, AZ 85054 — [Mayo Clinic Lab](https://mayo.edu/research/labs/mathematical-neuro-oncology)  
- **Why Perfect Fit:** Direct overlap on anisotropic diffusion tensors parameterized from DTI; inverse parameter estimation from longitudinal MRI; reaction-diffusion equation fundamentals matching your project CLAUDE.md model

---

**Tom Yankeelov (University of Texas at Austin, Oden Institute)**
- **Institution:** University of Texas at Austin, Oden Institute for Computational Engineering and Sciences  
- **Department:** W.A. "Tex" Moncrief, Jr. Chair in Computational Engineering and Sciences; Professor of Biomedical Engineering (Cockrell School); Director, Center for Computational Oncology; joint appointment in Dell Medical School  
- **Research Focus:** Patient-specific digital twins integrating routine clinical MRI into personalized tumor forecasting; TumorTwin Python framework for calibrating biophysical models to individual patient imaging; extended methodology to glioblastoma, breast, prostate, and head-neck cancers  
- **Recent Work:** Co-authored TumorTwin framework ([arXiv 2505.00670v1](https://arxiv.org/abs/2505.00670v1)) for automated parameter estimation and digital twin deployment; received **2026 President's Research Impact Award** for digital twin work; published 2022 Cancer Research study calibrating patient-specific models from **56 triple-negative breast cancer patients' MRI scans** demonstrating clinical translation pipeline  
- **Contact:** Oden Institute (oden.utexas.edu); Center for Computational Oncology director; institutional phone/email available through UT Austin directory  
- **Why Perfect Fit:** TumorTwin directly addresses inverse parameter estimation workflow; digital twin methodology bridges your computational models to clinical validation; collaborative with MD Anderson Cancer Center

---

**David Hormuth II (University of Texas at Austin, Oden Institute)**
- **Institution:** University of Texas at Austin, Oden Institute for Computational Engineering and Sciences  
- **Department:** Faculty researcher, Center for Computational Oncology  
- **Research Focus:** Develops predictive mathematical models of glioblastoma growth and treatment response from subject-specific MRI measurements; biophysical parameter estimation in murine glioblastoma models; co-developer TumorTwin framework  
- **Recent Work:** Lead author/co-author on TumorTwin paper ([arXiv 2505.00670v1](https://arxiv.org/abs/2505.00670v1)); validated biophysical models at preclinical level demonstrating reproducible parameter recovery from imaging  
- **Contact:** Oden Institute (oden.utexas.edu)  
- **Why Perfect Fit:** Direct glioblastoma focus with inverse estimation from non-invasive imaging; hands-on experience with TumorTwin parameter calibration pipeline

---

**Yang Kuang (Arizona State University)**
- **Institution:** Arizona State University, Department of Mathematical and Statistical Sciences  
- **Department:** Professor of Mathematics  
- **Research Focus:** Pioneer in delay differential equations applied to tumor growth and within-host disease dynamics; co-authored "Introduction to Mathematical Oncology" textbook with John D. Nagy  
- **Recent Work:** Established mathematical oncology curriculum at ASU bridging PDE and mechanistic biology; co-founder of interdisciplinary oncology research program  
- **Contact:** ASU Department of Mathematics  
- **Why Perfect Fit:** Expert in mechanistic mathematical biology; collaborates with Swanson at ASU hub; textbook work suggests pedagogical clarity valuable for translating computational models to clinical audiences

---

**Hermann B. Frieboes (University of Louisville)**
- **Institution:** University of Louisville, J.B. Speed School of Engineering  
- **Department:** Professor of Bioengineering  
- **Research Focus:** Physical Oncology framework—studying cancer as multiscale physical system using mathematics, physics, and computation; glioblastoma modeling integrating molecular signaling with spatiotemporal tumor dynamics; anisotropic heterogeneous environments  
- **Recent Work:** Integrated mathematical modeling, computational simulation, and experimental biology across molecular-to-organ scales; published on mechanistic coupling of nutrient/oxygen transport with reaction-diffusion tumor growth; emerging work on phenotype-structured multi-population models  
- **Contact:** University of Louisville Engineering directory  
- **Why Perfect Fit:** Physical oncology perspective complements your PDE approach; heterogeneous anisotropic environment modeling aligns with white-matter-dependent diffusion tensor work

---

**Wolfgang Bangerth (Texas A&M University)**
- **Institution:** Texas A&M University, Department of Mathematics  
- **Department:** Assistant Professor of Mathematics; maintains open-source deal.II finite-element software package widely used for inverse PDE problems  
- **Research Focus:** Numerical solution of inverse problems for PDEs with biomedical imaging applications, particularly optical tomography; PDE-constrained optimization with regularization; sparse localization in tumor growth models  
- **Recent Work:** Developed inverse solvers addressing ill-posed parameter estimation from imaging data; published on tumor growth inverse problems with bounded L-BFGS-B optimization ([arXiv 1907.06564](https://arxiv.org/pdf/1907.06564))  
- **Contact:** Texas A&M Mathematics department  
- **Why Perfect Fit:** Leading expert in inverse PDE problems and optimization algorithms essential for parameter estimation from imaging; deal.II toolkit applicable to glioblastoma PDE solver implementation

---

**Jacob Scott (Cleveland Clinic)**
- **Institution:** Cleveland Clinic, Department of Genomic Sciences & Systems Biology  
- **Department:** Radiation Oncology; Principal Investigator  
- **Research Focus:** Evolutionary game theory applied to adaptive therapy design; dynamic programming for optimal treatment scheduling; biological validation of mathematical models through clinical trials  
- **Recent Work:** Authored survey on adaptive therapy open questions ([PMC10036119](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10036119/)); active clinical trials implementing game-theoretic adaptive therapy for metastatic castration-resistant prostate cancer, ovarian cancer, BRAF-mutant melanoma; validated evolutionary dynamics predict treatment outcomes superior to maximum-tolerated-dose protocols  
- **Contact:** Cleveland Clinic Lerner Research Institute  
- **Why Perfect Fit:** Adaptive therapy framework directly applicable to glioblastoma; evolutionary population dynamics complement your reaction-diffusion PDE models; established clinical trial infrastructure

---

**Robert Gatenby & Renee Brady-Nicholls (Moffitt Cancer Center)**
- **Institution:** Moffitt Cancer Center, H. Lee Moffitt Cancer Center & Research Institute, Tampa, FL  
- **Department:** Mathematical Oncology; Department of Adaptive Therapy  
- **Research Focus:** Population ecology and evolutionary dynamics applied to cancer treatment; adaptive therapy pilot trials showing **>50% progression delay on half standard drug dose** in metastatic prostate cancer; range-bounded dosing optimization  
- **Recent Work:** Brady-Nicholls led clinical data publication on adaptive therapy efficacy ([PMC9657943](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9657943/)); Gatenby directs Integrated Mathematical Oncology (IMO) department; selected as core site in ARPA-H ADAPT consortium ($142 million, 2025–2029) for scaling adaptive approaches  
- **Contact:** Moffitt Cancer Center, IMO Department  
- **Why Perfect Fit:** Direct adaptive therapy clinical validation; evolutionary game theory connects to your robust MPC and RL adaptive steering work; proven clinical efficacy increases translational credibility

---

**Jasmine Foo & Kevin Leder (University of Minnesota)**
- **Institution:** University of Minnesota, School of Mathematics and Department of Industrial & Systems Engineering  
- **Department:** Jasmine Foo—McKnight Distinguished Professor, School of Mathematics; Kevin Leder—Full Professor, Industrial & Systems Engineering  
- **Research Focus:** Co-founded **Therapy Modeling & Design Center (TMDC)** integrating mechanistic mathematical oncology with machine learning; clonal evolutionary dynamics in leukemias and solid tumors; optimal treatment sequencing; personalized medicine  
- **Recent Work:** Foo received **April 2024 Distinguished McKnight University Professorship** for mathematical biology research on cancer evolution and treatment design; Leder develops novel computational approaches pairing cell population evolution models with ML for personalized therapy; both NSF-CAREER-funded; TMDC represents institutional consolidation of mathematical oncology + engineering  
- **Contact:** University of Minnesota School of Mathematics; TMDC official page  
- **Why Perfect Fit:** Mechanistic + ML integration directly matches your RL + robust MPC framework; clonal evolutionary modeling could enhance multi-population glioblastoma simulations; established NSF and McKnight funding track records

---

**Christopher L. Plaisier (Institute for Systems Biology, Seattle)**
- **Institution:** Institute for Systems Biology, Seattle, WA  
- **Department:** Systems biology researcher; SYGNAL pipeline developer  
- **Research Focus:** Causal inference in glioblastoma gene regulatory networks; SYGNAL (SYstems Genetics Network AnaLysis) pipeline integrating correlative, causal, and mechanistic inference to trace mutation → transcription factor dysregulation → gene module changes; multi-omics network pharmacology  
- **Recent Work:** Developed SYGNAL framework ([Cell Systems 2016](https://www.cell.com/cell-systems/fulltext/S2405-4712(16)30189-2)) predicting **112 somatically-mutated driver genes** operating through **74 transcription factors and 37 miRNAs** to dysregulate **237 co-regulated modules** in GBM; 67 TF-miRNA associations previously unknown; framework validated with CRISPR perturbation studies  
- **Contact:** Institute for Systems Biology, Seattle  
- **Why Perfect Fit:** SYGNAL causal network inference bridges your multi-omics Track A work with clinical drug targeting; validated glioblastoma-specific tool; connects single-cell biology to system-level therapeutic strategies

---

## Digital twins and image-based tumor forecasting

**Leading Hub: University of Texas at Austin, Oden Institute**

The **Oden Institute for Computational Engineering and Sciences at UT Austin has consolidated 40+ researchers** around digital twin development, establishing partnerships with **UT Health Dell Medical School, MD Anderson Cancer Center, and Texas Advanced Computing Center (TACC)**. Tom Yankeelov's **Center for Computational Oncology directs** this ecosystem, with co-investigator David Hormuth and co-developer Michael Kapteyn extending digital twin methodology across glioblastoma, breast, prostate, and head-neck cancers. The TumorTwin Python framework ([arXiv 2505.00670v1](https://arxiv.org/abs/2505.00670v1)) represents the most complete toolkit for automated patient-specific calibration from routine clinical MRI, reducing manual parameter inference burden. **Recent FDA grants and 2026 President's Research Impact Award recognition** signal mainstream clinical adoption momentum.

Additional digital twin researchers include **Kyle Lafata, Cristian Badea, and Paul Segars at Duke University's Department of Radiation Oncology**, who secured **$2.9 million NIH R01 funding (2025–2029)** for "Development of a virtual preclinical CT platform for advanced imaging and theranostics in head and neck cancer" ([Duke Radiation Oncology](https://radonc.duke.edu/news/lafata-badea-segars-awarded-2-9-million-r01-advance-digital-twin-models-cancer-research)). **Denis Wirtz at Johns Hopkins University** develops 3D computational tumor models from tissue biopsies using artificial intelligence and multiscale imaging, transforming biopsy data into navigable computational replicas ([Johns Hopkins Hub](https://hub.jhu.edu/2025/04/22/nih-funding-denis-wirtz-tumor-growth/)). **Hassan Fathallah-Shaykh (University of Alabama at Birmingham, MD, PhD Pure Mathematics)** integrates neuro-oncology clinical care with mathematical modeling of brain tumor growth and segmentation using deep learning approaches.

## Inverse parameter estimation and PDE-constrained optimization

**Leading Hub: Texas A&M University + UT Austin**

Parameter estimation from imaging data has emerged as the rate-limiting step in personalizing glioblastoma forecasts. **Wolfgang Bangerth at Texas A&M University** has published extensively on inverse PDE solvers and regularization strategies for ill-posed imaging inverse problems ([arXiv 1907.06564](https://arxiv.org/pdf/1907.06564)), with sparse localization techniques and inexact quasi-Newton methods combined with compressive sampling. His open-source **deal.II finite-element package** is widely adopted for PDE inverse problem implementations.

**Tom Yankeelov and David Hormuth** (UT Austin) shift the inverse estimation problem from theoretical to clinical by developing **dynamic contrast-enhanced MRI (DCE-MRI)–driven parameter recovery** for tumor proliferation rates (ρ) and infiltration diffusion coefficients (D). **Chengyue Wu at UT MD Anderson Cancer Center** integrates biomedical imaging with computational oncology, focusing on imaging-based parameter estimation with joint appointments in Biostatistics and Breast Imaging; received **H.D. Landahl Mathematical Biophysics Award (2023)** from Society for Mathematical Biology. **Chrysanthe Preza at University of Memphis** leads the **Computational Imaging Research Laboratory (CIRL)**, specializing in estimation theory and imaging science. **William Lotter at Dana-Farber Cancer Institute** develops computational methods for precision oncology image analysis.

**Key methodological advances**: Model order reduction (proper orthogonal decomposition, discrete empirical interpolation method) reduces computational burden of repeated PDE solves during optimization; Bayesian parameter fitting outperforms Levenberg-Marquardt in multi-modal MRI datasets; partial volume correction and joint bi-exponential fitting improve parameter identifiability in heterogeneous tumors.

## Adaptive therapy and reinforcement learning for dynamic dosing

**Leading Hub: Moffitt Cancer Center + University of Minnesota + Cleveland Clinic**

Adaptive therapy has transitioned from mathematical modeling to **clinical validation** with **Robert Gatenby and Renee Brady-Nicholls at Moffitt Cancer Center** demonstrating **>50% progression delay on half standard drug dose** in metastatic castration-resistant prostate cancer trials ([PMC9657943](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9657943/)). **Federal investment consolidation** around the **ARPA-H ADAPT program ($142 million, 2025–2029)** signals mainstream recognition; five institutions selected (Moffitt, UC San Diego, NC State, UTA, Cleveland Clinic) to develop "continuously adaptive treatment frameworks."

**University of Minnesota's Therapy Modeling & Design Center (TMDC)** represents the most mature institutional integration of mathematical oncology with machine learning. **Jasmine Foo (McKnight Distinguished Professor, 2024)** and **Kevin Leder (NSF-CAREER-funded)** combine clonal evolutionary dynamics with neural network optimization; **David Odde** (biomedical engineering) provides experimental validation. **Jacob Scott at Cleveland Clinic** applies evolutionary game theory and dynamic programming to adaptive treatment design with ongoing clinical trials in prostate cancer, ovarian cancer, and melanoma.

**International contributions**: **University of Waterloo** researchers **Brydon Eastman, Michelle Przedborski, and Mohammad Kohandel** published peer-reviewed evidence that **deep Q-learning agents generate chemotherapy schedules robust to unknown patient parameter variations** ([PMC8429726](https://pmc.ncbi.nlm.nih.gov/articles/PMC8429726/)), outperforming classical optimal control in sensitivity to model misspecification. **University of Milan-Bicocca** (Fabrizio Angaroni, Marco Antoniotti, Alex Graudenzi) developed **Control Theory for Therapy Design (CT4TD) framework** combining pharmacokinetic/pharmacodynamic models with optimal control for personalized toxicity-minimizing dosing.

## PDE oncology and reaction-diffusion tumor growth

**Geographic Clustering: Phoenix (Mayo + ASU) + Louisville**

**Kristin Swanson at Mayo Clinic and Arizona State University** pioneered the **"Swanson Reaction-Diffusion Model"** for glioblastoma, expressing normalized tumor cell density u(t,x) ∈ [0,1] with spatial diffusion D(x) and growth rate ρ. Her foundational work showed that **anisotropic diffusion modeling predicts higher subthreshold tumor burden** (below MRI detection limit) than isotropic models, with direct implications for radiation therapy planning ([arXiv 2210.01537](https://arxiv.org/pdf/2210.01537)). The model integrates **diffusion tensor imaging (DTI) data to parameterize white-matter-preferred tumor migration** directions.

**Arizona State University's Department of Mathematical and Statistical Sciences** houses a dedicated oncology research hub including **Yang Kuang (delay differential equations specialist)** and **John D. Nagy (evolutionary dynamics; Scottsdale Community College professor**). Together they authored **"Introduction to Mathematical Oncology" textbook** (Routledge), establishing educational infrastructure for the field. **Ashley Harrison**, an ASU PhD student, extends PDE approaches with analysis of mechanistic coupling and heterogeneous environments.

**Hermann Frieboes at University of Louisville** approaches tumor modeling through **Physical Oncology**, treating cancer as a multiscale physical system. His work integrates reaction-diffusion PDEs with **Cahn-Hilliard-type coupled systems** for nutrient/oxygen transport, tissue growth models incorporating advection from cell self-propulsion, and **multi-population models** (proliferating vs. quiescent cells). Recent work on **phenotype-structured reaction-diffusion models** captures heterogeneity at invasive tumor fronts.

**Core mathematical framework** (shared across groups): ∂u/∂t = ∇·(D(x)∇u) + ρu(1−u/K) − γC(t)u (anisotropic Fisher-Kolmogorov); finite element/difference solvers; multigrid techniques for efficiency; parameter estimation via L-BFGS-B bounded optimization and data assimilation from patient MRI.

## Multi-omics, gene regulatory networks, and causal inference

**Leading Hub: Institute for Systems Biology (Seattle) + UCSF + University of Michigan**

**Christopher L. Plaisier at Institute for Systems Biology, Seattle** developed the **SYGNAL (SYstems Genetics Network AnaLysis) pipeline**, a foundational causal inference framework validating that **112 somatically-mutated driver genes** operate through **74 transcription factors and 37 microRNAs** to dysregulate **237 co-regulated gene modules** in glioblastoma ([Cell Systems 2016](https://www.cell.com/cell-systems/fulltext/S2405-4712(16)30189-2)). Notably, **67 TF-miRNA associations were previously unknown**, and the framework includes CRISPR-Cas9 perturbation validation, establishing causal (not merely correlative) gene regulatory networks.

**UCSF Department of Neurological Surgery** conducted extensive single-cell profiling: a **single-cell RNA-seq + exome sequencing study** profiled **37,000+ cells** from **12 glioblastoma patients**, identifying hierarchical glioblastoma stem cell populations and progeny; a follow-up study examined **86 primary-recurrent patient-matched pairs** with **multi-modality single-nucleus RNA, single-cell open chromatin (snATAC-seq), DNA sequencing, and spatial transcriptomic/proteomic assays**, revealing **mesenchymal phenotype enrichment in recurrent tumors** after therapy ([PMC9767870](https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/)).

**University of Michigan** leads **gene therapy-immunotherapy translation**: **Maria G. Castro, R. C. Schneider Collegiate Professor of Neurosurgery**, and **Pedro Lowenstein** (collaborator) are conducting a **Phase 1 first-in-human trial of adenoviral vectors expressing HSV1-TK and Flt3L**, funded by **$4.5 million NIH grant** to understand tumor biology and develop novel treatments. **Nitin Baliga and the Baliga Lab at Institute for Systems Biology** apply systems biology frameworks to glioblastoma's regulatory architecture.

**Recent causal inference advances** (2024–2025): **Cell-type-specific causal inference** identified **fourteen cell-type-specific causal effects** in glioblastomagenesis, with three high-confidence genes (EGFR in astrocytes, CDKN2A in oligodendrocyte progenitor cells, JAK1 in excitatory neurons) ([2025 preprint](https://www.medrxiv.org/content/10.1101/2025.06.28.25330486v2.full)). **Mendelian randomization approaches** combined eQTL/pQTL analyses to identify **2,528 differentially expressed genes** and novel causal plasma proteins (AKR1C4); separate studies identified **GPX7 and CXCL10 as potential causal targets** with evidence across multiple brain regions ([PMC11780873](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11780873/), ([s12885-025-13979-3](https://link.springer.com/article/10.1186/s12885-025-13979-3))).

**Integration frameworks** (2024–2025): Comprehensive multi-omics machine learning frameworks combining **single-cell RNA-seq, snATAC-seq, and bulk RNA-seq** datasets identified core gene modules and therapeutic targets specific to glioblastoma subtypes, with **CEBPG as a key regulatory factor in mesenchymal GBM** ([s41598-025-09742-0](https://www.nature.com/articles/s41598-025-09742-0)). **Agent-based spatial computational models** predict combination chemotherapy, oncolytic virus, and immune checkpoint inhibitor efficacy by simulating individual cells and microenvironmental interactions ([s41540-024-00419-4](https://www.nature.com/articles/s41540-024-00419-4)).

## Contacting professors: Three actionable steps

### Step 1: Tier 1 research overlap briefing (1 week)

Prepare a **2–3 page research brief** positioning your work within each professor's framework:
- **For Swanson**: Anisotropic DTI parameterization advances; novel MRI features for D(x,y,z) tensor recovery
- **For Yankeelov/Hormuth**: TumorTwin integration pathway; inverse estimation validation on retrospective cohort
- **For Scott/Foo/Leder**: Adaptive therapy + robust MPC + RL chronotherapy combination; evolutionary game theory grounding
- **For Plaisier**: Multi-omics integration enabling drug target selection from causal GRNs
- **For Bangerth**: Novel regularization or surrogate model advances in inverse PDE solvers

Include: 1–2 preprints, 3–4 figures showing your anisotropic diffusion PDE results, clinical validation roadmap.

### Step 2: Institutional cluster outreach (2–4 weeks)

**Tier 1 cluster:** Mayo Clinic + Arizona State University (Swanson, Kuang, Nagy, Frieboes).
**Tier 2 cluster:** UT Austin Oden Institute (Yankeelov, Hormuth, Kapteyn) + UT MD Anderson + TACC collaborative partnership.
**Tier 3 cluster:** University of Minnesota (Foo, Leder, Odde) for adaptive therapy + mechanistic ML crossover.

Email sequence:
- **Week 1:** Brief inquiry to Tier 1 lead (Swanson) with 1-page research summary + preprint
- **Week 2:** Follow-up with Oden Institute point-of-contact (Center for Computational Oncology) for digital twin integration discussion
- **Week 3:** Reach out to Moffitt IMO (Gatenby/Brady-Nicholls) with adaptive therapy focus
- **Week 4:** Request informational 30-min call with one Tier 1 researcher; ask for lab collaborator introduction

### Step 3: Conference and summer school presence (3–6 months)

**Established venues** where these researchers convene:
- **Society for Mathematical Biology (SMB) annual meetings** — Swanson, Kuang, Nagy, Frieboes, Scott, Foo, Leder regularly present; Scott chairs adaptive therapy workshop
- **Mathematical Oncology course at Moffitt Cancer Center** — multi-week summer intensive; faculty from Cleveland Clinic, Mayo, UT Austin teach
- **Oden Institute summer workshops** on digital twins and inverse problems — Yankeelov, Hormuth direct; typically June–July
- **SIAM (Society for Industrial and Applied Mathematics) conferences** — inverse problems track (Bangerth, Yankeelov);  computational oncology mini-symposium
- **Glioblastoma Foundation Research Summit** — includes computational oncology sessions; opportunity to network with Swanson, Castro, clinical neuro-oncologists

**Presentation strategy:** Submit contributed talk (5 min) or poster on anisotropic diffusion + inverse parameter estimation + adaptive robust MPC integration; frame as methodological advance enabling clinical glioblastoma forecasting.

## Leading institutions for multi-track support

**Mayo Clinic (Phoenix, AZ)** – **Tier 1**
- **Mathematical Neuro-Oncology Laboratory** directed by Kristin Swanson; longest-established glioblastoma-specific program in US; 15+ years clinical translation experience; co-director Precision Neurotherapeutics Innovation Program
- Address: Support Services Building, Second Floor, 5777 E. Mayo Blvd., SSB 2-700, Phoenix, AZ 85054
- Phone: 480-342-1914
- **Strengths:** Anisotropic diffusion modeling, DTI-informed parameter estimation, clinical trial infrastructure, medical imaging integration

**Arizona State University (Tempe/Phoenix, AZ)** – **Tier 1**
- **Department of Mathematical and Statistical Sciences** + secondary appointment for Kristin Swanson; Yang Kuang (delay differential equations), John D. Nagy (evolutionary dynamics), Ashley Harrison (PDE analysis)
- **Strengths:** Mathematical oncology education, mechanistic PDE theory, delay dynamics in chemotherapy scheduling

**University of Texas at Austin, Oden Institute (Austin, TX)** – **Tier 1**
- **Center for Computational Oncology** directed by Tom Yankeelov; David Hormuth II, Michael Kapteyn; **TumorTwin framework** co-developers; partnerships with UT Health Dell Medical School, MD Anderson, TACC
- **Strengths:** Digital twin automation, inverse parameter estimation from MRI, multi-cancer validation, computational infrastructure (TACC supercomputing access)

**University of Minnesota (Minneapolis, MN)** – **Tier 1 (Adaptive Therapy Focus)**
- **Therapy Modeling & Design Center (TMDC)**: Jasmine Foo (McKnight Distinguished Professor), Kevin Leder (NSF-CAREER), David Odde
- **Strengths:** Mathematical + machine learning integration, clonal evolutionary dynamics, optimal treatment sequencing, NSF funding track record

**Moffitt Cancer Center (Tampa, FL)** – **Tier 1 (Clinical Validation Focus)**
- **Integrated Mathematical Oncology (IMO) Department**; Robert Gatenby, Renee Brady-Nicholls; part of ARPA-H ADAPT consortium ($142M, 2025–2029)
- **Strengths:** Adaptive therapy clinical trials, evolutionary game theory, clinical outcomes data, established translational infrastructure

**Cleveland Clinic (Cleveland, OH)**
- **Department of Genomic Sciences & Systems Biology**: Jacob Scott (MD, PhD)
- **Strengths:** Evolutionary game theory, dynamic programming, active clinical trials (prostate, ovarian, melanoma adaptive therapy)

**Institute for Systems Biology (Seattle, WA)**
- **Christopher L. Plaisier, Nitin Baliga Lab**: SYGNAL causal inference pipeline; glioblastoma gene regulatory networks; systems biology frameworks
- **Strengths:** Causal multi-omics inference, transcriptional regulatory network discovery, drug target prioritization

**Duke University School of Medicine (Durham, NC)**
- **Department of Radiation Oncology**: Kyle Lafata, Cristian Badea, Paul Segars; **$2.9M NIH R01 grant** (2025–2029) for digital twin CT imaging
- **Strengths:** Virtual imaging platforms, preclinical CT digital twins, head-and-neck cancer (transferable to CNS)

**University of Louisville (Louisville, KY)**
- **J.B. Speed School of Engineering, Department of Bioengineering**: Hermann Frieboes
- **Strengths:** Physical oncology framework, anisotropic heterogeneous environments, multi-scale mechanistic coupling (nutrient/oxygen transport + PDE)

**Johns Hopkins University (Baltimore, MD)**
- **Biomedical Engineering**: Denis Wirtz
- **Strengths:** 3D computational tumor models from biopsies, AI-enhanced multiscale imaging, tissue physics

**Texas A&M University (College Station, TX)**
- **Department of Mathematics**: Wolfgang Bangerth; deal.II finite-element software maintainer
- **Strengths:** Inverse PDE solvers, regularization theory, ill-posed problem numerical methods

**University of Alabama at Birmingham (Birmingham, AL)**
- **Department of Neurology and Neuro-Oncology**: Hassan Fathallah-Shaykh (MD, PhD Pure Mathematics)
- **Strengths:** Clinical neuro-oncology + mathematical modeling integration, brain tumor segmentation (deep learning), dynamical systems

**UCSF Department of Neurological Surgery (San Francisco, CA)**
- **Single-cell multi-omics profiling**: 37,000+ cells, primary-recurrent paired specimens, spatial transcriptomics/proteomics
- **Strengths:** Large-scale single-cell cohorts, spatial omics integration, clinical neurosurgery validation

**University of Michigan (Ann Arbor/Ann Arbor, MI)**
- **Department of Neurosurgery, Department of Cell and Developmental Biology**: Maria G. Castro (R. C. Schneider Collegiate Professor), Pedro Lowenstein
- **Strengths:** Gene therapy-immunotherapy combination, Phase 1 clinical trials, $4.5M NIH funding, translational integration

## Conclusion

US computational glioblastoma research has matured from mathematical novelty (1990s) to clinical translation infrastructure (2025) with **three geographic hubs** (Phoenix/Tempe, Austin, Twin Cities/Cleveland) housing researchers at all six technical frontiers your work spans: **digital twins, inverse parameter estimation, anisotropic diffusion tensors, adaptive therapy, reaction-diffusion PDEs, and multi-omics causal networks**. The **8–10 Tier 1 researchers** identified share overlapping journals (Cancer Research, Mathematical Biosciences, SIAM Multiscale Modeling & Simulation), conference venues (SMB, SIAM, Society for Mathematical Biology), and funding sources (NIH NCI, NSF CAREER, ARPA-H ADAPT). **Federal consolidation around the ARPA-H ADAPT program ($142 million, 2025–2029)** and **recent McKnight, NSF, and President's Research Impact Awards (2024–2026)** signal accelerating institutional recognition and career-stage support for computational oncology.

**Reaching out to Swanson, Yankeelov, Scott, Foo, or Plaisier with a 2–3 page research brief positioning your anisotropic inverse-parameter adaptive-control work will position you within the established translational pipeline.** Conference presentations at SMB (especially adaptive therapy workshop chaired by Scott) and summer schools at Moffitt or Oden Institute offer low-friction networking entry points. **Mayo Clinic and UT Austin offer the richest collaborative ecosystems** for your multi-track research; Cleveland Clinic and University of Minnesota offer strongest adaptive therapy mentorship. Approaching these researchers now—before hypothetical FDA regulatory submissions or clinical trials—allows you to influence methodological standards, access historical cohorts for validation, and establish publication partnerships on novel anisotropic diffusion or surrogate modeling advances.
