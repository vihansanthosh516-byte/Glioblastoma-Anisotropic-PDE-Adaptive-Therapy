# Inverse Parameter Estimation and Medical Imaging Researchers in Oncology

## Who leads labs in medical imaging inverse problems and parameter estimation?

### Takeaway
Leading research in inverse problems for medical imaging is concentrated at institutions like Texas A&M University, University of Texas at Austin, UT MD Anderson Cancer Center, and Dana-Farber Cancer Institute, with researchers developing mathematical frameworks and computational methods for parameter estimation in tumors.

### Cited Findings
- Wolfgang Bangerth, Assistant Professor of Mathematics at Texas A&M University's Department of Mathematics (since 2005), specializes in numerical solution of inverse problems for partial differential equations with applications to biomedical imaging, particularly optical tomography for tumor imaging — [Source](https://iamcs.tamu.edu/wolfgang-bangerth)
- Thomas Yankeelov is Director of the Center for Computational Oncology at the Oden Institute for Computational Engineering & Sciences at University of Texas at Austin, with positions as Professor of Biomedical Engineering in the Cockrell School of Engineering and of Diagnostic Medicine in Dell Medical School — [Source](https://oden.utexas.edu/research/centers-and-groups/center-for-computational-oncology)
- Chengyue Wu is Assistant Professor in the Department of Imaging Physics at University of Texas MD Anderson Cancer Center with joint appointments in Biostatistics and Breast Imaging, awarded the H.D. Landahl Mathematical Biophysics Award in 2023 from the Society of Mathematical Biology — [Source](https://faculty.mdanderson.org/profiles/chengyue_wu.html)
- William (Bill) Lotter is Assistant Professor of Pathology at Dana-Farber Cancer Institute leading an interdisciplinary research group developing computational methods for biomedical image analysis — [Source](https://biophysics.fas.harvard.edu/people/william-bill-lotter)
- Chrysanthe Preza leads the Computational Imaging Research Laboratory (CIRL) at University of Memphis, focusing on imaging science and estimation theory with applications in computational optical sensing, multidimensional light microscopy, and medical imaging — [Source](https://www.memphis.edu/cirl/about/index.php)

### Inferences
- University of Texas at Austin has become a major hub for computational oncology through its Oden Institute and Center for Computational Oncology, bridging mathematical modeling with clinical imaging applications
- Texas A&M University's mathematics department has established expertise in PDE-constrained inverse problems with direct biomedical applications
- Major cancer research centers like UT MD Anderson are embedding computational imaging researchers directly within clinical departments

### Gaps
- Limited information on specific labs at other major institutions (Johns Hopkins, Stanford, MIT, Harvard-affiliated institutions beyond Dana-Farber) working specifically on inverse parameter estimation in tumors
- Institutional funding and lab size information not available from searches

## Which researchers do image-to-parameter estimation for tumors or tissue?

### Takeaway
Multiple researchers across the US are actively developing methods to estimate tumor biology parameters (proliferation rates, diffusion coefficients, necrosis) directly from medical imaging data using PDE-constrained optimization, with particular focus on breast and brain tumors.

### Cited Findings
- Thomas Yankeelov's research develops quantitative medical imaging technologies and integrates imaging data into computational models to predict tumor response to therapy, creating "digital twins" of individual tumors from MRI scans to predict how they respond to different treatments — [Source](https://cco.oden.utexas.edu/?p=395)
- Chengyue Wu's research integrates emerging biomedical imaging techniques with computational modeling to improve cancer diagnosis, prognosis, and treatment, focusing on computational oncology — [Source](https://faculty.mdanderson.org/profiles/chengyue_wu.html)
- Wolfgang Bangerth and colleagues developed inverse solvers with sparse localization for tumor growth models, addressing the ill-posed problem of parameter estimation from imaging data — [Source](https://arxiv.org/pdf/1907.06564)
- Glioblastoma research demonstrates parameter estimation in Cahn-Hilliard type diffuse interface models using model order reduction (proper orthogonal decomposition) and neuroimaging data to define patient-specific diffusion tensors — [Source](https://mox.polimi.it/reports-and-books/publication-results/?id=836)
- DCE-MRI (dynamic contrast-enhanced MRI) driven approaches estimate tumor growth parameters (reaction coefficients for proliferation, diffusion coefficients for infiltration) from time series of imaging data — [Source](https://cris.fau.de/publications/290957046)
- Parameter estimation methods use PDE-constrained optimization with regularization based on biophysical motivation, employing inexact quasi-Newton methods combined with compressive sampling algorithms — [Source](https://arxiv.org/pdf/1408.6221)

### Inferences
- DCE-MRI and multimodal imaging (combining morphology, function, and molecular data) are preferred modalities for capturing tumor dynamics needed for parameter estimation
- Challenges remain in identifiability: growth and diffusion rates cannot always be distinguished from imaging data alone, limiting clinical predictive accuracy
- Model order reduction techniques (POD, DEIM) are becoming standard to reduce computational burden of repeated PDE solves during optimization

### Gaps
- Limited detail on validation against actual patient outcomes for these parameter estimation methods
- Specific quantitative accuracy metrics for parameter recovery in clinical datasets not consistently reported
- Information on comparative performance of different optimization approaches (Levenberg-Marquardt vs. Bayesian vs. deep learning) in tumor parameter estimation

## Who works on automated tumor segmentation and measurement from imaging?

### Takeaway
Automated tumor segmentation is a distributed research area across multiple institutions, with methods ranging from classical machine learning approaches to deep learning, often paired with measurement and parameter estimation for clinical characterization.

### Cited Findings
- Björn Menze (Technical University of Munich) is recognized for tumor segmentation research, leading work on automated segmentation and quantitative parameterization of brain tumors in MRI — [Source](https://www.ias.tum.de/en/ias/menze-bjoern)
- Multi-institutional efforts under The National Alliance for Medical Image Computing (NA-MIC) include projects on tumor modeling with three main aims: automated segmentation of tumors in multi-modal image datasets, image-based parameter estimation in reaction-diffusion models, and processing of magnetic resonance spectroscopic images — [Source](https://www.na-mic.org/wiki/Projects:TumorModeling)
- Segmentation methods employ iterative probabilistic voxel labeling (IPVL) with preprocessing, preliminary segmentation, classification, probability mapping, and final segmentation stages — [Source](https://patents.justia.com/patent/20170147908)
- Hidden Markov random fields (HMRF) and threshold methods applied to T2-weighted MRI images enable segmentation and volume estimation of brain tumors — [Source](https://researchconnect.suny.edu/en/publications/parameter-estimation-and-tissue-segmentation-from-multispectral-m/)
- Joint bi-exponential fitting algorithms provide fast and accurate parameter estimation while correcting for partial volume effects in MRI data — [Source](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10726552)
- Parameter fitting methods using Bayesian approaches yield superior image quality compared to Levenberg-Marquardt methods in Dynamic Contrast Enhanced (DCE) MRI parameter mapping — [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC6010170/)

### Inferences
- Tumor segmentation research has shifted from rule-based methods to statistical and machine learning approaches, with Bayesian methods outperforming traditional optimization approaches
- Multi-modal imaging is becoming standard, with 15+ MRI-based parameters (T1, T2, DWI, DTI, IVIM, DSC) combined with machine learning for robust characterization
- Partial volume correction and joint fitting algorithms are critical for accurate parameter estimation in heterogeneous tumors

### Gaps
- Limited specific information on deep learning approaches and convolutional neural networks for segmentation, despite their prevalence in recent literature
- Comparative benchmarks between different segmentation algorithms across standard datasets (BRATS, TCIA) not detailed in search results
- Clinical validation metrics (sensitivity, specificity) for segmentation algorithms not consistently reported

## What institutions have strong biomedical imaging and computational methods programs?

### Takeaway
Leading US institutions for biomedical imaging and computational methods combine strong engineering programs with clinical research infrastructure, with University of Texas at Austin, Texas A&M, UT MD Anderson, and Wisconsin-Madison being particularly prominent for inverse problems and computational oncology.

### Cited Findings
- University of Texas at Austin hosts the Oden Institute for Computational Engineering & Sciences with the Center for Computational Oncology directed by Thomas Yankeelov, integrating biomedical engineering, mathematics, and clinical oncology — [Source](https://oden.utexas.edu/research/centers-and-groups/center-for-computational-oncology)
- Texas A&M University's Department of Mathematics has Wolfgang Bangerth and maintains the open-source deal.II finite element software package used for inverse problems in computational science — [Source](https://iamcs.tamu.edu/wolfgang-bangerth)
- University of Wisconsin-Madison's Computational Imaging Group (CIG) develops theory and algorithms for imaging including MRI, CT, PET, microscopy, driven by challenges in biomedical imaging — [Source](https://cig.ece.wisc.edu)
- University of Memphis's Computational Imaging Research Laboratory (CIRL) led by Chrysanthe Preza focuses on imaging science, estimation theory, optical sensing, microscopy, and medical imaging — [Source](https://www.memphis.edu/cirl/about/index.php)
- UT MD Anderson Cancer Center hosts imaging physics researchers like Chengyue Wu with joint appointments spanning imaging, biostatistics, and breast imaging, creating integrated clinical-computational research — [Source](https://faculty.mdanderson.org/profiles/chengyue_wu.html)
- Dana-Farber Cancer Institute has established computational imaging methods research through faculty like William Lotter developing precision oncology approaches — [Source](https://biophysics.fas.harvard.edu/people/william-bill-lotter)
- Purdue University offers coursework in computational imaging that formulates signal processing as inverse problems with parameter estimation — [Source](https://engineering.purdue.edu/jump/aac74e)
- University of Central Florida's Computational Imaging Lab develops algorithms and system-level methods for medical imaging, autonomous navigation, and mixed/virtual reality — [Source](https://startupintros.com/orgs/computational-imaging-lab-at-ucf)

### Inferences
- Institutions with strongest programs combine three elements: (1) computational/mathematics departments with inverse problem expertise, (2) biomedical engineering schools with imaging focus, (3) affiliated or embedded clinical cancer research centers
- The University of Texas system (Austin + MD Anderson) appears uniquely positioned through institutional collaboration, with Yankeelov bridging Oden Institute mathematics with cancer center clinical imaging
- Midwest programs (Wisconsin, Illinois, Purdue) also maintain strong computational imaging traditions

### Gaps
- Limited information on industry partnerships (e.g., imaging device manufacturers, pharmaceutical companies) supporting academic inverse problem research
- Specific funding levels and grant sources for inverse problem research in oncology not detailed
- Information on European institutions with comparable programs (mentioned Turin Polytechnic, Technical University of Munich) not fully explored due to search focus on US institutions

---

## Summary of Key Researcher Directory (US Focus)

### Senior Faculty with Inverse Problem Expertise
- **Wolfgang Bangerth** — Texas A&M University, Department of Mathematics; inverse problems, optical tomography, PDE-constrained optimization
- **Thomas Yankeelov** — University of Texas at Austin, Oden Institute and Dell Medical School; computational oncology, digital twins, tumor forecasting
- **Chrysanthe Preza** — University of Memphis, Computational Imaging Research Laboratory; imaging science, estimation theory

### Clinical-Computational Researchers
- **Chengyue Wu** — UT MD Anderson Cancer Center, Department of Imaging Physics; computational oncology, imaging-based parameter estimation
- **William Lotter** — Dana-Farber Cancer Institute; computational image analysis, precision oncology

### Key Institutional Programs
- **Center for Computational Oncology** (UT Austin, Oden Institute)
- **Computational Imaging Research Laboratory** (University of Memphis)
- **Computational Imaging Group** (University of Wisconsin-Madison)
- **National Alliance for Medical Image Computing (NA-MIC)** — Multi-institutional tumor modeling projects

### References for Further Investigation
- Glioblastoma-specific research: Look into collaborations between tumor modeling groups and clinical neuroradiology departments
- Parameter estimation methods: Explore PDE-constrained optimization literature, especially work combining reduced-order models with imaging
- Segmentation and measurement: Track BRATS challenge participants (Brain Tumor Segmentation) for current best practices
