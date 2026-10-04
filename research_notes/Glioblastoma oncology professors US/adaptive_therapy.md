# Adaptive Therapy and Machine Learning for Cancer Treatment

## Who publishes on adaptive therapy, dynamic dosing, or treatment optimization for cancer?

### Takeaway
Adaptive therapy has emerged as a validated approach with multiple research groups publishing peer-reviewed work and clinical trials underway. Key researchers combine mathematical modeling with evolutionary principles to show that dynamic, response-based treatment is superior to fixed maximum-tolerated-dose protocols.

### Cited Findings

**Core Research Community:**
- Jacob Scott (Cleveland Clinic, Radiation Oncology) co-authored foundational work on adaptive drug therapy optimization using dynamic programming and evolutionary game theory, with an ongoing clinical trial framework for metastatic castration-resistant prostate cancer, ovarian cancer, and BRAF-mutant melanoma — [A survey of open questions in adaptive therapy: Bridging mathematics and clinical translation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10036119/)
- Robert Gatenby (Moffitt Cancer Center, Florida) leads mathematical oncology research applying population ecology and evolutionary dynamics to show that adaptive drug treatments based on tumor response are more effective than maximum-tolerated dose approaches — [Moffitt Cancer Center research on evolutionary principles](https://www.moffitt.org/newsroom/news-releases/moffitt-researchers-use-mathematical-modeling-and-evolutionary-principles-to-show-importance-of-basing-treatment-decisions-on-tumor-responses/)
- Renee Brady-Nicholls (Moffitt Cancer Center) published clinical data on range-bounded adaptive therapy in metastatic prostate cancer showing tumor progression delay of over 50% compared to conventional treatment while using approximately half the usual drug dose — [Range-Bounded Adaptive Therapy in Metastatic Prostate Cancer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9657943/)

**Clinical Translation:**
- ARPA-H announced the ADvanced Analysis for Precision cancer Therapy (ADAPT) program in May 2025 with up to $142 million in milestone-based funding to support institutions transforming cancer care into dynamic, continuously adaptive treatment frameworks — [ARPA-H ADAPT program 2025](https://www.discoveriesinhealthpolicy.com/2025/08/arpa-h-update-adaptive-treatment.html)
- Clinical trials are actively recruiting for adaptive therapy in prostate cancer (Moffitt Cancer Center), ovarian cancer, and melanoma with published evidence showing improved outcomes and reduced toxicity — [From food crops to cancer clinics article](https://innovations-report.com/?p=382583)

### Inferences
- Adaptive therapy research is transitioning from mathematical modeling validation to clinical implementation, with multiple ongoing trials suggesting this approach may reshape standard oncology practice
- The combination of evolutionary biology, population dynamics, and mathematical modeling is becoming standard methodology in adaptive therapy research

### Gaps
- Specific names of principal investigators at UC San Diego's adaptive therapy team (mentioned as part of ARPA-H ADAPT program) were not located in search results
- Complete publication lists for Brady-Nicholls and other Moffitt researchers were not retrieved

---

## Which labs apply reinforcement learning or machine learning to cancer treatment?

### Takeaway
Machine learning approaches to cancer treatment have split into two primary methodologies: reinforcement learning for dynamic dosing optimization, and mechanistic machine learning that integrates tumor growth models with neural networks. Waterloo and Bicocca lead RL work while US centers focus on precision medicine applications.

### Cited Findings

**Reinforcement Learning for Drug Dosing:**
- Brydon Eastman, Michelle Przedborski, and Mohammad Kohandel (Department of Applied Mathematics, University of Waterloo) demonstrated that deep double Q-learning agents trained on nominal patient parameters generate chemotherapy dosing schedules more robust to perturbations in unknown patient-specific parameters compared to classical optimal control — [Reinforcement learning derived chemotherapeutic schedules for robust patient-specific therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8429726/)
- MIT researchers developed reinforcement learning algorithms to design novel clinical trial protocols for cancer patient treatment, using agent-based methods to optimize dosing regimens balancing tumor size reduction against neutrophil maintenance — [Reinforcement learning for designing novel clinical trials for treating cancer patients](https://www.media.mit.edu/publications/reinforcement-learning-for-designing-novel-clinical-trials-for-treating-cancer-patients/)

**Mechanistic and Optimal Control Integration:**
- Fabrizio Angaroni, Alex Graudenzi, and Marco Antoniotti (University of Milan-Bicocca, with collaborators at Universitat Ulm, INFN, and Forschungszentrum Jülich) developed the Control Theory for Therapy Design (CT4TD) framework combining pharmacokinetic/pharmacodynamic models with optimal control theory to generate patient-specific dosing schedules minimizing toxicity — [An Optimal Control Framework for the Automated Design of Personalized Cancer Treatments](https://pmc.ncbi.nlm.nih.gov/articles/PMC7270334/)

**Precision Medicine ML:**
- UC San Diego researchers selected as part of ARPA-H ADAPT program to develop machine learning methods for making real-time therapy recommendations based on genetic analysis and tumor biology — [ARPA-H funding for precision cancer therapy](https://idekerlab.ucsd.edu/wp-content/uploads/2025/05/20250523-SanDiegoUnionTribune-ARPAH-Grant.pdf)
- NC State University researchers won TraCS Clinical and Translational Science Pilot funding to harness AI and deep neural networks of genetic and medication data for personalized cancer treatment optimization — [NC State precision medicine initiative 2025](https://news.ncsu.edu/page/477/?cat=health)
- University of Texas at Arlington received $450,000 NSF grant to apply statistical machine learning to predict which patients require additional treatment interventions — [UTA Researcher precision medicine grant](https://biontx.org/community-news/uta-researcher-aims-to-improve-precision-medicine)

### Inferences
- RL approaches are most developed for discrete dosing optimization problems (chemotherapy schedules) while mechanistic ML integration remains more complex for continuous treatment adaptation
- Federal funding (ARPA-H) is consolidating precision oncology ML efforts around integrated "continuously adaptive" frameworks rather than single-modality approaches

### Gaps
- Specific publication records and lab websites for UC San Diego, NC State, and UTA researchers on cancer ML were not retrieved in sufficient detail to name individual PIs
- Direct comparison of RL robustness vs. optimal control methods across the same tumor models was not found in accessible literature

---

## Who works on personalized medicine or precision oncology with computational methods?

### Takeaway
Personalized medicine for cancer has fractured into computational genomics (sequence-based drug response prediction), mechanistic modeling (PK/PD-based dosing), and RL-based dynamic optimization, with the strongest US centers at University of Minnesota, Cleveland Clinic, and Moffitt Cancer Center.

### Cited Findings

**University of Minnesota Precision Therapy Center:**
- Jasmine Foo (School of Mathematics) and Kevin Leder (Industrial & Systems Engineering) co-founded the Therapy Modeling & Design Center (TMDC) at University of Minnesota with David Odde (Biomedical Engineering), an interdisciplinary initiative that "leverages mechanistic mathematical modeling and engineering principles to optimize therapy development" — [Jasmine Foo Distinguished McKnight Award 2024](https://cse.umn.edu/node/164696/)
- Jasmine Foo received the 2014 NSF CAREER Award, 2013 McKnight Land Grant Professorship, and April 2024 Distinguished McKnight University Professorship for mathematical biology research on cancer evolution and treatment design — [Therapy Modeling & Design Center](https://cse.umn.edu/tmdc)
- Kevin Leder develops "novel computational approaches to address heterogeneity in personalized cancer therapy by pairing mathematical models of cell population evolution with machine learning methods," including models of clonal evolutionary dynamics in leukemias to predict optimal combination or sequential therapy — [Kevin Leder NSF-funded cancer research](https://cse.umn.edu/node/84426)

**Cleveland Clinic Evolution-Based Approach:**
- Jacob Scott (Radiation Oncology, Cleveland Clinic) applies "mathematical modeling and the biological and clinical validation of these models" to decompose cancer complexity through evolutionary game theory and dynamic programming, with research translated into ongoing clinical adaptive therapy trials — [Jacob Scott Cleveland Clinic laboratory](https://www.lerner.ccf.org/genomic-sciences-systems-biology/scott/)

**Moffitt Cancer Center:**
- Robert Gatenby and team apply population ecology and evolutionary dynamics, showing that adaptive therapy based on tumor response outperforms fixed protocols; clinical pilot at Moffitt with metastatic castration-resistant prostate cancer patients provided the calibration data for mathematical model validation — [Mathematical Oncology research profile](https://mathematical-oncology.org/people/Jacob-Scott)

**Computational Genomics and Drug Response:**
- Machine learning analysis of multi-omics and drug-screening datasets (e.g., CHIMERA system for breast cancer) combines mechanistic tumor growth modeling with neural networks to predict chemotherapy responses and optimize sequencing decisions — [arxiv pre-print on cancer machine learning](https://www.biorxiv.org/content/10.1101/2020.06.08.140756v1)

**University of California System:**
- UC Berkeley and UCSF have an Assistant Professor in Computational Precision Health developing machine learning methods for personalized medicine translation to clinical care — [Bakar Institute UCSF](https://bakarinstitute.ucsf.sf/)

### Inferences
- The strongest US programs combine mathematical oncology (ordinary differential equations, game theory) with machine learning (neural networks, RL agents), not either approach alone
- University of Minnesota's TMDC represents the most mature institutional integration of mathematical biology with engineering and computational precision medicine
- Federal support (NSF CAREER Awards, McKnight professorships, ARPA-H) is flowing toward groups bridging evolutionary biology and computational methods

### Gaps
- Specific names of UCSF/UC Berkeley Assistant Professor and their publication record were not retrieved
- Comprehensive publication count and citation metrics for named researchers were not gathered
- Contact information and current group composition at each institution was not verified

---

## What are their institutions and recent contributions?

### Takeaway
Leading US institutions are University of Minnesota (mathematical modeling + ML integration), Cleveland Clinic (evolutionary game theory clinical trials), Moffitt Cancer Center (adaptive therapy pilot trials), and emerging precision medicine centers at UC San Diego, NC State, and UT Arlington. Recent contributions (2024–2025) emphasize validated clinical trials and federal consolidation toward "continuously adaptive" treatment frameworks.

### Cited Findings

**University of Minnesota (Minneapolis, MN)**
- Institution: University of Minnesota School of Mathematics, Department of Industrial & Systems Engineering, Department of Biomedical Engineering
- Key PIs: Jasmine Foo (McKnight Distinguished Professor 2024), Kevin Leder (Full Professor, NSF-funded), David Odde (co-founder TMDC)
- Recent contribution: Founded Therapy Modeling & Design Center integrating mathematical oncology with engineering approaches for personalized therapy design — [TMDC official page](https://cse.umn.edu/tmdc)
- Recognition: Jasmine Foo received Distinguished McKnight University Professorship (April 2024) "for outstanding mathematical biology research and contributions to the mathematics community" — [McKnight Professorship award](https://scholarswalk.umn.edu/university-awards/mcknight-distinguished-professors/jasmine-foo)

**Cleveland Clinic (Cleveland, OH)**
- Institution: Cleveland Clinic, Department of Genomic Sciences & Systems Biology
- Key PI: Jacob Scott, MD, PhD (Radiation Oncology, Principal Investigator)
- Recent contribution: Active clinical trials implementing adaptive therapy based on evolutionary game-theoretic models; ongoing validation of dynamic programming approach for optimal treatment scheduling — [A survey of open questions in adaptive therapy](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10036119/)
- Recognition: Scott received NCI funding for evolution-based cancer research and was invited lecturer at major computational oncology conferences

**Moffitt Cancer Center (Tampa, FL)**
- Institution: Moffitt Cancer Center, H. Lee Moffitt Cancer Center & Research Institute
- Key PIs: Robert Gatenby (Mathematical Oncology), Renee Brady-Nicholls (Adaptive Therapy)
- Recent contribution: Completed pilot clinical trial showing adaptive therapy delays prostate cancer progression by >50% on half the standard drug dose — [Moffitt adaptive therapy prostate cancer trial](https://www.urologytimes.com/view/moffitt-study-shows-adaptive-therapy-improves-outcomes-reduces-care-costs-for-prostate-cancer-patients)
- Ongoing: Part of ARPA-H ADAPT consortium ($142 million) to scale adaptive approaches

**UC San Diego (San Diego, CA)**
- Institution: UC San Diego School of Medicine, Precision Medicine Program
- Key PI: Not specifically identified in available sources
- Recent contribution: Selected for ARPA-H ADAPT program to develop machine learning for real-time therapy recommendations based on genetic analysis and tumor biology — [ARPA-H five-organization consortium announcement](https://idekerlab.ucsd.edu/wp-content/uploads/2025/05/20250523-SanDiegoUnionTribune-ARPAH-Grant.pdf)

**NC State University (Raleigh, NC)**
- Institution: NC State University, College of Sciences / College of Engineering
- Key PI: Not specifically identified in available sources
- Recent contribution: Won TraCS Clinical and Translational Science Pilot award to develop deep neural networks integrating genetic and medication data for cancer treatment personalization (2025) — [NC State news announcement](https://news.ncsu.edu/page/477/?cat=health)

**University of Texas at Arlington (Arlington, TX)**
- Institution: University of Texas at Arlington, Department of Industrial & Manufacturing Systems Engineering
- Recent contribution: Received $450,000 NSF grant (NIGMS) to apply statistical machine learning predicting which patients require additional cancer treatment interventions, improving precision oncology outcomes — [UTA precision medicine research](https://biontx.org/community-news/uta-researcher-aims-to-improve-precision-medicine)

**University of Waterloo (Waterloo, Ontario, Canada)**
- Institution: Department of Applied Mathematics
- Key PIs: Brydon Eastman, Michelle Przedborski, Mohammad Kohandel
- Recent contribution: Published peer-reviewed work showing deep Q-learning generates chemotherapy schedules robust to unknown patient parameter variations, advancing RL methods for personalized dosing — [Scientific Reports 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8429726/)

**University of Milan-Bicocca (Italy, multi-institutional)**
- Institution: Department of Informatics, Systems and Communication (lead); collaborators at INFN, Universitat Ulm, Forschungszentrum Jülich
- Key PIs: Fabrizio Angaroni, Alex Graudenzi, Marco Antoniotti (senior author)
- Recent contribution: Developed Control Theory for Therapy Design (CT4TD) framework optimizing patient-specific dosing via pharmacokinetic/pharmacodynamic models and RedCRAB optimization algorithm, reducing toxicity while maximizing efficacy — [CT4TD framework paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7270334/)

### Inferences
- US leadership in adaptive therapy is concentrated at four institutions: University of Minnesota (mathematical theory), Cleveland Clinic (clinical trials), Moffitt (prostate cancer validation), and emerging precision medicine centers (ARPA-H funded)
- Recent federal funding (2025 ARPA-H ADAPT $142 million, 2024 NSF CAREER and McKnight awards) is rewarding groups that combine mathematical modeling with machine learning rather than single methodologies
- International contributions (Waterloo RL, Milan-Bicocca control theory) provide complementary optimization frameworks

### Gaps
- Individual faculty contact information and complete CVs for most researchers were not retrieved
- Graduate student and postdoctoral researcher names at each center were not systematically collected
- Publication counts, h-indices, and detailed citation metrics were not compiled
- Specific grant amounts for UMN, Cleveland Clinic, and Moffitt were not quantified beyond the ARPA-H $142 million consortium figure
