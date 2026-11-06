# Integrated Track A/B/C novelty update (2026-11-06)

## Status and scope

This is a targeted primary-source update, not an exhaustive systematic review. It does not establish that the project is first, best, or clinically useful. No patient data were read and no outcome benchmark was run for this update.

## Closest prior work

1. **MRI-informed glioma digital twin and treatment optimization.** Chaudhuri et al. developed a predictive digital twin for high-grade glioma using MRI-informed mechanistic modeling, uncertainty quantification, and risk-based radiotherapy-regimen optimization. This covers the broad claim that an MRI-calibrated tumor model can inform personalized therapy under uncertainty.
2. **Longitudinal multiparametric MRI and response forecasting.** A 21-patient high-grade-glioma data-assimilation study used diffusion MRI, T1, T2/FLAIR, and repeated treatment-time observations to predict chemoradiation response and discuss dynamic treatment planning. A 2026 JCO Clinical Cancer Informatics report also describes patient-specific, multiparametric-MRI-based forecasts of tumor-cell density during and after radiotherapy.
3. **Glioblastoma biology plus drug scheduling with reinforcement learning.** Liu et al. (2025) introduced M4RL: a multiscale glioblastoma tumor-cell/TAM model with molecular signaling and cytokine/drug diffusion, a Fokker–Planck/PINN surrogate, and actor-critic reinforcement learning for CSF1R-inhibitor and IGF1R-inhibitor schedules. The paper reports verification using spatial-transcriptomic data. This directly rules out a broad claim that combining GBM biology, mechanistic modeling, and AI drug scheduling is itself new.
4. **GBM world models that connect treatment and future imaging.** Wang et al. (2026) describe Brain-WM, an arXiv preprint that jointly predicts next-step treatment and future MRI from longitudinal GBM data. Its public repository lists LUMIERE, MU-Glioma-Post, and UCSF-ALPTDG as internal cohorts and UCSD-PTGBM as external validation. This makes a broad “integrated GBM model for treatment planning and MRI forecasting” claim especially weak. The paper is a preprint; its reported performance has not been independently reproduced in this project.
5. **Reusable patient-specific digital-twin software.** TumorTwin (2025) presents a public Python framework for composing patient data, models, solvers, and optimization, demonstrated on simulated high-grade-glioma growth and radiotherapy response. A 2026 systematic review also finds many neuro-oncology digital-twin studies, most focused on GBM/HGG. “Building a GBM digital twin” is therefore not a defensible novelty claim by itself.
6. **Robust treatment policies under imperfect biomarkers and model shift.** An October 2026 preprint develops robust contextual reinforcement learning for metastatic prostate-cancer adaptive therapy under measurement error and disease-dynamics shifts. An August 2026 preprint studies adaptive therapy under parameter, structural, and biomarker-measurement uncertainty, including cases where model misspecification changes conclusions. These are different diseases, but they substantially overlap the general robustness idea.

## What this means for this project

- The umbrella title can remain, but the current tracks do not yet form a patient-level multimodal model. Track A and Track B patient data are not linked to the Track C virtual patients. The new signal-to-resistance map is an analyst-chosen simulation assumption, not a relationship measured in those cohorts.
- No evidence currently supports “highest,” “best,” “first,” or a guaranteed positive result. A unit-test pass establishes implementation consistency, not a benefit from the policy.
- The strongest current candidate is a narrower methods question: **Across held-out simulated GBM worlds, can dependence-aware, uncertainty-aware fusion of molecular-like and longitudinal MRI/ADC-like signals improve a treatment policy's worst-endpoint utility over no-signal, single-signal, fixed-fusion, and robust-policy baselines when signal reliability and cross-modal dependence vary?** This remains a hypothesis, not a positive result or verified novelty claim. The new comparison must include Brain-WM conceptually and the robust-biomarker literature directly; a fair code-level comparison may be impossible if their input, treatment, or simulator definitions cannot be matched.
- The proposed method must do more than concatenate signals. A plausible candidate is a belief controller that represents uncertainty in each channel's reliability and their dependence, then discounts or declines a signal when its evidence conflicts across the prespecified observation-model set. This is not implemented or established as novel. It must be tested against the current fixed-likelihood fusion controller, independent-likelihood fusion, each single channel, no-signal control, robust MPC, Q-learning, and fixed schedules.
- The existing track data do not form a linked multimodal patient cohort. Track A can inform only an explicitly assumed molecular-like signal range; Track B can inform only permitted MRI/ADC observation characteristics; Track C supplies the assumed drug-suppression simulator. The molecular-to-resistance relationship remains simulated. The umbrella title can stay, but the protocol must say this clearly.
- To keep the result informative rather than positivity-seeking, evaluate all predeclared signal regimes, including independent, noisy, shuffled, wrong-sign, and mismatched-likelihood controls. Report both progression endpoints, treatment-course use, mean and lower-tail outcomes, and the full regime map. Do not choose a winning regime after seeing evaluation outcomes.

## Required next steps before outcome simulation

1. Freeze the data-generating signal regimes separately from the policy's fixed assumed likelihood; make the no-signal, molecular-only, MRI-only, and joint-signal arms use the same broad prior and patient worlds.
2. Freeze the shared training/evaluation split, controller and Q-learning training information, scan-noise streams, progression endpoints, treatment budget, primary utility, lower-tail metric, multiplicity rule, and seed assignments. Register development seeds before use; do not touch the reserved sealed block.
3. Wire the methods into one harness that evaluates them on the same held-out virtual patients and reports whether signal fusion helps, ties, or harms. Include simple fixed schedules and the pre-existing treatment rules as comparators.
4. Add and pass synthetic leakage, channel-removal, signal-shuffle, wrong-sign, signal-model-mismatch, and identical-observation tests. Then run development outcomes once under the frozen procedure; reserve any confirmation for a separate authorized protocol.
5. Revise the science-fair plan only after the benchmark protocol is reviewed. Describe the project as a synthetic methods study unless a defensible patient-linked biological validation is completed.

## Additional primary-source checks (2026-11-06)

- Brain-WM arXiv record: https://arxiv.org/abs/2603.07562. Its public code repository describes joint treatment recommendation and future MRI generation and names its internal and external datasets: https://github.com/thibault-wch/Brain-GBM-world-model.
- M4RL article and public code: https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/ and https://github.com/SunXQlab/M4RL.
- TumorTwin preprint: https://arxiv.org/abs/2505.00670.
- Robust contextual RL under measurement and dynamics uncertainty: https://arxiv.org/abs/2610.22336.
- Patient-specific HGG MRI data assimilation: https://pmc.ncbi.nlm.nih.gov/articles/PMC12308963/.
- 2026 neuro-oncology digital-twin review: https://pmc.ncbi.nlm.nih.gov/articles/PMC13036663/. This review is field context, not primary evidence of exact methodological overlap.

## Sources checked

- Chaudhuri et al. (2023), [Predictive digital twin for optimizing patient-specific radiotherapy regimens under uncertainty in high-grade gliomas](https://pmc.ncbi.nlm.nih.gov/articles/PMC10598726/).
- [A data assimilation framework for predicting the spatiotemporal response of high-grade gliomas to chemoradiation](https://pmc.ncbi.nlm.nih.gov/articles/PMC12308963/).
- Liu et al. (2025), [Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/), *Science Advances*, 11(32), eadv3316.
- [Patient-Specific Computational Framework to Spatiotemporally Forecast Treatment Response and Identify Biologically Informed Radiotherapy Targets in High-Grade Glioma](https://ascopubs.org/doi/abs/10.1200/CCI-25-00336), *JCO Clinical Cancer Informatics*.
- Xie et al. (2026), [Robust Reinforcement Learning for Dynamic Treatment Regimes under Model and Measurement Uncertainty](https://arxiv.org/abs/2610.22336), preprint posted October 15, 2026.
- Browning et al. (2026), [Adaptive therapy under parametric, structural, and measurement uncertainty](https://arxiv.org/abs/2608.18387), preprint posted August 2026.
