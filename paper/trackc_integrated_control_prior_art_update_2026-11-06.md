# Track C integrated prior-art and method update (2026-11-06)

## Status

The umbrella title remains **“Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization.”** A targeted primary-source review finds that this broad combination is established. The exact dependence-aware versus matched conditional-independence treatment-policy comparison remains unresolved, not proven novel. Track A and Duke Track B data are not linked at the patient level. Their relationship to resistance in Track C remains a simulator assumption.

## Relevant prior work

- An MRI-calibrated, biology-inspired high-grade glioma digital twin has been used to optimize radiotherapy under uncertainty in simulated patients ([Chaudhuri et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10598726/)).
- M4RL combines a multiscale glioblastoma tumor–macrophage model, a Fokker–Planck/PINN surrogate, and reinforcement learning for combination-inhibitor scheduling under modeled resistance ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)). Earlier work also applies Q-learning or deep Q-learning to simulated glioblastoma drug and chemoradiation schedules ([2020 temozolomide scheduling study](https://www.sciencedirect.com/science/article/pii/S0169260719316013); [Ebrahimi Zade et al., 2022](https://www.sciencedirect.com/science/article/pii/S1532046422000223)).
- A 2026 study combines longitudinal multiparametric MRI with a patient-specific reaction–diffusion model to identify candidate adaptive-radiotherapy targets ([Ty et al., 2026](https://ascopubs.org/doi/10.1200/CCI-25-00336)). A 2026 prospective pilot reports machine-learning-guided, personalized radiation-dose escalation from MRI-derived infiltration maps ([Nature Communications, 2026](https://www.nature.com/articles/s41467-026-72545-y)).
- Spatially matched MRI and biopsy/transcriptomic studies already connect imaging phenotypes with glioma biology ([Puchalski et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10539500/); [Plaisier et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760280/); [Zhou et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12685660/)). Longitudinal MRI recurrence localization also predates this project ([Chougule et al., 2022](https://pubmed.ncbi.nlm.nih.gov/34766380/)).
- Multimodal sensing, correlated observations, and robust policies are established in general POMDP/control research ([Choudhury et al., 2020](https://cdn.aaai.org/ojs/6645/6645-40-9874-1-10-20200521.pdf); [Bovy et al., 2025](https://elinebovy.github.io/assets/pdf/Multi-environment_POMDPs_discrete_model_uncertainty_under_partial_observability.pdf); [Zakeri et al., 2025](https://arxiv.org/abs/2503.12962)).

This is a targeted review, not a systematic search. It does not establish that the exact proposed policy contrast has never been studied. It does rule out novelty claims based on broad MRI/biology/PDE/RL integration alone.

## Narrow testable question

The most defensible remaining methods question is whether a dependence-aware controller changes simulated treatment decisions or outcomes relative to a matched conditional-independence controller when signal marginals, observation timing, treatment options, training access, and treatment-use comparison are held fixed. The frozen paired-signal predictive gate did not pass at the tested small cohort sizes, and existing identifiability work reports weak information about the dependence parameter. Any later policy result must therefore be conditional on the simulator assumptions; if dependence is not recoverable, report scenario-level policy sensitivity instead of attributing effects to one dependence coefficient.

The result could be a simulator-methods finding. It cannot validate a real molecular-to-resistance biomarker, establish patient-specific drug benefit, or show clinical value. Do not select parameters or seeds after viewing held-out outcomes to seek a positive result.

## Proposed future-scan sensitivity implementation

`RobustBeliefMPC` now accepts `assumed_future_scan_miss_probability`; `PolicyBenchmarkSettings` exposes it as `mpc_assumed_future_scan_miss_probability`. With the proposed `observed_two_scan` stop rule, MPC samples future scan availability under that fixed iid probability and uses the same availability draws for every candidate treatment pattern. A missed simulated scan resets two-scan confirmation. The planner never reads hidden scenario fields `miss_p` or `miss_p_hi`.

The default remains `0.0`, preserving the existing complete-follow-up assumption and seeded plans. Nonzero probabilities are allowed only with `stop_endpoint="observed_two_scan"`. This is an opt-in planner assumption, not a learned patient-specific missingness model. It does not represent state-dependent or informative scan loss inside MPC forecasts. Evaluation must predeclare assumed probabilities on development worlds and include a mismatch sensitivity against state-linked scan loss. This setting has not been used in a policy-outcome comparison and is not adopted as the primary protocol.

## Validation and next work

The scratch suite for the observed-stop, Q-learning, robust MPC, ambiguity-controller, and shared-policy benchmark modules passed 76 tests in 94.32 seconds under Python 3.14. These tests check settings validation, common future scan-availability draws, boundary probabilities, default seeded-plan preservation, and benchmark propagation. They establish code mechanics only. No research-scale multi-world policy outcome comparison, patient analysis, or sealed run occurred.

Before outcome scoring, freeze the stop rule, missing-scan sensitivity grid, utility, expected treatment-use comparison, evaluation-world block, precision target, interval method, and multiplicity rule. Then test parameter/set recovery and compare dependence-aware, matched conditional-independence, MRI-only, biology-only, no-signal, and equal-course fixed controls on paired independent development and held-out worlds. Report every predeclared regime, including ties and harm. Keep the umbrella title; make any result claim only about the exact frozen contrast.
