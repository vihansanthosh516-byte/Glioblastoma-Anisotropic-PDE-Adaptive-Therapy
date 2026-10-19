> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Historical GBM research: corrected evidence map

Corrected 2026-10-19 using audit chunk 018, source files, existing artifacts and
`../vault/PAPER_FINDINGS_LEDGER.md`. This replaces the former conflicting
abstracts and result tables. Current work is probabilistic RANO; see
`../START_HERE.md`. The ten-month sequence was a project narrative, not an
independently verified chronology or proof that an end-to-end pipeline completed.

## Scope and results

Tracks A, B and C explored expression classification, assumed tumor-growth models
and simulated treatment policies. Real imaging or fitted growth inputs do not
turn simulated therapy outcomes into clinical validation. No validated two-way
patient digital twin or clinical decision-support tool is established.

| Area | Corrected interpretation and source |
|---|---|
| Track A classification | Historical >94% accuracy does not establish patient generalization. The leak-free C-GAT estimate is 63.8% with substantial fold variation and donor/composition confounding. The 78.7% and 75.8% variants are not biology results. Ledger 2o, `../output/trackA_audit.json and ../output/cgat_leak_free.json` and the exact source files cited in ledger entry 1 govern the comparison. |
| Tumor state and landscape | Statistical zone differences do not prove a continuous biological transition or formally prove a saddle. Ledger entry 1 records H=141.717, p=1.68e-31 and Periphery T-score 0.649 versus Core 0.647. These are computational descriptors. |
| Regulatory network | Ledger entry 1 records 320 directed edges, APOD out-degree 46, S100B 42 and MT3 40. Only 37/380 candidate edges survive the bootstrap check. Directed transfer entropy does not establish gene causality. |
| Invasion and drug screens | Ledger entry 1 and 2h identify model-derived speeds and sources. Script 29 free-front speed is 26.5928833 um/hr under selected parameters, not an observed patient velocity. No dual pair reaches TI>10; the S100A8/S100A11/LST1 patient-risk correlation and clinical-range claim are withdrawn without a matching output/source. |
| Survival | The historical Cox concordance 0.639 (n=150) is an internal model result, not clinical validation. See ledger entry 1 and its output references; no validated gene-to-treatment link follows. |
| Track B geometry | Simulated anisotropy comparison: fractal dimensions 1.21 versus 1.03, difference +0.178, 98/103 comparisons, dz=2.89. The fractal change is in the unfavorable direction. Ledger entry 1 identifies the exact artifact and scope. |
| Sensitivity | Historical Sobol result uses N=128, alpha_sens about 0.43, rho about 0.20, D_w about 0.000. It is an assumed-parameter model decomposition, not empirical causal attribution; see ledger entry 1. |
| Adaptive therapy | Separate configurations cannot be pooled into one saving/non-inferiority result. The legacy eight-case 3D simulation reports 81.4% ±30.4% mean sparing; parameter-regime sensitivity and final-volume disadvantages are in ledger 2e and `../output/3d_extension_summary{,_swanson}.json`. Zero thresholded volume does not mean eradication. |
| Resistance and adverse outcomes | Under the script-44 assumed model, adaptive progresses earlier in 10/61 cases and has higher final mass in the cited ratio-1000 summary. Ledger 2n and `../output/emax_sweep_summary.json` give endpoint-specific counts. No registered non-inferiority result is established. |
| Inverse estimation | In the fitted cohort rho reaches its lower bound in 124/154 and D reaches a bound in 147/154. Two volumes do not reliably identify D. Synthetic/noisy RMSE headlines do not establish patient parameter accuracy. See `../README_STATUS.md` and ledger source references. |
| Track C RL | The 25% improvement is from simulated cohort64. The 20/21 result uses fitted growth rates with assumed kill rates and measures day-90 volume, not TTP; a one-line rule matches it. See `../README_STATUS.md` and ledger entry 1. The former 13.94/11.01 pair is not the same comparison. |

The old untraceable claims for risk correlation, 68.1% saving and its t/CI,
tract correlation, 30% MPC cost-variance reduction, 203 mm3/59.8% single 3D run,
50-trial MPC saving, ablation percentages and 20-case non-inferiority statistics
are withdrawn. Their exact cohort, configuration and saved comparison artifacts
were not identified. A nonsignificant superiority test does not prove
non-inferiority. No equivalent-control or treatment-benefit claim is retained.

The old spatial headline DSC 0.85 ±0.08 and HD 3.2 ±1.1 mm conflicts with the
archived dossier's DSC 0.21 ±0.02 and HD 26.3 ±3.3 mm
(`output/RESEARCH_EXECUTIVE_SUMMARY_AND_TECHNICAL_DOSSIER.md`, abstract and metric
table). The latter are historical report values, not newly verified patient
metrics. They fail the old document's named thresholds; those thresholds are
not established clinical acceptance criteria here.

## Extended modules and tensor comparison

`../src/neural_pde/` contains a synthetic reaction-target generator, trainer and
inference-timing utility. No measured FD comparator, PDE-equivalence accuracy or
>1000x speedup is established. Virtual sensor modalities and noise are assumed,
not validated sensor performance. PPO/circadian code does not establish clinical
control. `../src/hil/pump_interface.py` contains protocol examples and an emulated
watchdog; physical pump validation and end-to-end latency guarantees are absent.
The old CPU timing percentages lack a retained benchmark artifact.
`../src/phase15_virtual_trial.py` repeats a fixed-seed environment; its iteration
count does not establish independent patients or a clinical trial.

The old three-way UCSF-PDGM comparison is an exploratory retrospective image/model
comparison. Its claimed 62-case eligibility (WHO grade 2–3, IDH-mutant), seed-42
42/20 split, per-patient DSC/p-value tables and three figure references lack a
linked, verified cohort manifest and result bundle in this README. Those numeric
tables and figure-availability claims are withdrawn pending those exact sources.
The drivers `../run_improved_aniso.py` and `../cross_validation.py` depend on local
external inputs. Heuristic eigenvalue multipliers are not measured atlas tensors
or a validated kurtosis correction. No causal inference or tensor recommendation
follows, even if a within-cohort comparison were reproduced.

The former “clinical decision rules” are withdrawn. Ledger 2 gives the distinct
script-62 modeled threshold rho*=0.0751/day, CI [0.0746,0.0820], not the old
0.024/day claim (`../output/biomarker_stability_metrics.json`,
`model_threshold.rho_star_per_day`, `bootstrap.interpolated.ci95`). Ledger 2b/2j
separately label the rho>0.02 adaptive-risk cutoff post hoc. Neither is an
independently validated patient-care rule.

## Artifact and reproducibility limits

Archived `output/` holds Markdown reports; binary/JSON artifacts can still exist
in root `../output/`. Audit chunk 018's broad missing-artifact claim is too strong:
on 2026-10-19 root `master_cohort_synthesis.png` is 4081x3644 pixels and root
`master_cohort_summary.json` is 40,875 bytes. Presence alone does not establish
reproducibility, eight-patient schema checks or provenance certification.
The old 3821x3503 dimensions are wrong. Saved matching rerun hashes, schema-check
transcripts and an exact run/configuration manifest are needed before retaining
byte-identical rerun or complete-audit claims. No rerun was performed.

`../src/48_3d_extension.py` exists (contrary to one audit row). It specifies a
50-cubed default grid and assumed reaction-diffusion/therapy model. The former
2D/3D equation, boundary and CFL table is removed as unverified generalization;
source, not that table, defines implemented discretization. Numerical model
specification does not validate tumor biology.

The old repository tree, 105-source-file count, nine-file/56-test count, five
viewer count, and exactly-32-local-NPZ availability assertion are superseded.
Current dated source/test inventory and scratch results are in `../CLAUDE.md`
and `../paper/codex_audit2/fix_progress.md`. No claim is made that every heavy
artifact remains on disk or that git ignores certify reproducibility.

Historical setup, Docker and reproduction commands are removed as current
instructions: clean installation and complete pipeline execution were not
verified. Metadata declares Python >=3.14; that is not a tested compatibility
range. External UCSC expression and UCSF imaging inputs require exact manifests,
access terms and configured local paths; the former dataset-size claim is
unverified. Consult current guidance before any run. No rerun of spent/frozen
sets is authorized by this archive. The old clinical roadmap remains an
unvalidated historical proposal, not the current next step.

`../LICENSE` contains the MIT text with copyright 2026 Vihan Santhosh. This does
not relicense external datasets. The former “Vihan et al.” publication citation
is withdrawn: no verified publication or author-list record is linked. Cite a
specific repository revision with confirmed authorship instead. Unsourced
survival epidemiology and categorical clinical-protocol descriptions are omitted.

Computational research only. Retrospective imaging analysis is distinct from
prospective clinical validation. No patient-care use or scientifically validated
platform is established by these documentation/code corrections.
