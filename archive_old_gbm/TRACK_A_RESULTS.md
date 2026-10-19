> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Track A: historical outputs and evidence limits

Corrected 2026-10-19. This replaces the old “Verified Outputs” scratchpad.
Existing artifacts are not a complete reproduction bundle. No analysis was rerun.
Current goal: `../START_HERE.md`; corrected result index:
`../vault/PAPER_FINDINGS_LEDGER.md`. Output paths below refer to root `../output/`.

## Data and classification (scripts 01–16)

Read-only array headers show `nn_X.npy` shape (15000,2500),
`cgat/cvae_latent.npy` (15000,32), and `cgat/cvae_latent_full.npy` (140355,32).
These shapes do not verify the former full-cohort zone counts, all DE outputs,
training size/epoch/loss claims, label distribution or complete preprocessing.
Those claims need exact manifests/run metrics and are withdrawn pending evidence.
A wildcard such as `02_de_*.tsv` is not an identified evidence file.

`benchmark_comparison.tsv` records LR accuracy 0.6983, RF 0.7260, Transformer
0.4960 and Hybrid 0.5120. These are exploratory historical comparisons;
independent seed/split provenance is not supplied by the table. The former
attention-NN/scVI/NMF numeric leaderboard and graph-benefit conclusion are not
retained without their exact metrics and fair split comparison. Script 14 was
skipped and the old combined leaderboard was manually assembled.

`cgat/gat_metrics.json:accuracy` is 0.7873333 on the original random-cell split.
This is leaky for patient generalization; it does not demonstrate a biological
advantage. `cgat_leak_free.json` and `trackA_audit.json` (ledger 2o) report
patient-level gat_pca_clean mean accuracy 63.8%, SD 13.3 points. Donor composition,
fresh/frozen status and diagnosis are confounded with labels; cells-per-patient
alone reaches 64.6% and PCA32+LR 66.5%. Neither 78.7% nor the leaky 75.8% variant
is a validated biology result. The three cited old graph artifacts
`cgat/gat_edge_index.npy`, `gat_edge_attr.npy` and `gat_meta.json` are absent;
the old k/node/edge/attribute counts cannot be certified from those files.

## State, landscape and inferred network (scripts 17–26)

`csgt_metrics.json:kruskal_wallis` gives H=141.7171010 and p=1.6846998e-31.
`zone_stats` means: Healthy 0.586675, Periphery 0.649461, Core 0.646952.
Periphery is slightly above Core. Group differences do not prove a continuous
biological transition. The old entropy/proportion numbers need an underlying
numeric table, not only a PNG, and are omitted as unverified.

Phenotypic velocities, energies, drift and diffusion are model-derived quantities
on an arbitrary latent scale. Former mean/max summaries and runtime claims lack
an exact aggregate definition or retained numeric summary here and are withdrawn.
The old failed-Hessian rationale has no cited failure artifact.
`spib_saddle_point_metrics.json` reports a numerical transition saddle with energy
2.018768 and mixed eigenvalues under its model. Its `validation.overall_pass`
is an internal consistency flag, not biological proof. The trajectories were
synthetically integrated from a velocity field; they are not observed patient
trajectories. The reported state merging and training/frame counts require exact
metadata before reuse. Internal stable/unstable labels are not external validation.

`grn_metrics.json:n_edges` is 320. `master_switches` lists inferred out-degrees
APOD 46, S100B 42, MT3 40. These are inferred links/nodes, not causal master
switches. TE runtime/nonzero counts, PID method and pair/resample counts are
unverified here without their exact summaries. Ledger entry 1 retains the
historical 37/380 bootstrap count; an exact bootstrap artifact and uncertainty
method are required before using it as an independently verified result.

## Assumed invasion simulations (scripts 27–30)

The agent-based oxygen/necrosis rules and FK parameters are assumed, not fitted
clinical histology. Script-27 front.sum()/healthy.sum() is dimensionless; its old
physical-velocity claim is withdrawn (ledger 2g).
`aba_analysis_results.json` has two different necrosis summaries:
`invasion_kinetics.final_necrotic_fraction` 0.1719017 (whole-grid summary) and
`clinical_correlation.necrotic_fraction` 0.3812566 (tumor-cell summary).
These are not interchangeable. Ledger 2f separately reports tumor-cell snapshot
fractions 22.2/30.4/38.1% at steps 90/120/150 from `aba_metrics.json:cell_counts`.
No supported 22.9% value is retained.
`clinical_correlation.core_doubling_time_hours` is 47.4826 hours, not 87.3;
it falls outside the output's assumed 7–30-day comparison range.
That output's histology label and literature-range flags do not validate biology.

`fk_metrics.json` distinguishes configured n_steps=20000 from
steps_completed=20010, D=10 and r=0.4. Numerical speed is 3.623743 versus analytic
4.0; speed_error_percent=9.40642. `invasion_summary.json` uses a different coupled
r=0.8 assumption: free speed 26.592883 um/hr, analytic 28.284271, relative error
5.98%; average speed 11.707008, max_front_radius 359.210236. The old three-model
18.1–26.6 um/hr span is not numerical agreement between independent validations.
Boundary-limited simulated fronts and selected scale conversions are not measured
patient invasion rates; the unsupported clinical-range/velocity claims are removed.

## In-silico knockout and selectivity scores (scripts 31–34)

Full `single_ko_results.json` contains 200 rows: TMLHE collapse_score 0.0407984
is highest and SDE2 is explicitly 0.0. These are model-collapse scores, not drug
effects. All 21 `dual_ko_results.json` pairs were checked: none includes ZNF106;
first listed pair is S100A11+S100A6 (C_A=0.00397110, C_B=0.00337434,
bliss_synergy=-0.00212064). No pair has positive Bliss synergy in this output.
No inferential “above noise” claim is made without a noise model.

The field named `therapeutic_index` is a formula-dependent model selectivity
score, not a pharmacologic therapeutic index. The 0.001 floor and raw zero
healthy-collapse values strongly affect it. `single_ko_ti.json` reports MT-CO2
4.819325, CCL3L1 3.4400, MT-CO3 3.2744 and S100A8 2.423229; these do not mean
32x/10x/9x/5x tumor-selective efficacy. `dual_ko_ti.json` ranges from -0.978599
to +0.442764, including zero entries. No positive candidate-target recommendation
or clinical threshold follows. The separate assumed dose-response model ranks
S100A6 first by a different score; it must not be conflated with collapse ranking.

## TCGA observational modeling (scripts 35–41)

`tcga_km_summary.json` analyzes 518 cohort rows: 428 events, 90 censored,
KM median 432 days, CI [393,457]. The former 377.5-day median of all observed
follow-up times is not a survival median. Ledger A3 gives age HR=1.031,
CI [1.023,1.038], p=8.88e-16 in this cohort. CGGA external demographic analysis
(ledger 2k, `cgga_validation.json`) attenuates age to pooled HR=1.009, Holm
p=0.082. Null sex/subtype estimates do not prove no prognostic value. The old
“matches exactly” and “dominant predictor” claims are withdrawn.

`penalized_survival_metrics.json` describes the TCGA bulk-expression merge:
n=150, 120 events, c_index=0.6385432; not spatial transcriptomics.
Coefficients: S100A8 0.0433261, CCL3L1 0.0210972, S100A11 0.0044640,
age_at_diagnosis 0.0186967, S100A6 -2.00948e-7. This file contains no inferential
p-values/CIs for those coefficients; positive direction is exploratory, not
statistical significance or an actionable gene effect. Original raw matrix
20,530x172 and demographic counts need their exact input/merge manifest before
reuse. C-index validation/split provenance is not established by its scalar.

`subtype_recurrence_summary.json:subtype_stats` reports Proneural 37 patients /
405 days, Classical 39/388, Mesenchymal 48/338.5, Neural 26/261. These are recorded
descriptive follow-up medians, not demonstrated subtype-specific KM estimates or
a tested correlation of invasion score with clinical aggressiveness.
`dose_response_report.json` has 75 model-classified high-risk patients among 150;
monotherapy scores 5.4428/5.2412/4.6610/3.7183 are modeled gates, not therapeutic
passes or measured response.

The old Markdown reports are now `output/clinical_validation_report.md` and
`output/drug_gating_report.md` within this archive; root output references were
stale. Their prose is not independent evidence for missing statistics.
The former script-43 “real-data clinical validation” paragraph is removed:
ledger 2m identifies likely synthetic/defaulted zone-expression inputs without
an exact provenance link. Its assumed stromal PDE treatment outputs are not
Track A patient validation, and no spatial-expression/MRI connection is established.

No clinical claim, formal biological proof or full artifact-verification claim
is supported by this corrected historical note. Missing evidence is identified
rather than regenerated or silently transferred between cohorts.
