> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Historical model-output summary: evidence corrections

Corrected 2026-10-19. The old “Months 1–4 complete” matrix lists scripts 19–34
(16 entries), not a verified 16/39 execution manifest. Its July-15 single-session,
GPU/VRAM/CUDA details and runtime values lack saved environment/command logs and
are withdrawn. The data used are a 15000-cell subsample, not a complete spatial
matrix; full cohort dimensions need input metadata. These are legacy simulations,
not the current RANO goal (`../../START_HERE.md`).

The audit searched archived output only. Several corresponding root artifacts
do exist, so blanket “missing” claims would be wrong. Read-only root array
headers: phenotypic_velocity.npy (15000,32), drift_vectors.npy (15000,32),
diffusion_tensors.npy (15000,32,32), te_matrix.npy (100,100). Shapes do not prove
execution timing, biological validity or a complete reproduction chain.

| Root artifact (under ../../output/) | Corrected evidence and scope |
|---|---|
| `grn_metrics.json` | n_edges=320; inferred out-degrees APOD 46, S100B 42, MT3 40. “685 significant edges” and the old gene-degree list are not supported by this graph summary. No causal biology or unreported multiplicity control is inferred. |
| `grn_bootstrap_ci.json` | Present in root output. The old 32/380 count is superseded by ledger's historical 37/380; exact interval selection and resampling method must accompany any inferential use. |
| `fk_metrics.json` | Selected simulation: numerical speed 3.623743 versus analytical 4.0, error 9.40642%; configured 20000 steps, completed 20010. Old 15.6% and 31.6% estimates are unsupported for this saved artifact and not interchangeable. |
| `invasion_summary.json` | Selected coupled simulation: free speed 26.592883 versus analytical 28.284271 um/hr, 5.98% relative error. This is not a clinical wave-speed measurement. |
| `aba_metrics.json`, `aba_analysis_results.json` | Tumor-cell snapshot necrosis 22.2/30.4/38.1% at 90/120/150; separate whole-grid final fraction 17.19%. Old 6.7% is not a verified current calibration result. Denominators and timepoints matter. |
| `single_ko_results.json`, `dual_ko_results.json` | 200 single-gene rows, 21 dual pairs. Highest single collapse is TMLHE 0.0407984; no positive Bliss score in the dual output. These are latent-model perturbation scores, not drug response. |
| `single_ko_ti.json`, `dual_ko_ti.json` | MT-CO2 model score 4.819325; dual range -0.978599 to +0.442764. Floor-based score field named therapeutic_index is not a pharmacologic TI; old all-dual=1 and healthy-collapse=1 table is withdrawn. |

The old NEB energy/attractor assertions and “saddle proof” are not established by
this report. The distinct SPIB output has internal model-consistency checks,
not independent biological validation. See `../TRACK_A_RESULTS.md` for corrected
CSGT, latent/saddle and missing-provenance details.

The old final artifact tree is removed: its mix of exact files, intended outputs,
obsolete filenames and shape assertions did not establish completeness. Markdown
`drug_gating_report.md` exists under this archive's output directory; root output
contains the named machine-readable artifacts above. `invasion_metrics.json`
and old `saddle_point_metrics.json` must be individually checked before reuse;
neither is replaced by another run's similarly named summary. Artifact presence
alone is not a successful-execution or reproducibility certificate.

P0 proposals to tune wave speed, necrosis, collapse or a TI threshold were
historical modeling ideas. Their old numerical targets were not verified
measurements or clinical acceptance rules. No threshold tuning or rerun is
authorized here. The Month-5 one-week schedule, “Ready to Launch” and clinical
validation/manuscript plan are superseded as active next steps. Current ledger
and START_HERE govern; missing execution/provenance evidence is not regenerated.
