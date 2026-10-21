# Chunk 053 — script audit

Static review only; no script was executed. `python -m py_compile` passed for all eight listed scripts (cache directed to `paper/codex_audit2/.pycache`).

| Script | Check | Finding | Evidence | Severity |
|---|---|---|---|---|
| `src/143_anchored_trial_sim.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/143_anchored_trial_sim.py` | imports/outputs | numpy/pandas and repo script 142. Writes cells CSV and summary JSON; both exist. | Lines 18–22, 57–76; `output/anchored_sim/anchored_sim_cells.csv`, `anchored_sim_summary.json`. | low |
| `src/143_anchored_trial_sim.py` | claims/assumptions | Explicitly simulation only; anchor/conversion `s=RC%/1.96` from one published study, large lesions only. Code hard-codes four s values and N=50/100, D=8/12; summary labels not real trial evidence. | Lines 2–9, 32–33, 68–69. Outputs exist, but do not validate general clinical performance. | medium |
| `src/144_rano_timeline_viewer.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/144_rano_timeline_viewer.py` | paths/artifacts | Project-rooted output; viewer HTML and info JSON exist. It is a visualization, not an analysis model. | Lines 19, 150–162; `output/viewer_rano/` contains declared outputs. | low |
| `src/144_rano_timeline_viewer.py` | dependencies/data | numpy and stdlib; dynamically reads source data via its configured module. Seeded deterministic patient sampling (`random.Random(20261018)`). Exact cohort/freeze scope cannot be verified from this file alone. | Lines 10–18, 72–75, 140–162. | medium |
| `src/145_conversion_audit.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/145_conversion_audit.py` | paths/artifacts | Project-rooted; writes conversion table, cells, summary. All three exist. Reads script 142 and A19 cells. | Lines 22–26, 61, 75, 83, 111; `output/conversion_audit/` outputs exist; A19 CSV exists. | low |
| `src/145_conversion_audit.py` | claims/assumptions | Clearly labels simulation, assumed conversion from one repeatability study, independent normal error, volumetric component only, and not real-trial evidence. PU1–PU4 are computed comparisons against chosen grid and A19 values, not empirical clinical validation. | Lines 2–13, 77–101. | medium |
| `src/14_cgat_evaluate.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/14_cgat_evaluate.py` | dependencies/paths | Uses numpy/pandas/matplotlib/seaborn; `ROOT="."` makes inputs and outputs cwd-relative, violating `CLAUDE.md` path rule. | Lines 6–16, 18–21, 52–53, 79–80, 108–109, 127–138. | high |
| `src/14_cgat_evaluate.py` | artifacts/claims | Expected benchmark TSV, plots, and Markdown are absent. Summary writes fixed “12,000 train / 3,000 test” and 80/20 split; source code cannot establish those counts are the evaluated artifacts since required metrics inputs are also unavailable. | Lines 18–21, 52–53, 79–80, 108–109, 127–138; all named outputs absent. | high |
| `src/15_baseline_scvi.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/15_baseline_scvi.py` | dependencies/paths | Imports torch, scipy, sklearn, numpy; torch is not declared in root requirements/pyproject and Docker requirements explicitly exclude GPU torch. Writes relative `output/` paths. | Lines 1–10, 200–218; requirements-docker note at line 2. | high |
| `src/15_baseline_scvi.py` | split/seeding | Split uses `random_state=42`, but code trains/evaluates across epochs using test loader; static review indicates test-set monitoring during training, so reported test performance is not untouched. No global torch/numpy seed setup found. | Lines 121–127, 174–193; no seed call found. | high |
| `src/15_baseline_scvi.py` | artifacts | `output/scvi_latent.npy` exists; metrics output not verified in this pass. The latent is generated after evaluation from full dataset; it must not be treated as independent-test-only data. | Lines 205–218. | medium |
| `src/15_cgat_build_graph.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/15_cgat_build_graph.py` | dependencies/paths | sklearn/numpy; absolute WSL path `/mnt/c/...` hard-coded, not portable and not rooted from `__file__`. | Lines 7–21. | high |
| `src/15_cgat_build_graph.py` | artifacts/typed count | Requires full latent/label/patient/region arrays and writes graph arrays/meta; expected outputs absent. Comment asserts latent shape `(140355, 32)` but no validation evidence is bundled here. | Lines 31–34, 103–117; expected `output/cgat/*_full*` files absent. | medium |
| `src/16_baseline_nmf.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/16_baseline_nmf.py` | dependencies/paths | numpy/scipy/sklearn; inputs and outputs are cwd-relative `output/...`, contrary to `CLAUDE.md`. | Lines 1–8, 20–22, 81–83. | high |
| `src/16_baseline_nmf.py` | seed/artifacts | Uses random_state 42 for NMF, split and classifier. `output/nmf_fractions.npy` exists; metrics output not verified here. | Lines 35–60, 81–83. | low |
| `src/16_cgat_train.py` | py_compile | PASS. | Compile exit 0. | low |
| `src/16_cgat_train.py` | dependencies/paths | Imports torch and torch_geometric (not declared in root requirements/pyproject; Docker explicitly excludes them for Track A/GPU). Hard-coded WSL absolute root. | Lines 7–24; requirements-docker note line 2. | high |
| `src/16_cgat_train.py` | test leakage | Uses test split for per-epoch accuracy selection/checkpointing, despite logging these as validation metrics; this makes test set part of model selection. `random_state=42` controls split but no torch seed found. | Lines 71–78, 136–172. | high |
| `src/16_cgat_train.py` | artifacts | Requires graph/latent arrays and writes model, metrics, predictions and plots; all expected `output/cgat/*full*` outputs checked are absent. | Lines 26–35, 52–55, 166, 211–238. | medium |

No new clinical-result values were found hard-coded in scripts 143/145; their fixed numeric settings are simulation design assumptions, disclosed in labels/docstrings. Current-goal freeze relevance for legacy Track A scripts 14–16 is historical; their benchmark/test metrics should not be cited as fresh evaluation without independent artifact and split provenance.
