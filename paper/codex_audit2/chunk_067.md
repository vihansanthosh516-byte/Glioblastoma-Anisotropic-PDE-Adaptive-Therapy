# Chunk 067 — utility scripts

Static audit only. All eight files passed `python -m py_compile`; none was run. Cache directed to `paper/codex_audit2/.pycache`.

| Script | Check | Finding | Evidence | Severity |
|---|---|---|---|---|
| `src/create_master_dashboard.py` | compile/dependencies | PASS. Uses numpy/pandas/Plotly. | Compile exit 0; imports lines 13–20. | low |
| `src/create_master_dashboard.py` | paths/count/artifact | `BASE_OUTPUT=Path("output/batch_digital_twins")` and clinical path `Path("data/...")` are cwd-relative, contrary to `CLAUDE.md`. Header claims 183 twins; actual loaded number is determined dynamically from existing directories. Writes `master_cohort_dashboard.html` (optional PNG); checked expected HTML absent. | Lines 3–10, 23, 30, 91, 509–515; `output/batch_digital_twins/master_cohort_dashboard.html` missing. | high |
| `src/create_unified_viewer.py` | compile/dependencies | PASS. Uses numpy/pandas/Plotly. | Compile exit 0; imports lines 12–18. | low |
| `src/create_unified_viewer.py` | paths/count/artifact | cwd-relative `BASE_OUTPUT`; title claims 183 sims though actual count is runtime directory count. Writes HTML under batch output; expected `unified_4d_patient_viewer.html` absent. | Lines 3–8, 21, 28, 606–610. | high |
| `src/create_unified_viewer_simple.py` | compile/paths/artifact | PASS. `BASE_OUTPUT` is cwd-relative; expected HTML absent. | Compile exit 0; line 17; writes lines 605–609; expected `output/batch_digital_twins/unified_4d_patient_viewer.html` missing. | high |
| `src/create_unified_viewer_simple.py` | unsupported count | Docstring claims 190 simulations, while “unified viewer” count is obtained dynamically by globbing and no fixed count check exists. | Lines 3–6, 20–82. | medium |
| `src/debug_reward.py` | compile/scope | PASS. Debug helper imports `resistance_env` after adding cwd-relative `src` to sys.path; no declared output or seed. It is an interactive diagnostic, not a result-producing analysis. | Lines 1–3; `run_fixed` below. | medium |
| `src/dicom_loader.py` | compile/dependencies | PASS. Uses nibabel, numpy, scipy; nibabel declared in Requirements, scipy in pyproject. No output files hard-coded at module level; save functions accept caller paths. | Compile exit 0; imports 10–16; write functions lines 449, 485. | low |
| `src/dicom_loader.py` | claim/assumptions | Header says it supports BraTS, TCGA and standard neuroimaging; function has synthetic spherical seed fallback when imaging segmentation is absent. Such fallback must be identified as simulated input, not patient data. | Header lines 2–7; fallback comment around line 304. | medium |
| `src/dti_loader.py` | compile/dependencies | PASS. Imports numpy/scipy; DIPY imported lazily in tractography method and present in `Requirements.txt`. | Compile exit 0; lines 37–42, 446–461; requirements. | low |
| `src/dti_loader.py` | title/path | “Real-Patient DTI Integration” describes capability; it can also load/build synthetic inputs. CLI paths are supplied by caller, so actual patient provenance cannot be inferred from script. | Header lines 2–18; CLI lines 514–521. | medium |
| `src/forecast_metrics.py` | compile/dependencies | PASS. Pure metric functions use numpy/scipy, both declared. No outputs or result numbers. Empty-mask behavior is explicit in docstring and functions. | Compile exit 0; lines 1–10, 85. | low |
| `src/generate_summary_charts.py` | compile/path | PASS. Writes to cwd-relative `outputs/validation_summary_dashboard.png`, violating `CLAUDE.md`; artifact is absent. | Compile exit 0; lines 8, 55–57; expected file missing. | high |
| `src/generate_summary_charts.py` | hand-entered outcomes | Chart source labels a literal in-memory table “real validation metrics data”; seven patient IDs, day intervals, V0/V1 volumes and rho/D values are typed directly into code rather than read from a source file. Their provenance and accuracy could not be verified. | Lines 8–25 (data dict); path/output lines 55–57. | high |

No new locked-set timing or freeze logic appears in these utilities. Outputs were checked for presence only; no dashboard or chart was generated. The hard-coded chart values should be traced to source records before citing as measured validation results.
