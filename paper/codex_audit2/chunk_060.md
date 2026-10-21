# Chunk 060 — scripts 55–59

Static audit only. `python -m py_compile` passed for all eight files; none was executed. Bytecode cache directed to `paper/codex_audit2/.pycache`.

| Script | Check | Finding | Evidence | Severity |
|---|---|---|---|---|
| `src/55_dti_3d_ingestion.py` | compile/deps | PASS. Imports numpy, pandas, scipy, matplotlib; declared output artifacts include tensor profiles/geometry metrics, which exist. | Compile exit 0; imports lines 35–45; writes lines 1197–1251. | low |
| `src/55_dti_3d_ingestion.py` | path | Uses `Path("output")` in main and helper defaults, cwd-relative contrary to `CLAUDE.md`. | Lines 1178, 1197–1244. | high |
| `src/55_dti_3d_ingestion.py` | assumptions | PDE uses fixed grid, diffusivity, growth rate, Gaussian seed and simulated patient geometry. These are computational assumptions, not fitted patient parameters by default. | Lines 62–93, 537–543; header output description. | medium |
| `src/55_mechanics_hybrid_pde.py` | compile/deps | PASS. Uses numpy/scipy/matplotlib; phase2 metrics output absent. | Compile exit 0; imports lines 18–25; `output/phase2_mechanics_metrics.json` missing. | medium |
| `src/55_mechanics_hybrid_pde.py` | path/assumptions | Writes/reads `Path("output")` relative to cwd; fixed 128 grid, 100 mm domain, dt=.05, 30-day horizon, Young modulus 500 Pa and Poisson .45 are assumed simulation constants. | Lines 31–49, 62, 310–315, 344. | high |
| `src/56_anisotropic_3d_dti.py` | compile/path/artifact | PASS. Seeded NumPy orientation generation; writes `Path("output")` relative to cwd. Expected phase3 metrics JSON absent. | Compile exit 0; lines 136, 412–416, 465; `output/phase3_3d_dti_metrics.json` missing. | high |
| `src/56_anisotropic_3d_dti.py` | assumptions | Fixed tissue FA values, rho/diffusion ranges, 3D synthetic geometry and 30-day horizon; output is simulated anisotropic-vs-isotropic comparison, not patient validation. | Lines 35–41, 127–137, 253, 410–416. | medium |
| `src/56_anisotropic_dti_pde.py` | compile/path/artifact | PASS. Seed 42; writes relative `Path("output")`; expected phase3 metrics JSON absent. | Compile exit 0; lines 139, 447–451, 504; `output/phase3_dti_metrics.json` missing. | high |
| `src/56_anisotropic_dti_pde.py` | assumptions | Fixed carrying capacity, DTI guidance, tissue FA values and 30-day synthetic run; no independent patient fit evidenced. | Lines 35–45, 125–150. | medium |
| `src/56_multi_scale_framework.py` | compile/path/artifact | PASS. Writes summary under `Path("output")`; expected JSON absent. | Compile exit 0; lines 427–431; `output/multi_scale_framework_summary.json` missing. | high |
| `src/56_multi_scale_framework.py` | assumptions | Simulation defaults and first-orientation choice are hard-coded; results depend on synthetic scale/geometry setup, not patient-specific validation. | Lines 234, 263, 399–411. | medium |
| `src/57_virtual_therapy_pde.py` | compile/path | PASS. Writes to relative `Path("output")`; expected `phase4_virtual_therapy_metrics.json` absent. | Compile exit 0; lines 415–419, 457. | high |
| `src/57_virtual_therapy_pde.py` | header/assumptions | Header says fixes were applied and cites standard protocol sources, but code uses assumed radiation/chemo schedule and normalized doses, not patient-specific dosing or clinical validation. | Header lines 2–12; constants lines 35–55; output absent. | high |
| `src/58_rl_adaptive_steering.py` | compile/deps/path | PASS. `OUTPUT_DIR=Path("output")` cwd-relative; expected metrics JSON and plot exist. Imports Gymnasium/SB3/torch stack; those are not declared in root `Requirements.txt` or pyproject. | Compile exit 0; line 71; metrics writes 586–591; requirements-docker note excludes GPU torch. | high |
| `src/58_rl_adaptive_steering.py` | header claims/assumptions | “Final Optimized”, 40 episodes and <2 min are header assertions not verified by static check; fixed reward, PDE parameters and schedules make this a simulation result. No independent patient data evidence. | Header lines 2–9; constants lines 50–69. | medium |
| `src/59_sensitivity_analysis.py` | compile/deps/path | PASS. Project-rooted output and outputs exist; imports scipy and custom loader. | Compile exit 0; lines 33, 71, 83–84, 796–804; outputs present. | low |
| `src/59_sensitivity_analysis.py` | scope/assumptions | Header says 30 parameter sets × 2 protocols; sensitivity sampling seed 42. Results are conditional on fixed PDE/therapy model assumptions, and script itself does not establish independent validation. | Header lines 2–15; constants lines 45–83; sampler lines 102–108. | medium |

No script was run, so performance timings, numerical outputs and dynamic path behavior were not reproduced. These are historical simulated GBM/PDE/therapy scripts; held-out test or freeze provenance is not evidenced within this chunk.
