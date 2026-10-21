# Chunk 054 — scripts 17–22

Static audit only. `python -m py_compile` passed for all eight files; no scripts were run. Bytecode cache was directed inside `paper/codex_audit2/.pycache`.

| Script | Check | Finding | Evidence | Severity |
|---|---|---|---|---|
| `src/17_gradient_diagnostic.py` | py_compile/deps | PASS. Imports numpy, matplotlib, umap; `umap` is not declared in root Requirements.txt or pyproject.toml. | Compile exit 0; lines 1–7. | medium |
| `src/17_gradient_diagnostic.py` | paths/artifact | All data/output paths are cwd-relative, contrary to `CLAUDE.md`; its plot and three input arrays exist. | Lines 17–19, 85; `output/gradient_failure_analysis.png` exists. | high |
| `src/18_csgt_framework.py` | py_compile/deps | PASS. Imports numpy/matplotlib/umap/scipy; `umap` not declared in root manifests. | Compile exit 0; lines 1–9. | medium |
| `src/18_csgt_framework.py` | paths/artifact | Uses cwd-relative paths, contrary to `CLAUDE.md`; declared plot, metrics and transition-score outputs exist. | Lines 19–21, 116–120, 149–150; output files present. | high |
| `src/18_csgt_framework.py` | assumptions | Fixed `alpha=beta=0.5`, `sigma=1.0` parameters define model; they are choices, not measured effects. | Lines 31–38. | medium |
| `src/19_phenotypic_velocity.py` | py_compile/deps | PASS. Imports torch, numpy, matplotlib; these are not all declared in root manifests (torch absent). | Compile exit 0; lines 17–20; requirements-docker line 2 excludes GPU torch. | high |
| `src/19_phenotypic_velocity.py` | header/code/randomness | Header says 25x25 UMAP grid; implementation projects via PCA and bins to grid. Sampling uses unseeded `np.random.choice`, so plot sampling is not reproducible despite seeded PCA. | Header lines 4–8; lines 135, 215, 275; no numpy seed. | medium |
| `src/19_phenotypic_velocity.py` | paths/artifact | Uses cwd-relative paths contrary to `CLAUDE.md`; declared PNG/NPY outputs exist. | Lines 26–28, 279–290; outputs checked present. | high |
| `src/20_fokker_planck_solver.py` | py_compile/deps | PASS. Imports torch, numpy, matplotlib, sklearn; torch not declared in root manifests. | Compile exit 0; lines 19–23; requirements-docker line 2. | high |
| `src/20_fokker_planck_solver.py` | header/assumptions | Header calls landscape “CALIBRATED” and describes desired topology; priors/weights and KDE bandwidth are hard-coded design values. Static code alone does not establish calibration to measured biology or prove topology. | Header lines 2–9; function defaults lines 96–102; main overrides lines 260–269. | high |
| `src/20_fokker_planck_solver.py` | paths/artifacts | All data/output paths are cwd-relative contrary to `CLAUDE.md`; declared output arrays/plot exist. | Lines 29–57, 184, 241–243; output files present. | high |
| `src/21_drift_diffusion_analysis.py` | py_compile/deps | PASS. Imports torch and numpy; torch absent from root manifests. | Compile exit 0; lines 14–16. | high |
| `src/21_drift_diffusion_analysis.py` | paths/artifact | Cwd-relative inputs/outputs violate `CLAUDE.md`; drift and diffusion outputs exist. | Lines 22–25, 145–148; all four output arrays present. | high |
| `src/21_drift_diffusion_analysis.py` | randomness | No random seed found; stochastic/downstream behavior cannot be assumed reproducible. | Full-file search for `random_state`, `manual_seed`, `np.random`, `random.` returned no match. | medium |
| `src/22_saddle_point_proof.py` | py_compile/deps | PASS. Imports torch, sklearn, numpy; torch absent from root manifests. | Compile exit 0; lines 19–25. | high |
| `src/22_saddle_point_proof.py` | method/claims | Header says “saddle point proof”; code uses PCA 5D, selected kNN/bandwidth and eigenvalue noise thresholds. This is a computational diagnostic under chosen assumptions, not proof of a biological saddle. | Header lines 2–17; constants lines 29–33; output JSON `output/saddle_point_metrics.json` exists. | high |
| `src/22_saddle_point_proof.py` | paths/artifact | Reads/writes cwd-relative `output/...`, contrary to `CLAUDE.md`; metrics JSON exists. | Lines 385 and 503–506; output present. | high |
| `src/22_saddle_point_proof_backup.py` | py_compile/deps | PASS. Imports torch, sklearn, numpy; torch absent from root manifests. | Compile exit 0; lines 12–19. | high |
| `src/22_saddle_point_proof_backup.py` | status/path | Backup implementation still writes same `output/saddle_point_metrics.json` relative to cwd. It is an alternate legacy solver; if included in runnable inventory, duplicate output could overwrite the primary script result. | Lines 392–406 and 503–506. | medium |
| `src/22_spib_saddle_proof.py` | py_compile/deps | PASS. Imports torch, sklearn, numpy; torch absent from root manifests. | Compile exit 0; lines 11–16. | high |
| `src/22_spib_saddle_proof.py` | claims/artifacts | Header claims SPIB 3D kinetic coordinates; code loads `output/spib_rc_3d.npy` and writes metrics. Both input and `output/spib_saddle_point_metrics.json` exist. “Saddle point proof” remains stronger than what a chosen-coordinate Hessian classification establishes. | Header lines 2–8; lines 36–39, 290–292; metric JSON exists. | medium |
| `src/22_spib_saddle_proof.py` | paths/freeze | Uses cwd-relative outputs contrary to `CLAUDE.md`; consumes precomputed SPIB coordinates and labels. Their split/freeze provenance is not recoverable from this script. | Lines 36–39, 290–292. | high |

All scripts in this chunk belong to historical Track A. No held-out/locked data timing can be established from the scripts beyond their listed precomputed inputs; do not treat legacy descriptive plots or solver outputs as current-goal validation. Output presence confirms files exist only, not that rerunning the script would reproduce them.
