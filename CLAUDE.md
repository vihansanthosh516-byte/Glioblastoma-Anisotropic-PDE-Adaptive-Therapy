# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Multi-track computational oncology research platform for glioblastoma treatment planning. Four tracks, 65+ numbered pipeline scripts:

| Track | Scripts | Domain |
|-------|---------|--------|
| **A — MSOS** | `src/01–41` | Single-cell multi-omics → causal GRN → virtual drug screening |
| **B — PDE Cohort** | `src/42–48` | Anisotropic diffusion PDE cohort, adaptive therapy, sensitivity analysis |
| **C — Digital Twin** | `src/49–65` | Inverse parameter estimation, robust MPC, RL adaptive steering, virtual cohort |
| **C-Extended** | `src/neural_pde/`, `src/rl/`, `src/sensing/`, `src/hil/` | FNO surrogate, virtual biosensors, PPO chronotherapy, HIL pump interface |

A real-patient validation study (`run_improved_aniso.py`, `cross_validation.py`, `isef_figures.py`) tests anisotropic vs. isotropic tensor modeling against 62 UCSF-PDGM patients.

## Environment Setup

```bash
python -m venv venv
source venv/bin/activate            # Linux/macOS
# .\venv\Scripts\Activate.ps1      # Windows PowerShell

pip install numpy scipy matplotlib pillow pandas
# Full pipeline additionally needs:
# scanpy anndata umap-learn torch torch-geometric gymnasium stable-baselines3
# networkx seaborn pymc pytensor nibabel dipy SALib plotly h5py scikit-learn
```

**Path configuration:** Most scripts read `GBM_PROJECT_ROOT` (defaults to `/mnt/c/Users/vihan/20206 science fair` on WSL, `C:/Users/vihan/20206 science fair` on Windows). Track A scripts also use `GBM_DATA_DIR` for the external multiomic-gbm dataset.

## Running Scripts

**Run a single numbered script:**
```bash
python -m src.cli 53          # runs src/53_spatial_metrics.py
python -m src.cli --list      # list all pipeline scripts
python -m src.cli --run-all   # run all sequentially
```

**Track B cohort pipeline (scripts 42→48→45):**
```bash
bash run_all.sh               # full chain
bash run_all.sh --month10     # only script 45 (assumes 42–44 already done)
# Windows PowerShell:
.\run_all.ps1
.\run_all.ps1 -Month10
```

**Track C key scripts:**
```bash
python src/51_inverse_parameter_estimation.py --test
python src/52_robust_mpc_controller.py --benchmark --n-mc 50
python src/58_rl_adaptive_steering.py
python src/64_virtual_cohort_simulation.py
python src/65_generate_final_report.py
```

**Validation study:**
```bash
python run_improved_aniso.py
python isef_figures.py
```

**Track C-Extended (Phase 15 quick test):**
```bash
python src/phase15_virtual_trial.py \
    --model-path output/phase13_ppo_chronotherapy/ppo_chronotherapy_final.zip \
    --n-patients 5 --max-episode-hours 12
```

## Testing

```bash
pytest tests/
```

9 test files, 56 tests covering: `test_inverse_estimation`, `test_joint_inverse_estimation`, `test_robust_mpc`, `test_spatial_metrics`, `test_multiomic_fusion`, `test_hybrid_controller`, `test_mu_glioma_loader`, `test_retrospective_validation`, `test_treatment_models`.

## Docker (Track B + C)

```bash
make build && make run        # full benchmark
make test                     # pytest inside container
make shell                    # interactive bash with data/output mounted
make multiomic                # ElasticNet multi-omic training
make uq                       # FNO-ensemble UQ trajectory
make clean                    # remove output/*
```

## Path Convention

All scripts must derive paths from `__file__`, not the working directory:

```python
from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
```

Never use relative paths like `Path("output")` or `"output/file.json"` — they break when scripts are run from a different working directory. This was the root cause of repeated failures in scripts 27–29 and is the most common regression pattern in this repo.

## Architecture Notes

**Core PDE:** 2D/3D anisotropic Fisher-Kolmogorov equation — `∂u/∂t = ∇·(D(x)∇u) + ρu(1−u/K) − γC(t)u`. The diffusion tensor D encodes white-matter tract anisotropy (10× parallel vs. perpendicular in Track B synthetic cohort).

**Key modules in `src/`:**
- `51_inverse_parameter_estimation.py` / `51_joint_inverse_estimation.py` — L-BFGS-B bounded optimization to estimate ρ and D from longitudinal imaging volumes
- `52_robust_mpc_controller.py` — uncertainty-aware MPC over confidence intervals (ρ±15%, D±15%), dynamic horizon 7–21 days
- `53_spatial_metrics.py` — DSC and Hausdorff Distance between simulated and observed tumor masks
- `58_rl_adaptive_steering.py` — RL agent (behavioral cloning from optimal controller)
- `neural_pde/fno_solver.py` — Fourier Neural Operator surrogate (>1000× faster than PDE solver)
- `sensing/virtual_sensor.py` — 4-modality virtual biosensor (MRI, PET, ctDNA, ICP)
- `hil/pump_interface.py` — hardware-in-the-loop pump interface with safety watchdog

**Output directory:** All generated artifacts (JSON, PNG, MD, NPZ) go to `output/`. Heavy `.npz` binary arrays are git-ignored; the JSON + PNG evidence trail is tracked.

**Script numbering collisions:** Several script numbers have multiple variants (e.g., `51_inverse_parameter_estimation.py` and `51_joint_inverse_estimation.py`, `52_dti_anisotropic_solver.py` and `52_robust_mpc_controller.py`). The CLI and `run_all.sh` use the canonical versions listed in their respective script arrays.

**Track A data dependency:** Scripts `01`–`41` require the external UCSC Cell Browser `multiomic-gbm` dataset (~1.4 GB) pointed to via `GBM_DATA_DIR`. Without it, Track A scripts will fail at data loading.

**`Requirements.txt` has unresolved merge conflict markers** — use the dependency list from `pyproject.toml` or the README installation instructions instead.

## Verification Protocol

After any script runs, verify the output actually exists on disk — do not trust console output alone:

1. Check the file exists: `Get-ChildItem output/<expected_file>` (PowerShell) or `ls output/<expected_file>` (bash)
2. Check the size is non-zero
3. For JSON outputs, read and confirm key values are in expected ranges

Common failure mode: a script prints "SUCCESS" but writes nothing — due to a relative path resolving to the wrong directory, a swallowed exception, or an early return. If `output/<expected>` is missing after a run, the script failed regardless of what it printed. Fix the path or exception, then rerun from scratch.

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).

## Outside review (ChatGPT checks this work)
The author gives Claude's overviews, numbers and claims to ChatGPT for an independent rating. Assume anything Claude writes will be checked.
- Every number must cite a file in `output/` or `vault/PAPER_FINDINGS_LEDGER.md`. No number from memory.
- State the weakest point of each result next to the result. Do not leave a limit out of a summary.
- Do not overclaim: say "simulation", "assumed", "post hoc", "exploratory" where they apply.
- When an overview is written for review, include the negative findings, as `paper/overview_for_review.md` does.
- When ChatGPT feedback comes back, check each claim against the files before acting. Do not agree or disagree by deference. Record what is wrong in `paper/chatgpt_review_response.md` and what is right too.
- Interim numbers (jobs still running) must be labelled interim. A reviewer may quote them.
