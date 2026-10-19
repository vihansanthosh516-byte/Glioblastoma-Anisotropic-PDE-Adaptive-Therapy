# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Read `START_HERE.md` first.** CURRENT PROJECT (2026-10-18): "probabilistic RANO" (measurement noise of RANO 2.0 volume lines and its effect on PFS and trial power; see `START_HERE.md`, `analysis_plan_v1.md` Amendments 14-21, `paper/phd_level_plan_probabilistic_rano.md`). The four tracks below are the OLDER GBM platform: background, not the goal. Old files that say "first", "novel" or "revolutionary" are not supported.

Older platform: four tracks, 159 numbered top-level source files (maximum prefix 145, with collisions), 208 top-level and 217 recursive source `.py` files, counted 2026-10-19; file counts do not establish safe runnable pipelines:

| Track | Scripts | Domain |
|-------|---------|--------|
| **A — MSOS** | `src/01–41` | Single-cell multi-omics → causal GRN → virtual drug screening |
| **B — PDE Cohort** | `src/42–48` | Anisotropic diffusion PDE cohort, adaptive therapy, sensitivity analysis |
| **C — Digital Twin** | `src/49–65` | Inverse parameter estimation, robust MPC, RL adaptive steering, virtual cohort |
| **C-Extended** | `src/neural_pde/`, `src/rl/`, `src/sensing/`, `src/hil/` | FNO surrogate, virtual biosensors, PPO chronotherapy, HIL pump interface |

Historical modeling scripts (`run_improved_aniso.py`, `cross_validation.py`, `isef_figures.py`) compare tensor assumptions on imaging-derived inputs. The former 62-patient validation claim has no direct result-key citation here and is withdrawn as a validation statement.

## Environment Setup

Historical setup examples below are not a verified clean environment. No supported Python/package lock is established here; ledger 2bo reports a partial-data Python 3.11 scratch environment. Resolve optional dependencies before any authorized run.

```bash
python -m venv venv
source venv/bin/activate            # Linux/macOS
# .\venv\Scripts\Activate.ps1      # Windows PowerShell

pip install numpy scipy matplotlib pillow pandas
# Full pipeline additionally needs:
# scanpy anndata umap-learn torch torch-geometric gymnasium stable-baselines3
# networkx seaborn pymc pytensor nibabel dipy SALib plotly h5py scikit-learn
```

**Path configuration:** `GBM_PROJECT_ROOT` and `GBM_DATA_DIR` are intended overrides in parts of the older code; support/defaults vary by script and have not been exhaustively verified. The required path convention remains derivation from `__file__`.

## Running Scripts

The following CLI, Track B/C, validation and Phase 15 commands are historical examples, not verified current pipelines. Their paths/options/input artifacts have not all been checked. Do not run project scripts during audit correction; many write output or train models. `--run-all` is not a safe audit command. Current RANO scripts are 140-145; frozen/spent analyses require their plan rules, not these examples.

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
pytest tests/ -q  # scratch copy only; never the working tree
```

Static inventory checked 2026-10-19: 16 test files and 108 test functions (including calculator tests test1-test7). This is not a fresh test run. Ledger 2bo is an author-reported historical 101-test result assembled across multiple runs with partial data; its transcripts are absent. The prior audit reports 108 passed/13 warnings in a Git-backed scratch clone (`paper/codex_audit_report.md`, section E). Modules include: `test_inverse_estimation`, `test_joint_inverse_estimation`, `test_robust_mpc`, `test_spatial_metrics`, `test_multiomic_fusion`, `test_hybrid_controller`, `test_mu_glioma_loader`, `test_retrospective_validation`, `test_treatment_models`, `test_baselines`, `test_metrics_stats`, `test_no_leakage`, `test_run_manifest`, `test_solver`, `test_solver_monotone`, `test_rano_calc`.

## Docker (Track B + C)

Historical target descriptions below are unverified against the current Makefile/compose definitions; they are not authorization to run benchmarks or delete output.

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

**Core PDE:** 2D/3D anisotropic Fisher-Kolmogorov equation — `∂u/∂t = ∇·(D(x)∇u) + ρu(1−u/K) − γC(t)u`. The diffusion tensor D encodes white-matter tract anisotropy (assumed 10× parallel vs. perpendicular in the historical Track B synthetic configuration; see `src/42_anisotropic_pde.py`, not measured patient anisotropy).

**Key modules in `src/`:**
- `51_inverse_parameter_estimation.py` / `51_joint_inverse_estimation.py` — L-BFGS-B bounded optimization to estimate ρ and D from longitudinal imaging volumes
- `52_robust_mpc_controller.py` — uncertainty-aware MPC over confidence intervals (ρ±15%, D±15%), dynamic horizon 7–21 days
- `53_spatial_metrics.py` — DSC and Hausdorff Distance between simulated and observed tumor masks
- `58_rl_adaptive_steering.py` — RL agent (behavioral cloning from optimal controller)
- `neural_pde/fno_solver.py` — Fourier Neural Operator surrogate; the uncited >1000× speed claim is withdrawn
- `sensing/virtual_sensor.py` — 4-modality virtual biosensor (MRI, PET, ctDNA, ICP)
- `hil/pump_interface.py` — hardware-in-the-loop pump interface; watchdog behavior requires implementation-specific tests and is not a clinical safety guarantee

**Output directory:** The convention is to write generated artifacts under `output/` using `__file__`-derived paths. Historical scripts 15, 18, 20, 49 and 58 have known cwd-relative exceptions. Tracking/ignore state is artifact-specific; existence does not establish provenance or validation.

**Script numbering collisions:** Several script numbers have multiple variants (e.g., `51_inverse_parameter_estimation.py` and `51_joint_inverse_estimation.py`, `52_dti_anisotropic_solver.py` and `52_robust_mpc_controller.py`). `src/cli.py:SCRIPTS` and the `run_all.sh` arrays select specific variants; these selections are historical and do not designate validated current pipelines.

**Track A data dependency:** Some Track A scripts require the external UCSC Cell Browser `multiomic-gbm` dataset via `GBM_DATA_DIR`; exact requirements vary. The former blanket 01-41 dependency and uncited dataset size are omitted. Inspect each script and input manifest before an authorized run.

**`Requirements.txt` had merge conflict markers; resolved 2026-10-18 by keeping the union of both sides (SALib, nibabel, dipy, scikit-learn, h5py).** `pyproject.toml` does not list these packages yet, so use `Requirements.txt` together with the install lines under Environment Setup.

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

## Past work is not automatically good (user rule, 2026-10-18)
- Do not assume something is good, correct or worth keeping because it was already built, committed, or reported in the ledger. Judge it on evidence and on whether it serves the current goal.
- Discarding the existing work and starting over is acceptable to the user. Say so plainly when it is the better choice. Do not defend sunk cost.
- Current goal (set 2026-10-18): a completely new, strong project in computational biology with substantial mathematics, for the Greater Austin Regional Science and Engineering Fair (GARSEF), with a February 2027 deadline. Chosen project: probabilistic RANO (see `START_HERE.md`). Findings must be backed by solid science. Claims of novelty must say how they were checked and what is unverified.
- Presentation rule from the user: show best real results, do not advertise failures, never fabricate. Answer questions on limits fully when asked.
