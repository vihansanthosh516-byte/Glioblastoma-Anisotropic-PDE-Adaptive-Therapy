"""Run the frozen Track C full-harness timing-only pilot from a scratch source copy.

The benchmark result object is discarded. This script records no policy outcomes.
"""
from __future__ import annotations

import gc
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import os
import platform
import subprocess
import sys
import time
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from pathlib import Path

import trackc_calibration as calibration
import trackc_sim as tc
from trackc_integrated_signal_model import F0_REFERENCE_BOUNDS, SignalModelParameters, generate_integrated_signals


SOURCE_ROOT = Path(__file__).resolve().parent
AMENDMENT_RELATIVE = Path("paper/trackc_robust_policy_timing_amendment_2026-11-06.md")
OUTPUT_RELATIVE = Path("output/trackc_robust_policy_timing_pilot_20261106_followup/manifest.json")
SIGNAL_PARAMETERS = SignalModelParameters(q=0.7, a=0.75, b=0.75, d=0.75)
SEED_FAMILIES = (
    {
        "response_family": "single",
        "world_seed": 20261952,
        "scenario_and_signal_seeds": (20261953, 20261954, 20261955, 20261956, 20261957),
        "policy_seeds": (20261958, 20261959, 20261960, 20261961, 20261962),
    },
    {
        "response_family": "mixture",
        "world_seed": 20261963,
        "scenario_and_signal_seeds": (20261964, 20261965, 20261966, 20261967, 20261968),
        "policy_seeds": (20261969, 20261970, 20261971, 20261972, 20261973),
    },
)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): sha256_path(path)
        for path in sorted(root.rglob("*.py"))
    }


def load_source_module(name: str, filename: str):
    path = SOURCE_ROOT / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {filename}: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


benchmark = load_source_module("trackc_policy_benchmark_for_timing_followup", "211_trackc_policy_benchmark.py")
grid_module = load_source_module("trackc_signal_grid_for_timing_followup", "212_trackc_signal_fisher_grid.py")
MODEL_GRID = tuple(SignalModelParameters(*row) for row in grid_module.parameter_grid())
D0_GRID = benchmark.conditional_independence_grid(MODEL_GRID)
SETTINGS = benchmark.PolicyBenchmarkSettings(signal_f0_bounds=F0_REFERENCE_BOUNDS)


def _package_version(name: str):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _manifest(repo_root: Path, amendment_path: Path, script_path: Path) -> dict:
    commit = subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    return {
        "record_type": "timing_only_synthetic_pilot_followup",
        "status": "running",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_title": "Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization",
        "git_commit": commit,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "processor": platform.processor(),
            "packages": {name: _package_version(name) for name in ("numpy", "scipy", "scikit-learn", "torch")},
            "thread_environment": {key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
        },
        "amendment_sha256": sha256_path(amendment_path),
        "runner_script_sha256": sha256_path(script_path),
        "source_sha256": source_hashes(repo_root / "src"),
        "cohort_size_each": 2,
        "cohorts_per_world": ["selection", "calibration", "evaluation"],
        "signal_generation_parameters": {
            "q": SIGNAL_PARAMETERS.q,
            "a": SIGNAL_PARAMETERS.a,
            "b": SIGNAL_PARAMETERS.b,
            "d": SIGNAL_PARAMETERS.d,
        },
        "signal_f0_bounds": list(F0_REFERENCE_BOUNDS),
        "candidate_model_grid_size": len(MODEL_GRID),
        "matched_d0_model_grid_size": len(D0_GRID),
        "settings": {
            key: (list(value) if isinstance(value, tuple) else value)
            for key, value in vars(SETTINGS).items()
        },
        "families": [
            {
                "response_family": row["response_family"],
                "world_seed": row["world_seed"],
                "scenario_and_signal_seeds": list(row["scenario_and_signal_seeds"]),
                "policy_seeds": list(row["policy_seeds"]),
            }
            for row in SEED_FAMILIES
        ],
        "include_ambiguity_controller": True,
        "outcome_scores_saved": False,
        "outcome_scores_inspected": False,
        "returned_outcome_objects_discarded": True,
        "timing_seconds": {},
    }


def _run_family(row: dict) -> float:
    family = row["response_family"]
    world = calibration.draw_worlds(1, seed=row["world_seed"], family=family)[0]
    scenario_seeds = row["scenario_and_signal_seeds"]
    policy_seeds = row["policy_seeds"]
    selection_scenarios = calibration.draw_patients(world, 2, seed=[scenario_seeds[0], 1])
    calibration_scenarios = calibration.draw_patients(world, 2, seed=[scenario_seeds[1], 1])
    evaluation_scenarios = calibration.draw_patients(world, 2, seed=[scenario_seeds[2], 1])
    evaluation_observations, _ = generate_integrated_signals(
        2, q=SIGNAL_PARAMETERS.q, a=SIGNAL_PARAMETERS.a, b=SIGNAL_PARAMETERS.b, d=SIGNAL_PARAMETERS.d,
        seed=scenario_seeds[3], f0=evaluation_scenarios["f0"],
    )
    calibration_observations, _ = generate_integrated_signals(
        2, q=SIGNAL_PARAMETERS.q, a=SIGNAL_PARAMETERS.a, b=SIGNAL_PARAMETERS.b, d=SIGNAL_PARAMETERS.d,
        seed=scenario_seeds[4], f0=calibration_scenarios["f0"],
    )
    world_id = f"timing-followup-{family}"
    ids = {
        role: tuple(f"{world_id}-{role}-{index:02d}" for index in range(2))
        for role in ("selection", "calibration", "evaluation")
    }
    call_kwargs = {
        "world_id": world_id,
        "selection_scenarios": selection_scenarios,
        "calibration_scenarios": calibration_scenarios,
        "evaluation_scenarios": evaluation_scenarios,
        "selection_patient_ids": ids["selection"],
        "calibration_patient_ids": ids["calibration"],
        "evaluation_patient_ids": ids["evaluation"],
        "calibration_observations": calibration_observations,
        "evaluation_observations": evaluation_observations,
        "model_parameters": MODEL_GRID,
        "signal_generation_parameters": SIGNAL_PARAMETERS,
        "signal_generation_seed": scenario_seeds[3],
        "calibration_signal_generation_seed": scenario_seeds[4],
        "signal_generation_mri_template_sha256": None,
        "selection_scenario_seed": scenario_seeds[0],
        "calibration_scenario_seed": scenario_seeds[1],
        "evaluation_scenario_seed": scenario_seeds[2],
        "selection_seed": [policy_seeds[0]],
        "q_training_seed": [policy_seeds[1]],
        "calibration_seed": [policy_seeds[2]],
        "evaluation_seed": [policy_seeds[3]],
        "schedule_mixture_seed": [policy_seeds[4]],
        "settings": SETTINGS,
        "include_ambiguity_controller": True,
    }
    start = time.perf_counter()
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        result = benchmark.run_shared_policy_world(**call_kwargs)
    elapsed = time.perf_counter() - start
    del result, call_kwargs
    gc.collect()
    return elapsed


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        raise SystemExit("usage: 218_trackc_robust_policy_timing_pilot.py REPOSITORY_ROOT")
    repo_root = Path(argv[0]).resolve()
    source_root = repo_root / "src"
    output_path = repo_root / OUTPUT_RELATIVE
    amendment_path = repo_root / AMENDMENT_RELATIVE
    if output_path.exists():
        raise SystemExit(f"refusing to overwrite existing manifest: {output_path}")
    if source_hashes(SOURCE_ROOT) != source_hashes(source_root):
        raise SystemExit("scratch src copy does not match the repository src tree")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = _manifest(repo_root, amendment_path, Path(__file__).resolve())
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    try:
        for row in SEED_FAMILIES:
            elapsed = _run_family(row)
            manifest["timing_seconds"][row["response_family"]] = elapsed
            manifest["updated_utc"] = datetime.now(timezone.utc).isoformat()
            output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            print(f"{row['response_family']}: {elapsed:.3f} seconds")
        manifest["status"] = "completed"
        manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    except BaseException as exc:
        manifest["status"] = "failed"
        manifest["failure"] = {"type": type(exc).__name__, "message": str(exc)}
        manifest["failed_utc"] = datetime.now(timezone.utc).isoformat()
        output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"pilot failed: {type(exc).__name__}: {exc}")
        raise
    finally:
        output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"manifest: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
