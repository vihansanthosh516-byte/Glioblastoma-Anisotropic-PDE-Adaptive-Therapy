"""Result writer for the LUMIERE synthetic protocol implementation (draft amendment A-2026-11-02-T; not adopted).

Writes one JSON record of a `run_lumiere_protocol` result: code and protocol-document hashes, settings, cohort and
fold identifiers, gate status, all four analyses, missing-arm reasons and fallback counts. It never loads patient
data and stores no voxel arrays. It refuses to overwrite an existing file.
"""
from __future__ import annotations

import hashlib
import json
import platform
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CODE_FILES = (
    "src/fk_lumiere_runner.py", "src/fk_lumiere_protocol.py", "src/fk_lumiere_pde_provider.py",
    "src/fk_lumiere_stats.py", "src/fk_lumiere_result_writer.py", "src/fk_time_aware_baselines.py",
    "src/fk_hybrid_forecast.py", "src/fk_learned_baselines.py", "src/fk_fixed_growth.py",
    "src/fk_mu_scoring.py", "src/fk_patient_cv.py", "src/fk_boundary_calibration.py",
    "src/fk_mu_training_nuisance.py", "src/fk_mu_patient_runner.py", "src/fk_longitudinal_prep.py",
)
PROTOCOL_FILES = (
    "paper/track_b_lumiere_time_aware_amendment_draft_2026-11-02.md",
    "paper/track_b_lumiere_scoring_amendment_draft_2026-11-01.md",
    "paper/track_b_no_contact_dataset_plan_2026-10-31.md",
    "paper/research_plan_complete_project_appendix.md",
)
ANALYSES = ("primary", "one_decade_narrower", "one_decade_wider", "central_gap")


def _sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def _json_safe(x):
    if isinstance(x, Mapping):
        return {str(k): _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [_json_safe(v) for v in (sorted(x) if isinstance(x, set) else x)]
    if isinstance(x, np.ndarray):
        raise TypeError("the result writer stores no arrays")
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    return x


def _summarize_run(run: Mapping) -> dict:
    out = {"status": run.get("status"), "analysis": run.get("analysis"), "flow": run.get("flow")}
    if "dice" not in run:
        return out
    folds = run["folds"]
    pde_rows = [r for f in folds.values() for r in f["pde"].values()]
    missing = {}
    for f, info in folds.items():
        if info["film_status"] != "ok":
            missing.setdefault("film_unet", []).append(f"fold {f}: FiLM U-Net {info['film_status']}")
        if info["pde_only_cutoff"].get("cutoff") is None:
            missing.setdefault("pde_only", []).append(f"fold {f}: {info['pde_only_cutoff'].get('reason')}")
    out.update({
        "patient_dice": run["dice"],
        "patient_macro_dice": run["patient_macro_dice"],
        "pass_test": run["pass_test"],
        "comparison_complete": run["comparison_complete"],
        "missing_arm_reasons": missing,
        "fallback_counts": {"pde_forecast_failure_fallbacks": sum(r["forecast_failure_fallback"] for r in pde_rows),
                            "pde_only_empty_masks": sum(bool(r["pde_only_mask_empty"]) for r in pde_rows),
                            "folds_pde_not_ready": [f for f, i in folds.items() if not i["pde_model_ready"]]},
        "folds": {f: {"training_ids": i["training_ids"], "test_ids": i["test_ids"],
                      "time_aware_rate_mm_per_day": i["rate"]["rate_mm_per_day"],
                      "pde_only_cutoff": i["pde_only_cutoff"].get("cutoff"),
                      "pde_only_cutoff_boundary_hit": i["pde_only_cutoff"].get("boundary_hit"),
                      "film_status": i["film_status"], "film_reused": i["film_reused"],
                      "film_fingerprint": i["film_fingerprint"], "count_model_status": i["count_model_status"],
                      "pde_model_ready": i["pde_model_ready"], "pde_unavailable_reason": i["pde_unavailable_reason"],
                      "pde_per_patient": i["pde"]} for f, i in folds.items()},
        "count_arm": run["count_arm"],
    })
    return out


def write_lumiere_protocol_result(result: Mapping, path, *, cases: Sequence, fold_by_patient: Mapping[str, int],
                                  gates: Mapping[str, bool] | None = None) -> dict:
    from src import fk_lumiere_protocol as lp
    from src import fk_lumiere_runner as lr
    from src import fk_time_aware_baselines as ta
    from src.fk_learned_baselines import UNET_RECIPE
    from src.fk_lumiere_stats import BOOTSTRAP_SEED, LEVEL, N_BOOTSTRAP
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    synthetic = all(c.synthetic for c in cases)
    ids = sorted(c.patient_id for c in cases)
    runs = dict(result.get("runs", {}))
    if "central_gap_run" in result:
        runs["central_gap"] = result["central_gap_run"]
    record = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "label": "synthetic protocol implementation" if synthetic else "LUMIERE protocol run",
        "amendment_status": "A-2026-11-02-T and the LUMIERE scoring draft are NOT adopted; patient scoring blocked",
        "overall_status": result.get("status"),
        "code_sha256": {f: _sha(PROJECT_ROOT / f) for f in CODE_FILES},
        "protocol_sha256": {f: _sha(PROJECT_ROOT / f) for f in PROTOCOL_FILES},
        "settings": {
            "prior_scenarios": lp.PRIOR_SCENARIOS, "central_gap_weeks": lp.CENTRAL_GAP_WEEKS,
            "calibration_minimum_patients": lp.CALIBRATION_MINIMUM_PATIENTS, "calibration_folds": lp.CALIBRATION_FOLDS,
            "initial_boundary_covariance": lp.INITIAL_BOUNDARY_COVARIANCE, "posterior": lp.POSTERIOR_SETTINGS,
            "seeds": {"forecast": lp.FORECAST_SEED, "boundary_calibration": lp.BOUNDARY_CALIBRATION_SEED,
                      "bootstrap": BOOTSTRAP_SEED},
            "bootstrap": {"n": N_BOOTSTRAP, "level": LEVEL},
            "pde_only_cutoff_grid": lr.PDE_ONLY_CUTOFF_GRID, "time_aware_rate_grid_mm_per_day": ta.RATE_GRID_MM_PER_DAY,
            "film_config": ta.FILM_CONFIG, "unet_recipe": asdict(UNET_RECIPE),
            "minimums": lr.REQUIRED_MINIMUM, "comparators": lr.COMPARATORS,
        },
        "cohort": {"synthetic": synthetic, "n_input": len(ids), "patient_ids": ids,
                   "patient_ids_sha256": hashlib.sha256("\0".join(ids).encode()).hexdigest(),
                   "fold_by_patient": {p: int(fold_by_patient[p]) for p in ids}},
        "gates": {g: bool((gates or {}).get(g) is True) for g in lr.LUMIERE_GATES},
        "analyses": {name: (_summarize_run(runs[name]) if name in runs else {"status": "not run"}) for name in ANALYSES},
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
    }
    if "prior_support_sensitivity" in result:
        record["prior_support_sensitivity"] = result["prior_support_sensitivity"]
    safe = _json_safe(record)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(safe, indent=1))
    return safe
