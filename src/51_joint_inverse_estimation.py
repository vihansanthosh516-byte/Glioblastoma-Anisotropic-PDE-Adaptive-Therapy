"""Joint longitudinal inverse estimation for treatment-aware FK surrogates."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Sequence
import numpy as np
from scipy.optimize import minimize

try:
    from .treatment_aware_pde import TreatmentSchedule, treatment_aware_ode_model
except ImportError:  # direct script or import-by-path execution
    _treatment_path = Path(__file__).with_name("treatment_aware_pde.py")
    _spec = importlib.util.spec_from_file_location("treatment_aware_pde", _treatment_path)
    _treatment_module = importlib.util.module_from_spec(_spec)
    sys.modules["treatment_aware_pde"] = _treatment_module
    assert _spec.loader is not None
    _spec.loader.exec_module(_treatment_module)
    TreatmentSchedule = _treatment_module.TreatmentSchedule
    treatment_aware_ode_model = _treatment_module.treatment_aware_ode_model


def estimate_joint_parameters(
    volumes: Sequence[float],
    days: Sequence[float],
    schedule: TreatmentSchedule | None = None,
    initial_guess: tuple[float, float] = (0.02, 0.013),
    bounds: tuple[tuple[float, float], tuple[float, float]] = ((0.005, 0.1), (0.001, 0.05)),
    prior: tuple[float, float] = (0.02, 0.013),
    prior_scale: tuple[float, float] = (0.02, 0.01),
    regularization: float = 0.01,
) -> dict[str, Any]:
    """Fit rho and D to all adjacent intervals, penalizing implausible extremes."""
    v = np.asarray(volumes, dtype=float)
    t = np.asarray(days, dtype=float)
    if v.size < 3 or t.size != v.size or np.any(v <= 0) or np.any(np.diff(t) <= 0):
        raise ValueError("need >=3 positive volumes at strictly increasing days")

    def objective(params: np.ndarray) -> float:
        rho, diffusion = params
        predictions = [v[0]]
        for i in range(len(v) - 1):
            predictions.append(treatment_aware_ode_model(
                rho, diffusion, predictions[-1], t[i + 1] - t[i], schedule,
                start_day=float(t[i]),
            ))
        residual = (np.asarray(predictions[1:]) - v[1:]) / np.maximum(v[1:], 1.0)
        penalty = regularization * np.sum(((params - prior) / np.maximum(prior_scale, 1e-12)) ** 2)
        return float(np.mean(residual ** 2) + penalty)

    result = minimize(objective, np.asarray(initial_guess, dtype=float), method="L-BFGS-B", bounds=bounds)
    rho, diffusion = result.x
    return {
        "rho": float(rho), "D": float(diffusion), "objective": float(result.fun),
        "convergence": bool(result.success), "n_intervals": int(len(v) - 1),
    }
