"""Protocol-complete LUMIERE wrapper (draft amendment A-2026-11-02-T §5g; not adopted). Synthetic use only.

Three prior-support runs on the same eligible patients and folds, sharing one FiLM reuse cache, plus a
central-gap refit (forecast gap 12-14 week labels) with its own empty FiLM cache. Edge-gap patients are left
out of every fit and score in the refit. If the primary run stops (fewer than 25 eligible), every run stops.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence

import numpy as np

from src import fk_lumiere_runner as lr
from src.fk_learned_baselines import UNET_RECIPE

PRIOR_SCENARIOS = {
    "primary": ((1e-5, 1e-1), (1e-5, 10.0)),
    "one_decade_narrower": ((1e-4, 1e-2), (1e-4, 1.0)),
    "one_decade_wider": ((1e-6, 1.0), (1e-6, 100.0)),
}
CENTRAL_GAP_WEEKS = (12, 13, 14)
CALIBRATION_MINIMUM_PATIENTS = 8
CALIBRATION_FOLDS = 5
INITIAL_BOUNDARY_COVARIANCE = {"sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0}
POSTERIOR_SETTINGS = {"n_coarse": 21, "mass_fraction": 0.9, "n_refine": 11, "n_draws": 200, "n_anchors": 256}
FORECAST_SEED, BOUNDARY_CALIBRATION_SEED = 20261028, 20261029


def default_provider_factory(rho_bounds, D_bounds):
    from src.fk_lumiere_pde_provider import LumierePDEProvider
    return LumierePDEProvider(
        settings={**POSTERIOR_SETTINGS, "rho_bounds": rho_bounds, "D_bounds": D_bounds},
        initial_covariance=INITIAL_BOUNDARY_COVARIANCE, cutoff_grid=lr.PDE_ONLY_CUTOFF_GRID,
        seed=FORECAST_SEED, boundary_calibration_seed=BOUNDARY_CALIBRATION_SEED,
        calibration_folds=CALIBRATION_FOLDS, calibration_minimum_patients=CALIBRATION_MINIMUM_PATIENTS)


def run_lumiere_protocol(cases: Sequence[lr.LumiereCase], fold_by_patient: Mapping[str, int], *,
                         provider_factory: Callable = default_provider_factory, gates=None, recipe=UNET_RECIPE,
                         n_bootstrap: int | None = None, minimum: int = lr.MIN_ELIGIBLE_PATIENTS) -> dict:
    prior_cache: dict = {}
    runs = {}
    for name, (rho_b, d_b) in PRIOR_SCENARIOS.items():
        runs[name] = lr.run_lumiere_cv(cases, fold_by_patient, pde_provider=provider_factory(rho_b, d_b),
                                       gates=gates, recipe=recipe, film_reuse_cache=prior_cache,
                                       minimum=minimum, n_bootstrap=n_bootstrap)
        if runs[name]["flow"]["stopped"]:
            return {"status": f"stopped in {name} run: {runs[name]['status']}", "runs": runs}
    sets = [set(r["dice"]["candidate"]) for r in runs.values()]
    if any(s != sets[0] for s in sets[1:]):
        raise RuntimeError("prior-support runs did not score the same patient set")
    for r in runs.values():
        for f, info in r["folds"].items():
            if info["test_ids"] != runs["primary"]["folds"][f]["test_ids"]:
                raise RuntimeError("prior-support runs did not keep the same folds")
    central = [c for c in cases if c.forecast_weeks in CENTRAL_GAP_WEEKS]
    gap_run = lr.run_lumiere_cv(central, fold_by_patient, pde_provider=provider_factory(*PRIOR_SCENARIOS["primary"]),
                                gates=gates, recipe=recipe, film_reuse_cache={}, minimum=1, n_bootstrap=n_bootstrap)
    base = {p: v for p, v in runs["primary"]["dice"]["candidate"].items()}
    sensitivity = {}
    for name, r in runs.items():
        sensitivity[name] = {
            "rho_bounds": list(PRIOR_SCENARIOS[name][0]), "D_bounds": list(PRIOR_SCENARIOS[name][1]),
            "pass_test_status": r["pass_test"]["status"],
            "mean_candidate_dice_change_from_primary": float(np.mean(
                [r["dice"]["candidate"][p] - base[p] for p in sorted(base)])),
            "n_pde_fallbacks": sum(v["forecast_failure_fallback"] for i in r["folds"].values() for v in i["pde"].values()),
            "film_reused_folds": sum(i["film_reused"] for i in r["folds"].values()),
        }
    return {"status": "completed (synthetic)" if all(c.synthetic for c in cases) else "completed",
            "runs": runs, "prior_support_sensitivity": sensitivity,
            "central_gap_run": gap_run, "central_gap_weeks": list(CENTRAL_GAP_WEEKS),
            "central_gap_excluded_ids": sorted(c.patient_id for c in cases if c.forecast_weeks not in CENTRAL_GAP_WEEKS),
            "interpretation_limit": "exploratory LUMIERE; not confirmation, novelty or clinical evidence"}
