"""PDE provider for the guarded LUMIERE runner (draft amendment A-2026-11-02-T). Synthetic use only for now.

Mirrors the frozen MU outer-nuisance rules (`src/fk_mu_patient_runner.py`), which stay unchanged:
- volume nuisance from outer-training patients only (`fit_mu_training_volume_nuisance`);
- inner patient-held-out boundary calibration on outer-training patients only; any inner fallback or fewer
  than two valid residuals -> model not ready -> every held-out forecast uses the persistence fallback;
- expected model failures and empty history masks raise ForecastFailure (persistence fallback).
LUMIERE differences: days are nominal (0, 7 x history weeks, + 7 x forecast weeks), so the inner calibration
forecasts each training patient at its own nominal horizon; every PDE call carries the PDE-only cutoff grid, so
the inner calibration runs supply each training patient's raw cutoff masks (no extra solver run).
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from src.fk_boundary_calibration import run_cross_fitted_boundary_calibration
from src.fk_hybrid_forecast import hybrid_posterior_predictive_forecast
from src.fk_mu_patient_runner import _EXPECTED_MODEL_FAILURES, _stable_seed
from src.fk_mu_training_nuisance import fit_mu_training_volume_nuisance
from src.fk_patient_cv import ForecastFailure, PatientSeries

H_MM = 2.0


class LumierePDEProvider:
    def __init__(self, *, settings: Mapping, initial_covariance: Mapping, cutoff_grid: Sequence[float],
                 seed: int = 20261102, boundary_calibration_seed: int = 20261029,
                 calibration_folds: int = 5, calibration_minimum_patients: int = 8,
                 forecast_fn=hybrid_posterior_predictive_forecast):
        self.settings = dict(settings)
        self.initial_covariance = dict(initial_covariance)
        self.cutoff_grid = tuple(float(c) for c in cutoff_grid)
        self.seed, self.calibration_seed = int(seed), int(boundary_calibration_seed)
        self.calibration_folds, self.calibration_minimum = int(calibration_folds), int(calibration_minimum_patients)
        self.forecast_fn = forecast_fn
        self.calls = []  # (patient_id, calibration_only) for audit and tests

    @staticmethod
    def _days(case):
        d_h = 7.0 * case.history_weeks
        return np.array([0.0, d_h, d_h + 7.0 * case.forecast_weeks])

    def _call(self, patient_id, masks, days, target_day, domain, covariance, nuisance, *, calibration_only, seed):
        self.calls.append((patient_id, calibration_only))
        s = self.settings
        try:
            return self.forecast_fn(
                masks=list(masks), days=np.asarray(days, float), target_day=float(target_day),
                prior={"distribution": "log_uniform"}, h=H_MM, boundary_covariance=dict(covariance),
                population_log_volume_slope_per_day=nuisance["population_log_volume_slope_per_day"],
                shrinkage_time_constant_days=nuisance["shrinkage_time_constant_days"],
                volume_variance_fit=nuisance["volume_variance_fit"], domain=domain,
                rho_bounds=tuple(s["rho_bounds"]), D_bounds=tuple(s["D_bounds"]), n_coarse=int(s["n_coarse"]),
                mass_fraction=float(s["mass_fraction"]), n_refine=int(s["n_refine"]), n_draws=int(s["n_draws"]),
                n_anchors=int(s["n_anchors"]), seed=_stable_seed(seed, patient_id),
                calibration_only=calibration_only, pde_only_cutoff_grid=self.cutoff_grid)
        except ValueError as exc:
            if str(exc) in _EXPECTED_MODEL_FAILURES:
                raise ForecastFailure(str(exc)) from exc
            raise

    def fit_fold(self, training_view) -> dict:
        """training_view: [{'case': LumiereCase, 'target': mask}] for outer-training patients only."""
        domains = {v["case"].patient_id: np.asarray(v["case"].support, bool) for v in training_view}
        records = tuple(PatientSeries(v["case"].patient_id, self._days(v["case"]),
                                      (np.asarray(v["case"].visit1, bool) & domains[v["case"].patient_id],
                                       np.asarray(v["case"].visit2, bool) & domains[v["case"].patient_id],
                                       np.asarray(v["target"], bool) & domains[v["case"].patient_id]))
                        for v in training_view)
        try:
            volume = fit_mu_training_volume_nuisance(records)
        except ValueError as exc:
            return {"model_ready": False, "unavailable_reason": str(exc), "pde_only_training_rows": {}}
        cal_records = tuple(r for r in records if all(m.any() and not m.all() for m in r.masks))
        if len(cal_records) < self.calibration_minimum:
            return {"model_ready": False, "pde_only_training_rows": {},
                    "unavailable_reason": "too few training patients with three usable masks for boundary calibration"}
        captured, failed = {}, []

        def inner_forecast(history, inner_nuisance):
            out = self._call(history.patient_id, history.masks, history.days, history.forecast_day,
                             domains[history.patient_id], self.initial_covariance, inner_nuisance,
                             calibration_only=True, seed=self.calibration_seed)
            captured[history.patient_id] = out["pde_only_inclusion_probability_by_cutoff"]
            return out

        def inner_fallback(history, _n):
            failed.append(history.patient_id)
            return {"pre_boundary_probability": np.asarray(history.masks[-1], np.float32), "_calibration_fallback": True}

        try:
            cal = run_cross_fitted_boundary_calibration(
                cal_records, fit_nuisance=fit_mu_training_volume_nuisance, forecast=inner_forecast,
                n_folds=min(self.calibration_folds, len(cal_records)), seed=self.calibration_seed,
                minimum_history=2, spacing_mm=H_MM, n_anchors=int(self.settings["n_anchors"]),
                ell_candidates_mm=(5.0, 10.0, 20.0), forecast_horizon_days=None,
                fallback_forecast=inner_fallback, allow_empty_history_masks=True)
        except ValueError as exc:
            return {"model_ready": False, "unavailable_reason": f"boundary calibration unavailable: {exc}",
                    "pde_only_training_rows": captured, "calibration_failed_ids": failed}
        base = {"pde_only_training_rows": captured, "calibration_failed_ids": failed,
                "n_calibration_records": len(cal_records)}
        if sum(row["fallback_count"] for row in cal["fold_audit"]):
            return {**base, "model_ready": False,
                    "unavailable_reason": "inner cross-fitted calibration used persistence fallback"}
        if cal["n_residual_cases"] < 2:
            return {**base, "model_ready": False, "unavailable_reason": "fewer than two valid boundary residuals"}
        cov = cal["covariance_fit"]
        return {**base, **volume, "model_ready": True,
                "boundary_covariance": {k: cov[k] for k in ("sigma_b_mm", "tau_b_mm", "ell_mm")}}

    def forecast(self, case, nuisance) -> dict:
        if not nuisance.get("model_ready", False):
            raise ForecastFailure(str(nuisance.get("unavailable_reason", "outer-fold model unavailable")))
        sup = np.asarray(case.support, bool)
        masks = (np.asarray(case.visit1, bool) & sup, np.asarray(case.visit2, bool) & sup)
        if any(not m.any() for m in masks):
            raise ForecastFailure("empty history mask requires the prespecified persistence fallback")
        days = self._days(case)
        return self._call(case.patient_id, masks, days[:2], days[2], sup, nuisance["boundary_covariance"],
                          nuisance, calibration_only=False, seed=self.seed)
