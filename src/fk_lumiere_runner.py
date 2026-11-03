"""Guarded LUMIERE runner under DRAFT amendment A-2026-11-02-T (not adopted). Synthetic use only for now.

Required comparators: copy-last, time-aware fixed growth, PDE-only and the FiLM time-aware U-Net. The PDE
candidate and PDE-only forecasts are injected callables; this module does not define the PDE-only mask.
The MU/Duke runner (`src/fk_mu_patient_runner.py`) is not used or changed.

Fail-closed rules:
- Real (non-synthetic) cases are refused unless every gate in LUMIERE_GATES is recorded True.
- Eligibility runs before any target read: QC failures and empty visit-2 cores are excluded from every arm,
  counted, and the run stops (no target read, no scoring) if fewer than MIN_ELIGIBLE_PATIENTS remain.
- Targets are lazy (`load_target`); they are read only for eligible patients after the stop check.
- Every arm uses the same eligible patients, the same fold per patient, the same support and the same target.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from src.fk_learned_baselines import UNET_RECIPE, unet_inputs
from src.fk_mu_scoring import dice_with_empty_rule
from src import fk_time_aware_baselines as ta

MIN_ELIGIBLE_PATIENTS = 25
DAYS_PER_WEEK = ta.DAYS_PER_WEEK
COMPARATORS = ("copy_last", "time_aware_fixed_growth", "pde_only", "film_unet")
LUMIERE_GATES = (
    "license_derived_outputs_stated",
    "deepbratumia_version_stated",
    "label0_meaning_stated",
    "school_written_approval_saved",
    "preflight_read_scope_approved",
    "registration_qc_run",
    "support_check_passed",
    "amendment_A_2026_11_02_T_adopted",
    "scoring_draft_adopted",
    "protocol_frozen_with_hashes",
    "scoring_authorized_by_user",
)


@dataclass(frozen=True)
class LumiereCase:
    patient_id: str
    visit1: np.ndarray
    visit2: np.ndarray
    support: np.ndarray
    forecast_weeks: int          # week(v3) - week(v2), published labels
    history_weeks: int           # week(v2) - week(v1), published labels
    qc_passed: bool
    load_target: Callable[[], np.ndarray]
    spacing_mm: tuple[float, float, float] = (2.0, 2.0, 2.0)
    synthetic: bool = False


def nominal_days(case: LumiereCase) -> dict:
    """7 x week-label difference: an approximation, not an exact interval."""
    d_h = DAYS_PER_WEEK * case.history_weeks
    return {"days": (0.0, d_h), "target_day": d_h + DAYS_PER_WEEK * case.forecast_weeks,
            "forecast_nominal_days": DAYS_PER_WEEK * case.forecast_weeks, "history_nominal_days": d_h}


def check_gates(cases: Sequence[LumiereCase], gates: Mapping[str, bool]) -> None:
    if all(c.synthetic for c in cases):
        return
    if any(c.synthetic for c in cases):
        raise PermissionError("refused: synthetic and real cases cannot be mixed")
    missing = [g for g in LUMIERE_GATES if gates.get(g) is not True]
    if missing:
        raise PermissionError(f"scoring refused: gates not recorded as passed: {missing}")


def eligibility(cases: Sequence[LumiereCase], *, minimum: int = MIN_ELIGIBLE_PATIENTS) -> dict:
    """Exclude QC failures, then empty visit-2 cores, before any target read. Stop below ``minimum``."""
    excluded, eligible = {}, []
    for c in cases:
        if not c.qc_passed:
            excluded[c.patient_id] = "registration QC failed"
        elif not np.asarray(c.visit2, bool).any():
            excluded[c.patient_id] = "empty visit-2 core"
        else:
            eligible.append(c)
    reasons = list(excluded.values())
    return {"eligible": eligible, "excluded": excluded, "n_input": len(cases), "n_eligible": len(eligible),
            "n_qc_failed": reasons.count("registration QC failed"),
            "n_empty_visit2": reasons.count("empty visit-2 core"),
            "stopped": len(eligible) < minimum, "minimum": minimum}


def run_lumiere_cv(cases: Sequence[LumiereCase], fold_by_patient: Mapping[str, int], *,
                   pde_only_forecast: Callable, candidate_forecast: Callable,
                   gates: Mapping[str, bool] | None = None, recipe=UNET_RECIPE,
                   film_reuse_cache: dict | None = None, minimum: int = MIN_ELIGIBLE_PATIENTS) -> dict:
    """Outer patient CV with the four proposed required comparators and the PDE candidate.

    ``pde_only_forecast`` and ``candidate_forecast`` are called as f(case, timing, training_cases) and must
    return a boolean mask on the case grid. ``training_cases`` carry loaded targets of outer-training patients
    only, for any training-only nuisance fit.
    """
    check_gates(cases, gates or {})
    flow = eligibility(cases, minimum=minimum)
    result = {"flow": {k: v for k, v in flow.items() if k != "eligible"}, "comparators": COMPARATORS}
    if flow["stopped"]:
        result["status"] = f"stopped: {flow['n_eligible']} eligible < {minimum}; no target read, no scoring"
        return result
    elig = flow["eligible"]
    missing = [c.patient_id for c in elig if c.patient_id not in fold_by_patient]
    if missing:
        raise ValueError(f"patients without a fold: {missing}")
    targets = {c.patient_id: np.asarray(c.load_target(), bool) for c in elig}
    for c in elig:
        if targets[c.patient_id].shape != np.asarray(c.support).shape:
            raise ValueError(f"{c.patient_id}: target and support grids differ")
    cache = {} if film_reuse_cache is None else film_reuse_cache
    dice = {name: {} for name in COMPARATORS + ("candidate",)}
    folds = {}
    for f in sorted({fold_by_patient[c.patient_id] for c in elig}):
        train = [c for c in elig if fold_by_patient[c.patient_id] != f]
        test = [c for c in elig if fold_by_patient[c.patient_id] == f]
        if not test:
            continue
        timed = [ta.TimedGrowthCase(c.patient_id, c.visit2, c.support, c.spacing_mm, c.forecast_weeks,
                                    targets[c.patient_id]) for c in train]
        sel = ta.select_rate(timed)
        pairs = [(unet_inputs(c.visit1, c.visit2, c.support), targets[c.patient_id] & np.asarray(c.support, bool))
                 for c in train]
        intervals = [(c.forecast_weeks, c.history_weeks) for c in train]
        ids = [c.patient_id for c in train]
        fp = ta.film_training_fingerprint(ids, pairs, intervals, fold=f, recipe=recipe)
        reused = fp in cache
        if not reused:
            cache[fp] = ta.train_film_models(pairs, intervals, fold=f, recipe=recipe)
        film = cache[fp]
        training_view = [dict(case=c, target=targets[c.patient_id]) for c in train]
        for c in test:
            sup = np.asarray(c.support, bool)
            tgt = targets[c.patient_id]
            timing = nominal_days(c)
            preds = {
                "copy_last": np.asarray(c.visit2, bool) & sup,
                "time_aware_fixed_growth": ta.fixed_growth_forecast(
                    c.visit2, sup, spacing_mm=c.spacing_mm,
                    radius_mm=ta.timed_radius_mm(sel["rate_mm_per_day"], c.forecast_weeks)),
                "pde_only": np.asarray(pde_only_forecast(c, timing, training_view), bool) & sup,
                "film_unet": (ta.predict_film(film, unet_inputs(c.visit1, c.visit2, sup),
                                              (c.forecast_weeks, c.history_weeks), domain=sup)
                              if film.get("status") == "ok" else None),
                "candidate": np.asarray(candidate_forecast(c, timing, training_view), bool) & sup,
            }
            for name, p in preds.items():
                dice[name][c.patient_id] = None if p is None else float(dice_with_empty_rule(p, tgt))
        folds[f] = {"training_ids": ids, "test_ids": [c.patient_id for c in test], "rate": sel,
                    "film_fingerprint": fp, "film_reused": reused, "film_status": film.get("status")}
    complete = all(v is not None for v in dice["film_unet"].values())
    result.update({
        "status": "scored (synthetic)" if all(c.synthetic for c in cases) else "scored",
        "folds": folds, "dice": dice,
        "patient_macro_dice": {k: (float(np.mean(list(v.values()))) if v and all(x is not None for x in v.values())
                                   else None) for k, v in dice.items()},
        "comparison_complete": complete,
    })
    return result
