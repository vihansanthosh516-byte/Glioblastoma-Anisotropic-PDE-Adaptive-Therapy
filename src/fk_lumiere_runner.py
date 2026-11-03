"""Guarded LUMIERE runner under DRAFT amendment A-2026-11-02-T (not adopted). Synthetic use only for now.

Required comparators: copy-last, time-aware fixed growth, PDE-only and the FiLM time-aware U-Net. One injected
PDE call per patient serves both the candidate and PDE-only (shared posterior draws and solver run).
PDE-only = PDE posterior forecast with a training-selected concentration cutoff, without forecast-volume
matching or boundary perturbation: per outer fold one cutoff from PDE_ONLY_CUTOFF_GRID is chosen on
outer-training patients only (their calibration-only raw fields), ties to the larger cutoff, then applied to the
held-out patient's raw draws with the candidate's 0.5 inclusion vote and the common support. An empty
PDE-only forecast stays empty. A PDE ForecastFailure uses the frozen persistence fallback and is recorded.
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
from src.fk_patient_cv import ForecastFailure
from src import fk_time_aware_baselines as ta

MIN_ELIGIBLE_PATIENTS = 25
CENTRAL_GAP_MINIMUM_PATIENTS = 10  # draft §5i: at least 2 per outer fold on average
# Real-data runs must use these minimums; only synthetic runs may pass a smaller ``minimum``.
REQUIRED_MINIMUM = {"primary": MIN_ELIGIBLE_PATIENTS, "central_gap": CENTRAL_GAP_MINIMUM_PATIENTS}
DAYS_PER_WEEK = ta.DAYS_PER_WEEK
PDE_ONLY_CUTOFF_GRID = tuple(round(0.05 * i, 2) for i in range(1, 11))  # 0.05 ... 0.50
CUTOFF_TIE_TOLERANCE = 1e-9
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


def select_pde_only_cutoff(training: Sequence[tuple], *, grid=PDE_ONLY_CUTOFF_GRID) -> dict:
    """training: (probability_by_cutoff, support, target) per outer-training patient. Ties go to the larger cutoff."""
    grid = sorted(float(c) for c in grid)
    if not training:
        return {"cutoff": None, "reason": "no training patient with usable raw PDE fields", "n_training": 0}
    means = []
    for c in grid:
        means.append(float(np.mean([dice_with_empty_rule((prob[c] >= 0.5) & sup, tgt) for prob, sup, tgt in training])))
    best = max(means)
    chosen = max(c for c, m in zip(grid, means) if m >= best - CUTOFF_TIE_TOLERANCE)
    return {"cutoff": chosen, "grid": grid, "training_mean_dice": means, "n_training": len(training),
            "boundary_hit": chosen in (grid[0], grid[-1])}


def run_lumiere_cv(cases: Sequence[LumiereCase], fold_by_patient: Mapping[str, int], *,
                   pde_provider,
                   gates: Mapping[str, bool] | None = None, recipe=UNET_RECIPE,
                   film_reuse_cache: dict | None = None, minimum: int | None = None,
                   n_bootstrap: int | None = None, analysis: str = "primary") -> dict:
    """Outer patient CV with the four proposed required comparators, the PDE candidate, a separately named
    matched-volume no-PDE ablation and the frozen Poisson count arm (count metrics only, never a mask).

    ``pde_provider.fit_fold(training_view)`` sees outer-training patients (with targets) only and returns the
    fold nuisance, including ``pde_only_training_rows`` {patient: probability_by_cutoff} from its inner
    calibration runs. ``pde_provider.forecast(case, nuisance)`` returns one hybrid-forecast dict (candidate
    draws, no-PDE draws and the PDE-only cutoff grid from one solver call) or raises ForecastFailure.
    """
    from src.fk_learned_baselines import count_metrics, count_rate_per_day, fit_count_model, new_label_pieces
    from src.fk_lumiere_stats import N_BOOTSTRAP, max_t_pass_test
    check_gates(cases, gates or {})
    if analysis not in REQUIRED_MINIMUM:
        raise ValueError(f"analysis must be one of {sorted(REQUIRED_MINIMUM)}")
    if minimum is None:
        minimum = REQUIRED_MINIMUM[analysis]
    if not all(c.synthetic for c in cases) and minimum != REQUIRED_MINIMUM[analysis]:
        raise PermissionError(f"real-data {analysis} runs must use the {REQUIRED_MINIMUM[analysis]}-patient minimum")
    flow = eligibility(cases, minimum=minimum)
    result = {"flow": {k: v for k, v in flow.items() if k != "eligible"}, "comparators": COMPARATORS,
              "analysis": analysis}
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
    dice = {name: {} for name in COMPARATORS + ("candidate", "matched_volume_no_pde_ablation")}
    count_rows = {}
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
        nuisance = pde_provider.fit_fold(training_view)
        raw_rows = nuisance.get("pde_only_training_rows", {})
        cut = select_pde_only_cutoff([(raw_rows[c.patient_id], np.asarray(c.support, bool), targets[c.patient_id])
                                      for c in train if c.patient_id in raw_rows])
        cut["training_ids_without_raw_fields"] = [c.patient_id for c in train if c.patient_id not in raw_rows]
        # Frozen Poisson count arm: outer-training patients only, nominal-interval offset; count metrics only.
        vox_mm3 = float(np.prod(train[0].spacing_mm))
        count_model = fit_count_model(
            [new_label_pieces(c.visit2, targets[c.patient_id], spacing_mm=c.spacing_mm, domain=c.support)
             for c in train],
            [np.log(np.count_nonzero(np.asarray(c.visit2, bool) & np.asarray(c.support, bool)) * vox_mm3)
             for c in train],
            [DAYS_PER_WEEK * c.forecast_weeks for c in train])
        pde_info = {}
        for c in test:
            sup = np.asarray(c.support, bool)
            tgt = targets[c.patient_id]
            last = np.asarray(c.visit2, bool) & sup
            try:
                full = pde_provider.forecast(c, nuisance)
                cand = (np.asarray(full["mask_draws"]).mean(axis=0) >= 0.5) & sup
                no_pde = (np.asarray(full["no_pde_ranker_mask_draws"]).mean(axis=0) >= 0.5) & sup
                if cut["cutoff"] is None:
                    pde_only, empty_draws = None, None
                else:
                    pde_only = (full["pde_only_inclusion_probability_by_cutoff"][cut["cutoff"]] >= 0.5) & sup
                    empty_draws = int(full["pde_only_empty_draws_by_cutoff"][cut["cutoff"]])
                failed, reason = False, None
            except ForecastFailure as exc:
                cand, pde_only, no_pde = last.copy(), last.copy(), last.copy()  # frozen persistence fallback
                empty_draws, failed, reason = None, True, str(exc)
            pde_info[c.patient_id] = {"forecast_failure_fallback": failed, "failure_reason": reason,
                                      "empty_draws_at_cutoff": empty_draws,
                                      "pde_only_mask_empty": None if pde_only is None else not pde_only.any()}
            preds = {
                "copy_last": last,
                "time_aware_fixed_growth": ta.fixed_growth_forecast(
                    c.visit2, sup, spacing_mm=c.spacing_mm,
                    radius_mm=ta.timed_radius_mm(sel["rate_mm_per_day"], c.forecast_weeks)),
                "pde_only": pde_only,
                "film_unet": (ta.predict_film(film, unet_inputs(c.visit1, c.visit2, sup),
                                              (c.forecast_weeks, c.history_weeks), domain=sup)
                              if film.get("status") == "ok" else None),
                "candidate": cand,
                "matched_volume_no_pde_ablation": no_pde,
            }
            for name, p in preds.items():
                dice[name][c.patient_id] = None if p is None else float(dice_with_empty_rule(p, tgt))
            mu = count_rate_per_day(count_model, float(np.log(last.sum() * vox_mm3))) * DAYS_PER_WEEK * c.forecast_weeks
            y = new_label_pieces(c.visit2, tgt, spacing_mm=c.spacing_mm, domain=sup)
            count_rows[c.patient_id] = {"predicted_count": mu, "observed_count": int(y), **count_metrics(mu, y)}
        folds[f] = {"training_ids": ids, "test_ids": [c.patient_id for c in test], "rate": sel,
                    "film_fingerprint": fp, "film_reused": reused, "film_status": film.get("status"),
                    "pde_only_cutoff": cut, "pde": pde_info, "count_model_status": count_model["status"],
                    "pde_model_ready": bool(nuisance.get("model_ready", False)),
                    "pde_unavailable_reason": nuisance.get("unavailable_reason")}
    stats = max_t_pass_test(dice, COMPARATORS, n_bootstrap=n_bootstrap or N_BOOTSTRAP)
    result.update({
        "status": "scored (synthetic)" if all(c.synthetic for c in cases) else "scored",
        "folds": folds, "dice": dice,
        "patient_macro_dice": {k: (float(np.mean(list(v.values()))) if v and all(x is not None for x in v.values())
                                   else None) for k, v in dice.items()},
        "comparison_complete": stats["status"] != "incomplete",
        "pass_test": stats,
        "count_arm": {"rows": count_rows, "scope": "count metrics only; never a Dice mask",
                      **{f"mean_{k}": float(np.mean([r[k] for r in count_rows.values()]))
                         for k in ("poisson_deviance", "abs_error", "brier_any_new")}},
    })
    return result
