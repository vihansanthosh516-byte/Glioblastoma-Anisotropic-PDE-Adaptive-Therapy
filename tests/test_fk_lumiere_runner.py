"""Synthetic tests for the guarded LUMIERE runner (draft amendment A-2026-11-02-T). No patient data."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.ndimage import distance_transform_edt

from src import fk_lumiere_runner as lr
from src import fk_time_aware_baselines as ta
from src.fk_learned_baselines import UNetRecipe

TINY = UNetRecipe(widths=(4, 8), epochs=1, n_seeds=1)
SHAPE = (16, 16, 16)
FOLDS = {f"s{i:02d}": i % 5 for i in range(26)}


def _ball(r, c=(8, 8, 8)):
    g = np.ogrid[tuple(slice(0, s) for s in SHAPE)]
    return sum(((g[i] - c[i]) * 2.0) ** 2 for i in range(3)) <= r ** 2


def _cases(n=26, loads=None, **kw):
    out = []
    for i in range(n):
        fw, hw = 11 + i % 5, 6 + i % 7
        def load(i=i, fw=fw):
            if loads is not None:
                loads.append(f"s{i:02d}")
            return _ball(6 + 0.05 * 7 * fw)
        out.append(lr.LumiereCase(f"s{i:02d}", _ball(4), _ball(6), np.ones(SHAPE, bool), fw, hw, True, load,
                                  synthetic=True, **kw))
    return out


def _folds(cases):
    return {c.patient_id: FOLDS[c.patient_id] for c in cases}


def _grid(c, empty=False):
    d = distance_transform_edt(~np.asarray(c.visit2, bool), sampling=2.0)
    u = np.zeros(SHAPE, np.float32) if empty else np.exp(-d / 4.0).astype(np.float32)
    return ({k: (u >= k).astype(np.float32) for k in lr.PDE_ONLY_CUTOFF_GRID},
            {k: int(not (u >= k).any()) for k in lr.PDE_ONLY_CUTOFF_GRID})


class FakeProvider:
    """Checks fold isolation on every call; one forecast call per held-out patient."""

    def __init__(self, empty_ids=(), fail_ids=()):
        self.empty_ids, self.fail_ids, self.fit_views, self.forecasts = set(empty_ids), set(fail_ids), [], []

    def fit_fold(self, view):
        ids = [v["case"].patient_id for v in view]
        assert len({FOLDS[i] for i in ids}) == 4  # the held-out fold is absent
        self.fit_views.append(ids)
        return {"model_ready": True, "training_ids": set(ids),
                "pde_only_training_rows": {v["case"].patient_id: _grid(v["case"], v["case"].patient_id in self.empty_ids)[0]
                                           for v in view}}

    def forecast(self, case, nuisance):
        assert case.patient_id not in nuisance["training_ids"]
        self.forecasts.append(case.patient_id)
        if case.patient_id in self.fail_ids:
            raise lr.ForecastFailure("synthetic failure")
        prob, empty = _grid(case, case.patient_id in self.empty_ids)
        v2 = np.asarray(case.visit2, bool)
        return {"mask_draws": np.stack([v2] * 2), "no_pde_ranker_mask_draws": np.stack([_ball(7)] * 2),
                "pde_only_inclusion_probability_by_cutoff": prob, "pde_only_empty_draws_by_cutoff": empty}


def _run(cases, provider=None, **kw):
    return lr.run_lumiere_cv(cases, _folds(cases), pde_provider=provider or FakeProvider(), recipe=TINY,
                             n_bootstrap=200, **kw)


def test_real_cases_refused_until_every_gate_passes():
    real = [replace(c, synthetic=False) for c in _cases()]
    with pytest.raises(PermissionError, match="gates not recorded"):
        _run(real)
    almost = {g: True for g in lr.LUMIERE_GATES[:-1]}
    with pytest.raises(PermissionError, match="scoring_authorized_by_user"):
        _run(real, gates=almost)
    with pytest.raises(PermissionError, match="mixed"):
        _run(_cases()[:-1] + real[-1:])


def test_empty_visit2_and_qc_failures_excluded_before_target_reads_and_stop_rule():
    loads = []
    cases = _cases(loads=loads)
    cases[0] = replace(cases[0], visit2=np.zeros(SHAPE, bool))
    cases[1] = replace(cases[1], qc_passed=False)
    res = _run(cases)
    fl = res["flow"]
    assert fl["n_empty_visit2"] == 1 and fl["n_qc_failed"] == 1 and fl["n_eligible"] == 24 and fl["stopped"]
    assert res["status"].startswith("stopped") and loads == [] and "dice" not in res


def test_common_patients_folds_shared_pde_call_count_arm_and_ablation():
    loads = []
    cases = _cases(loads=loads)
    cases[0] = replace(cases[0], visit2=np.zeros(SHAPE, bool))
    prov = FakeProvider()
    res = _run(cases, prov)
    assert res["status"] == "scored (synthetic)" and res["flow"]["n_eligible"] == 25
    assert "s00" not in loads and sorted(loads) == sorted(c.patient_id for c in cases[1:])
    ids = {c.patient_id for c in cases[1:]}
    for name in lr.COMPARATORS + ("candidate", "matched_volume_no_pde_ablation"):
        assert set(res["dice"][name]) == ids
    assert sorted(prov.forecasts) == sorted(ids)  # exactly one PDE call per held-out patient
    assert len(prov.fit_views) == 5
    for f, info in res["folds"].items():
        assert not set(info["training_ids"]) & set(info["test_ids"])
        assert set(info["test_ids"]) == {p for p in ids if FOLDS[p] == f}
        assert info["rate"]["n_training"] == len(info["training_ids"])
        assert info["pde_only_cutoff"]["cutoff"] in lr.PDE_ONLY_CUTOFF_GRID
        assert info["pde_only_cutoff"]["n_training"] == len(info["training_ids"])
    assert res["dice"]["copy_last"] == res["dice"]["candidate"]  # fake candidate draws = visit-2 core
    assert res["dice"]["matched_volume_no_pde_ablation"] != res["dice"]["pde_only"]
    assert set(res["count_arm"]["rows"]) == ids and "count" not in " ".join(res["dice"])
    row = next(iter(res["count_arm"]["rows"].values()))
    assert {"poisson_deviance", "abs_error", "brier_any_new", "predicted_count"} <= set(row)
    assert "pass_test" in res and res["pass_test"]["status"] in {"pass", "not passed", "incomplete"}


def test_empty_pde_only_stays_empty_and_failures_use_recorded_fallback():
    res = _run(_cases(), FakeProvider(empty_ids={"s03"}, fail_ids={"s04"}))
    p3 = res["folds"][FOLDS["s03"]]["pde"]["s03"]
    assert p3["pde_only_mask_empty"] and not p3["forecast_failure_fallback"]
    assert res["dice"]["pde_only"]["s03"] == 0.0  # frozen empty-mask rule against a non-empty target
    assert res["dice"]["copy_last"]["s03"] > 0.0
    p4 = res["folds"][FOLDS["s04"]]["pde"]["s04"]
    assert p4["forecast_failure_fallback"] and p4["failure_reason"] == "synthetic failure"
    for arm in ("pde_only", "candidate", "matched_volume_no_pde_ablation"):
        assert res["dice"][arm]["s04"] == res["dice"]["copy_last"]["s04"]


def test_available_empty_film_mask_is_scored_but_unavailable_model_is_missing(monkeypatch):
    monkeypatch.setattr(ta, "train_film_models", lambda pairs, intervals, fold, recipe: {"status": "ok"})
    monkeypatch.setattr(ta, "predict_film", lambda trained, x, iv, domain: np.zeros(SHAPE, bool))
    res = _run(_cases())
    assert all(v == 0.0 for v in res["dice"]["film_unet"].values())  # empty mask vs non-empty target
    # fake candidate == copy-last gives a zero-SE contrast; the point is that no arm is reported missing
    assert res["pass_test"].get("reason") != "missing arm for at least one scored patient"
    monkeypatch.setattr(ta, "train_film_models",
                        lambda pairs, intervals, fold, recipe: {"status": "unavailable" if fold == 3 else "ok"})
    res = _run(_cases())
    assert {p for p, v in res["dice"]["film_unet"].items() if v is None} == {p for p in FOLDS if FOLDS[p] == 3}
    assert res["pass_test"]["status"] == "incomplete" and res["folds"][3]["film_status"] == "unavailable"


def test_cutoff_selection_uses_training_dice_and_ties_go_larger():
    sup = np.ones(SHAPE, bool)
    tgt = _ball(6)
    same = {k: _ball(6).astype(np.float32) for k in lr.PDE_ONLY_CUTOFF_GRID}
    sel = lr.select_pde_only_cutoff([(same, sup, tgt)])
    assert sel["cutoff"] == 0.5 and sel["boundary_hit"]
    graded = {k: (_ball(6) if k == 0.2 else _ball(3)).astype(np.float32) for k in lr.PDE_ONLY_CUTOFF_GRID}
    assert lr.select_pde_only_cutoff([(graded, sup, tgt)])["cutoff"] == 0.2
    assert lr.select_pde_only_cutoff([])["cutoff"] is None


def test_missing_cutoff_makes_comparison_incomplete():
    class NoRows(FakeProvider):
        def fit_fold(self, view):
            return {**super().fit_fold(view), "pde_only_training_rows": {}}
    res = _run(_cases(), NoRows())
    assert res["pass_test"]["status"] == "incomplete" and not res["comparison_complete"]


def test_film_reuse_only_on_full_fingerprint_match():
    cases = _cases()
    cache = {}
    first = _run(cases, film_reuse_cache=cache)
    assert not any(f["film_reused"] for f in first["folds"].values())
    second = _run(cases, film_reuse_cache=cache)
    assert all(f["film_reused"] for f in second["folds"].values())
    assert second["dice"]["film_unet"] == first["dice"]["film_unet"]
    shifted = [replace(c, history_weeks=c.history_weeks + 1) for c in cases]
    third = _run(shifted, film_reuse_cache=cache)
    assert not any(f["film_reused"] for f in third["folds"].values())


def test_film_fingerprint_covers_recipe_config_and_intervals(monkeypatch):
    x = np.zeros((3, 4, 4, 4), np.float32)
    y = np.zeros((4, 4, 4), bool)
    env = {"device": "cpu"}
    base = ta.film_training_fingerprint(["a"], [(x, y)], [(13, 13)], fold=0, recipe=TINY, environment=env)
    assert base != ta.film_training_fingerprint(["a"], [(x, y)], [(14, 13)], fold=0, recipe=TINY, environment=env)
    assert base != ta.film_training_fingerprint(["a"], [(x, y)], [(13, 12)], fold=0, recipe=TINY, environment=env)
    assert base != ta.film_training_fingerprint(["a"], [(x, y)], [(13, 13)], fold=0,
                                                recipe=replace(TINY, epochs=2), environment=env)
    monkeypatch.setitem(ta.FILM_CONFIG, "film_hidden", 32)
    assert base != ta.film_training_fingerprint(["a"], [(x, y)], [(13, 13)], fold=0, recipe=TINY, environment=env)


def test_provider_wraps_only_typed_expected_model_failures():
    from src.fk_hybrid_forecast import ForecastVolumeCapacityError
    from src.fk_lumiere_pde_provider import LumierePDEProvider

    settings = {"rho_bounds": (1e-5, 1e-1), "D_bounds": (1e-5, 10.0), "n_coarse": 5,
                "mass_fraction": 0.9, "n_refine": 3, "n_draws": 2, "n_anchors": 8}
    covariance = {"sigma_b_mm": 0.1, "tau_b_mm": 0.0, "ell_mm": 2.0}
    nuisance = {"population_log_volume_slope_per_day": 0.0,
                "shrinkage_time_constant_days": 90.0,
                "volume_variance_fit": {"intercept": 0.0, "slope_per_day": 0.0}}
    masks = [np.ones((4, 4, 4), dtype=bool), np.ones((4, 4, 4), dtype=bool)]
    domain = np.ones((4, 4, 4), dtype=bool)

    def capacity_failure(**_kwargs):
        raise ForecastVolumeCapacityError("synthetic capacity failure")

    provider = LumierePDEProvider(settings=settings, initial_covariance=covariance,
                                  cutoff_grid=(0.1,), forecast_fn=capacity_failure)
    with pytest.raises(lr.ForecastFailure, match="capacity failure"):
        provider._call("p", masks, [0.0, 10.0], 20.0, domain, covariance, nuisance,
                       calibration_only=False, seed=1)

    def unexpected_error(**_kwargs):
        raise ValueError("unexpected programming error")

    provider = LumierePDEProvider(settings=settings, initial_covariance=covariance,
                                  cutoff_grid=(0.1,), forecast_fn=unexpected_error)
    with pytest.raises(ValueError, match="unexpected programming error"):
        provider._call("p", masks, [0.0, 10.0], 20.0, domain, covariance, nuisance,
                       calibration_only=False, seed=1)


def test_real_pde_provider_one_shared_call_per_heldout_patient_and_training_only_calibration():
    from src.fk_lumiere_pde_provider import LumierePDEProvider
    g = np.ogrid[tuple(slice(0, s) for s in SHAPE)]
    r2 = sum(((g[i] - 8) * 2.0) ** 2 for i in range(3))
    sup = r2 <= 15 ** 2
    cases = []
    for i in range(10):
        fw, hw = 11 + i % 5, 6 + i % 3
        r1, rr2 = 5.0 + 0.1 * i, 6.0 + 0.1 * i
        cases.append(lr.LumiereCase(f"s{i:02d}", r2 <= r1 ** 2, r2 <= rr2 ** 2, sup, fw, hw, True,
                                    lambda rr2=rr2: r2 <= (rr2 + 1.0) ** 2, synthetic=True))
    prov = LumierePDEProvider(
        settings={"rho_bounds": (1e-5, 1e-1), "D_bounds": (1e-5, 10.0), "n_coarse": 5, "mass_fraction": 0.9,
                  "n_refine": 3, "n_draws": 6, "n_anchors": 16},
        initial_covariance={"sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0},
        cutoff_grid=lr.PDE_ONLY_CUTOFF_GRID, calibration_minimum_patients=5)
    res = lr.run_lumiere_cv(cases, _folds(cases), pde_provider=prov, recipe=TINY, minimum=10, n_bootstrap=200)
    full = [p for p, cal in prov.calls if not cal]
    cal = [p for p, cal in prov.calls if cal]
    ready = {f: info["pde_model_ready"] for f, info in res["folds"].items()}
    for f, info in res["folds"].items():
        for p, row in info["pde"].items():
            assert row["forecast_failure_fallback"] or full.count(p) == 1
    # each patient is calibrated only while it is an outer-training patient (4 of 5 folds), never when held out
    for c in cases:
        assert cal.count(c.patient_id) <= 4
    assert len(full) == sum(1 for c in cases if ready[FOLDS[c.patient_id]])
    assert set(res["dice"]["pde_only"]) == {c.patient_id for c in cases}
