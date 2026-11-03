"""Synthetic tests for the guarded LUMIERE runner (draft amendment A-2026-11-02-T). No patient data."""
from dataclasses import replace

import numpy as np
import pytest

from src import fk_lumiere_runner as lr
from src import fk_time_aware_baselines as ta
from src.fk_learned_baselines import UNetRecipe

TINY = UNetRecipe(widths=(4, 8), epochs=1, n_seeds=1)
SHAPE = (16, 16, 16)


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
    return {c.patient_id: i % 5 for i, c in enumerate(cases)}


FOLDS = {f"s{i:02d}": i % 5 for i in range(26)}
CALLS = []


def _fake_result(c, calibration_only, empty=False):
    from scipy.ndimage import distance_transform_edt
    d = distance_transform_edt(~np.asarray(c.visit2, bool), sampling=2.0)
    u = np.zeros(SHAPE, np.float32) if empty else np.exp(-d / 4.0).astype(np.float32)
    out = {"pde_only_inclusion_probability_by_cutoff": {k: (u >= k).astype(np.float32) for k in lr.PDE_ONLY_CUTOFF_GRID},
           "pde_only_empty_draws_by_cutoff": {k: int(not (u >= k).any()) for k in lr.PDE_ONLY_CUTOFF_GRID}}
    if not calibration_only:
        out["mask_draws"] = np.stack([np.asarray(c.visit2, bool)] * 2)
    return out


def _pde(c, timing, view, calibration_only, empty_ids=(), fail_ids=()):
    ids = {v["case"].patient_id for v in view}
    assert c.patient_id not in ids
    assert timing["forecast_nominal_days"] == 7 * c.forecast_weeks
    assert timing["target_day"] == 7 * (c.history_weeks + c.forecast_weeks)
    present = {FOLDS[i] for i in ids}
    if calibration_only:
        present.add(FOLDS[c.patient_id])
        assert len(present) <= 4  # one outer fold (the held-out one) is always absent
    else:
        assert FOLDS[c.patient_id] not in present
    CALLS.append((c.patient_id, calibration_only))
    if c.patient_id in fail_ids and not calibration_only:
        raise lr.ForecastFailure("synthetic failure")
    return _fake_result(c, calibration_only, empty=c.patient_id in empty_ids)


def test_real_cases_refused_until_every_gate_passes():
    real = [replace(c, synthetic=False) for c in _cases()]
    with pytest.raises(PermissionError, match="gates not recorded"):
        lr.run_lumiere_cv(real, _folds(real), pde_forecast=_pde)
    almost = {g: True for g in lr.LUMIERE_GATES[:-1]}
    with pytest.raises(PermissionError, match="scoring_authorized_by_user"):
        lr.run_lumiere_cv(real, _folds(real), pde_forecast=_pde, gates=almost)
    mixed = _cases()[:-1] + real[-1:]
    with pytest.raises(PermissionError, match="mixed"):
        lr.run_lumiere_cv(mixed, _folds(mixed), pde_forecast=_pde)


def test_empty_visit2_and_qc_failures_excluded_before_target_reads_and_stop_rule():
    loads = []
    cases = _cases(loads=loads)
    cases[0] = replace(cases[0], visit2=np.zeros(SHAPE, bool))
    cases[1] = replace(cases[1], qc_passed=False)
    res = lr.run_lumiere_cv(cases, _folds(cases), pde_forecast=_pde)
    fl = res["flow"]
    assert fl["n_empty_visit2"] == 1 and fl["n_qc_failed"] == 1 and fl["n_eligible"] == 24 and fl["stopped"]
    assert res["status"].startswith("stopped") and loads == [] and "dice" not in res


def test_full_synthetic_run_scores_all_arms_on_common_patients(monkeypatch):
    loads = []
    cases = _cases(loads=loads)
    cases[0] = replace(cases[0], visit2=np.zeros(SHAPE, bool))
    res = lr.run_lumiere_cv(cases, _folds(cases), pde_forecast=_pde,
                            recipe=TINY)
    assert res["status"] == "scored (synthetic)" and res["flow"]["n_eligible"] == 25
    assert "s00" not in loads and sorted(loads) == sorted(c.patient_id for c in cases[1:])
    ids = {c.patient_id for c in cases[1:]}
    for name in lr.COMPARATORS + ("candidate",):
        assert set(res["dice"][name]) == ids
    for f, info in res["folds"].items():
        assert not set(info["training_ids"]) & set(info["test_ids"])
        assert set(info["test_ids"]) == {p for p in ids if _folds(cases)[p] == f}
        assert info["rate"]["n_training"] == len(info["training_ids"])
    assert res["dice"]["copy_last"] == res["dice"]["candidate"]  # fake candidate draws = visit-2 core
    for info in res["folds"].values():
        assert info["pde_only_cutoff"]["cutoff"] in lr.PDE_ONLY_CUTOFF_GRID
        assert info["pde_only_cutoff"]["n_training"] == len(info["training_ids"])


def test_empty_pde_only_stays_empty_and_failures_use_recorded_fallback():
    cases = _cases()
    f = lambda c, t, v, cal: _pde(c, t, v, cal, empty_ids={"s03"}, fail_ids={"s04"})
    res = lr.run_lumiere_cv(cases, _folds(cases), pde_forecast=f, recipe=TINY)
    p3 = res["folds"][FOLDS["s03"]]["pde"]["s03"]
    assert p3["pde_only_mask_empty"] and not p3["forecast_failure_fallback"]
    assert res["dice"]["pde_only"]["s03"] == 0.0  # frozen empty-mask rule against a non-empty target
    assert res["dice"]["copy_last"]["s03"] > 0.0
    p4 = res["folds"][FOLDS["s04"]]["pde"]["s04"]
    assert p4["forecast_failure_fallback"] and res["dice"]["pde_only"]["s04"] == res["dice"]["copy_last"]["s04"]


def test_cutoff_selection_uses_training_dice_and_ties_go_larger():
    sup = np.ones(SHAPE, bool)
    tgt = _ball(6)
    same = {k: _ball(6).astype(np.float32) for k in lr.PDE_ONLY_CUTOFF_GRID}
    sel = lr.select_pde_only_cutoff([(same, sup, tgt)])
    assert sel["cutoff"] == 0.5 and sel["boundary_hit"]
    graded = {k: (_ball(6) if k == 0.2 else _ball(3)).astype(np.float32) for k in lr.PDE_ONLY_CUTOFF_GRID}
    assert lr.select_pde_only_cutoff([(graded, sup, tgt)])["cutoff"] == 0.2
    assert lr.select_pde_only_cutoff([])["cutoff"] is None


def test_film_reuse_only_on_full_fingerprint_match():
    cases = _cases()
    cache = {}
    first = lr.run_lumiere_cv(cases, _folds(cases), pde_forecast=_pde,
                              recipe=TINY, film_reuse_cache=cache)
    assert not any(f["film_reused"] for f in first["folds"].values())
    second = lr.run_lumiere_cv(cases, _folds(cases), pde_forecast=_pde,
                               recipe=TINY, film_reuse_cache=cache)
    assert all(f["film_reused"] for f in second["folds"].values())
    assert second["dice"]["film_unet"] == first["dice"]["film_unet"]
    shifted = [replace(c, history_weeks=c.history_weeks + 1) for c in cases]
    third = lr.run_lumiere_cv(shifted, _folds(shifted), pde_forecast=_pde,
                              recipe=TINY, film_reuse_cache=cache)
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
