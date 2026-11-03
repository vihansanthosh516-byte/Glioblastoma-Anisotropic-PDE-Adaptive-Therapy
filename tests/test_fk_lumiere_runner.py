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


def _pde(c, timing, train):
    assert c.patient_id not in {t["case"].patient_id for t in train}
    assert timing["forecast_nominal_days"] == 7 * c.forecast_weeks
    assert timing["target_day"] == 7 * (c.history_weeks + c.forecast_weeks)
    return np.asarray(c.visit2, bool)


def test_real_cases_refused_until_every_gate_passes():
    real = [replace(c, synthetic=False) for c in _cases()]
    with pytest.raises(PermissionError, match="gates not recorded"):
        lr.run_lumiere_cv(real, _folds(real), pde_only_forecast=_pde, candidate_forecast=_pde)
    almost = {g: True for g in lr.LUMIERE_GATES[:-1]}
    with pytest.raises(PermissionError, match="scoring_authorized_by_user"):
        lr.run_lumiere_cv(real, _folds(real), pde_only_forecast=_pde, candidate_forecast=_pde, gates=almost)
    mixed = _cases()[:-1] + real[-1:]
    with pytest.raises(PermissionError, match="mixed"):
        lr.run_lumiere_cv(mixed, _folds(mixed), pde_only_forecast=_pde, candidate_forecast=_pde)


def test_empty_visit2_and_qc_failures_excluded_before_target_reads_and_stop_rule():
    loads = []
    cases = _cases(loads=loads)
    cases[0] = replace(cases[0], visit2=np.zeros(SHAPE, bool))
    cases[1] = replace(cases[1], qc_passed=False)
    res = lr.run_lumiere_cv(cases, _folds(cases), pde_only_forecast=_pde, candidate_forecast=_pde)
    fl = res["flow"]
    assert fl["n_empty_visit2"] == 1 and fl["n_qc_failed"] == 1 and fl["n_eligible"] == 24 and fl["stopped"]
    assert res["status"].startswith("stopped") and loads == [] and "dice" not in res


def test_full_synthetic_run_scores_all_arms_on_common_patients(monkeypatch):
    loads = []
    cases = _cases(loads=loads)
    cases[0] = replace(cases[0], visit2=np.zeros(SHAPE, bool))
    res = lr.run_lumiere_cv(cases, _folds(cases), pde_only_forecast=_pde, candidate_forecast=_pde,
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
    assert res["dice"]["copy_last"] == res["dice"]["pde_only"]  # both are the visit-2 core here


def test_film_reuse_only_on_full_fingerprint_match():
    cases = _cases()
    cache = {}
    first = lr.run_lumiere_cv(cases, _folds(cases), pde_only_forecast=_pde, candidate_forecast=_pde,
                              recipe=TINY, film_reuse_cache=cache)
    assert not any(f["film_reused"] for f in first["folds"].values())
    second = lr.run_lumiere_cv(cases, _folds(cases), pde_only_forecast=_pde, candidate_forecast=_pde,
                               recipe=TINY, film_reuse_cache=cache)
    assert all(f["film_reused"] for f in second["folds"].values())
    assert second["dice"]["film_unet"] == first["dice"]["film_unet"]
    shifted = [replace(c, history_weeks=c.history_weeks + 1) for c in cases]
    third = lr.run_lumiere_cv(shifted, _folds(shifted), pde_only_forecast=_pde, candidate_forecast=_pde,
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
