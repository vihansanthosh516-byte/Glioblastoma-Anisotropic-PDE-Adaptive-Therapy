"""Synthetic tests: unobserved voxels (outside geometric coverage) cannot affect fits, selection or scores.

Supported labels stay fixed in every test; only labels or predictions outside support change. No patient data.
"""
from dataclasses import replace

import numpy as np
import pytest
import torch

from src import fk_time_aware_baselines as ta
from tests.test_fk_lumiere_coverage_scoring import _cov, _run, _with_outside_lesion
from tests.test_fk_lumiere_runner import SHAPE, TINY, _ball, _cases

RNG = np.random.default_rng(7)


def _valid():
    v = np.ones(SHAPE, bool)
    v[12:] = False                      # slices 12-15 unobserved
    return v


def _noise_outside(a, v, seed):
    out = np.asarray(a).copy()
    r = np.random.default_rng(seed).random(out.shape)
    out[~v] = (r[~v] > 0.5) if out.dtype == bool else r[~v] * 20 - 10
    return out


# ------------------------------------------------------------------ masked loss
def test_masked_loss_ignores_labels_and_predictions_outside_support():
    v = _valid()
    target = _ball(4).astype(np.float32)
    logits = RNG.normal(size=SHAPE).astype(np.float32)
    w = torch.from_numpy(v.astype(np.float32))
    base = ta.masked_loss(torch.from_numpy(logits), torch.from_numpy(target), w, 1.0)
    t2 = _noise_outside(target, v, 1).astype(np.float32)
    l2 = _noise_outside(logits, v, 2).astype(np.float32)
    assert np.array_equal(t2[v], target[v]) and np.array_equal(l2[v], logits[v])
    again = ta.masked_loss(torch.from_numpy(l2), torch.from_numpy(t2), w, 1.0)
    assert float(again) == pytest.approx(float(base), abs=1e-6)


def test_masked_loss_all_ones_equals_frozen_loss():
    target = torch.from_numpy(_ball(4).astype(np.float32))
    logits = torch.from_numpy(RNG.normal(size=SHAPE).astype(np.float32))
    ones = torch.ones(SHAPE)
    from src.fk_learned_baselines import _loss
    frozen = _loss(logits, target, 1.0)
    assert float(ta.masked_loss(logits, target, ones, 1.0)) == pytest.approx(float(frozen), rel=1e-5)


# ------------------------------------------------------------------ growth-rate selection
def _timed(targets, v):
    return [ta.TimedGrowthCase(f"p{i}", _ball(3), np.ones(SHAPE, bool), (2.0, 2.0, 2.0), 8, t, v)
            for i, t in enumerate(targets)]


def test_selected_growth_rate_ignores_labels_outside_support():
    v = _valid()
    targets = [_ball(r) for r in (3, 4, 5)]
    base = ta.select_rate(_timed(targets, v))
    altered = [_noise_outside(t, v, i) for i, t in enumerate(targets)]
    for a, t in zip(altered, targets):
        assert np.array_equal(a[v], t[v]) and not np.array_equal(a, t)
    again = ta.select_rate(_timed(altered, v))
    assert again["rate_mm_per_day"] == base["rate_mm_per_day"]
    assert again["training_mean_dice"] == base["training_mean_dice"]


def test_selection_scores_ignore_forecast_outside_support(monkeypatch):
    v = _valid()
    targets = [_ball(r) for r in (3, 4, 5)]
    base = ta.select_rate(_timed(targets, v))
    real = ta.fixed_growth_forecast

    def noisy(*a, **k):
        return _noise_outside(real(*a, **k), v, 11)
    monkeypatch.setattr(ta, "fixed_growth_forecast", noisy)
    again = ta.select_rate(_timed(targets, v))
    assert again["training_mean_dice"] == base["training_mean_dice"]


# ------------------------------------------------------------------ FiLM weights
def _film(labels, v):
    x = np.stack([_ball(3), _ball(3), np.ones(SHAPE)]).astype(np.float32)
    pairs = [(x, y) for y in labels]
    return ta.train_film_models(pairs, [(8, 8)] * len(pairs), fold=0, recipe=TINY, valid_masks=[v] * len(pairs))


def test_film_weights_ignore_labels_outside_support():
    v = _valid()
    labels = [_ball(4), _ball(5)]
    a = _film(labels, v)
    b = _film([_noise_outside(y, v, i + 20) for i, y in enumerate(labels)], v)
    assert a["status"] == b["status"]
    for ma, mb in zip(a["models"], b["models"]):
        for (ka, ta_), (kb, tb) in zip(ma.state_dict().items(), mb.state_dict().items()):
            assert ka == kb and torch.equal(ta_, tb), ka


def test_film_fingerprint_covers_validity_masks():
    x = np.zeros((3, 4, 4, 4), np.float32)
    y = np.zeros((4, 4, 4), bool)
    v1, v2 = np.ones((4, 4, 4), bool), np.ones((4, 4, 4), bool)
    v2[0] = False
    fp = lambda vm: ta.film_training_fingerprint(["a"], [(x, y)], [(13, 13)], fold=0, recipe=TINY,
                                                 environment={"device": "cpu"}, valid_masks=vm)
    assert fp([v1]) != fp([v2]) and fp([v1]) != fp(None)


# ------------------------------------------------------------------ runner end to end
def test_runner_partial_coverage_fails_closed_for_pde_and_count():
    cases = [replace(c, coverage=_cov()) for c in _cases(12)]
    res = _run(cases, require_coverage=True)
    assert res["coverage_scoring"]["partial_coverage"] is True
    for arm in ("candidate", "pde_only", "matched_volume_no_pde_ablation"):
        assert all(d is None for d in res["dice"][arm].values()), arm
    assert res["count_arm"]["status"] == "not_scored_partial_coverage"
    assert all(r == {"status": "not_scored_partial_coverage"} for r in res["count_arm"]["rows"].values())
    for arm in ("copy_last", "time_aware_fixed_growth"):     # tiny FiLM may be unavailable (None) by design
        assert all(d is not None for d in res["dice"][arm].values()), arm


def test_runner_outside_labels_change_no_fit_selection_or_score():
    cases = [replace(c, coverage=_cov()) for c in _cases(12)]
    base = _run(cases, require_coverage=True)
    alt = _run(_with_outside_lesion(cases), require_coverage=True)
    assert base["dice"] == alt["dice"]
    assert base["count_arm"] == alt["count_arm"]
    for f in base["folds"]:
        assert repr(base["folds"][f]) == repr(alt["folds"][f]), f


def test_runner_full_coverage_runs_every_arm():
    cases = [replace(c, coverage=_cov(shift_target=0)) for c in _cases(12)]
    res = _run(cases, require_coverage=True)
    assert res["coverage_scoring"]["partial_coverage"] is False
    assert "mean_poisson_deviance" in res["count_arm"]
    assert any(d is not None for d in res["dice"]["pde_only"].values())


# ------------------------------------------------------------------ padding and status (Codex review 2)
def test_full_coverage_film_training_is_exactly_the_frozen_padded_path():
    shape = (15, 13, 11)                                     # needs padding for TINY (2 levels)
    rng = np.random.default_rng(3)
    x = rng.random((3,) + shape).astype(np.float32)
    ys = [rng.random(shape) > 0.7, rng.random(shape) > 0.6]
    pairs = [(x, y) for y in ys]
    frozen = ta.train_film_models(pairs, [(8, 8)] * 2, fold=0, recipe=TINY)
    masked = ta.train_film_models(pairs, [(8, 8)] * 2, fold=0, recipe=TINY,
                                  valid_masks=[np.ones(shape, bool)] * 2)
    assert frozen["status"] == masked["status"]
    for ma, mb in zip(frozen["models"], masked["models"]):
        for (ka, a), (kb, b) in zip(ma.state_dict().items(), mb.state_dict().items()):
            assert ka == kb and torch.equal(a, b), ka


def test_partial_film_padding_voxels_stay_background_and_unobserved_are_ignored():
    shape = (15, 13, 11)
    rng = np.random.default_rng(4)
    x = rng.random((3,) + shape).astype(np.float32)
    v = np.ones(shape, bool)
    v[12:] = False
    y = rng.random(shape) > 0.7
    a = ta.train_film_models([(x, y)], [(8, 8)], fold=0, recipe=TINY, valid_masks=[v])
    y2 = y.copy()
    y2[~v] = ~y2[~v]
    b = ta.train_film_models([(x, y2)], [(8, 8)], fold=0, recipe=TINY, valid_masks=[v])
    for ma, mb in zip(a["models"], b["models"]):
        for (k, t1), (_, t2) in zip(ma.state_dict().items(), mb.state_dict().items()):
            assert torch.equal(t1, t2), k


class _NoCallProvider:
    def fit_fold(self, view):
        raise AssertionError("fit_fold called under partial coverage")

    def forecast(self, case, nuisance):
        raise AssertionError("forecast called under partial coverage")


def test_partial_coverage_never_calls_pde_provider_and_is_incomplete():
    from src import fk_lumiere_runner as lr
    from tests.test_fk_lumiere_runner import _folds
    cases = [replace(c, coverage=_cov()) for c in _cases(12)]
    res = lr.run_lumiere_cv(cases, _folds(cases), pde_provider=_NoCallProvider(), recipe=TINY,
                            minimum=10, n_bootstrap=50, require_coverage=True)
    assert res["status"].startswith("incomplete") and "scored" != res["status"]
    assert res["comparison_complete"] is False
    assert res["count_arm"]["status"] == "not_scored_partial_coverage" and res["count_arm"]["reason"]
    for fold in res["folds"].values():
        assert fold["pde_only_cutoff"]["cutoff"] is None and fold["pde_only_cutoff"]["reason"]


def test_full_coverage_status_is_scored():
    res = _run([replace(c, coverage=_cov(shift_target=0)) for c in _cases(12)], require_coverage=True)
    assert res["status"] == "scored (synthetic)"
