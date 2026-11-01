"""Synthetic tests for the Track B fold-trained U-Net and Poisson count baselines. No patient data."""
import dataclasses

import numpy as np
import pytest

from src import fk_learned_baselines as lb
from src.fk_mu_scoring import make_first_third_score_callback
from src.fk_patient_cv import ForecastFailure, PatientSeries, run_first_third_90day_patient_cv

SHAPE = (12, 12, 8)
FAST = dataclasses.replace(lb.UNET_RECIPE, widths=(4, 8, 8, 8), epochs=2)


def _ball(r, c=(6, 6, 4), shape=SHAPE, spacing=2.0):
    g = np.indices(shape).astype(float)
    return sum(((g[i] - c[i]) * spacing) ** 2 for i in range(3)) <= r ** 2


def test_frozen_recipe_values():
    r = lb.UNET_RECIPE
    assert r.widths == (16, 32, 64, 128) and r.epochs == 300 and r.learning_rate == 1e-3
    assert r.betas == (0.9, 0.999) and r.weight_decay == 0.0 and r.n_seeds == 3
    assert r.seed_offsets == (0, 1000, 2000) and r.base_seed == 20261102 and r.threshold == 0.5


def test_runner_uses_frozen_recipe():
    import inspect
    from src import fk_mu_patient_runner as runner
    assert runner.UNET_RECIPE is lb.UNET_RECIPE
    assert "recipe=UNET_RECIPE" in inspect.getsource(runner.run_mu_first_two_to_third_development_cv)


def test_unet_shapes_four_levels_and_padding():
    import torch
    model = lb.UNet3D(lb.UNET_RECIPE.widths)
    assert len(model.down) == 4 and len(model.up) == 3
    assert lb._padded_shape((13, 9, 7), 4) == (16, 16, 8)
    out = lb.UNet3D((4, 8, 8, 8))(torch.zeros(1, 3, 16, 16, 8))
    assert out.shape == (1, 1, 16, 16, 8)


def test_inputs_are_history_only():
    x = lb.unet_inputs(_ball(4), _ball(6), np.ones(SHAPE, bool))
    assert x.shape == (3,) + SHAPE and set(np.unique(x)) <= {0.0, 1.0}


def test_three_seeds_are_averaged_and_thresholded():
    pairs = [(lb.unet_inputs(_ball(4), _ball(6), np.ones(SHAPE, bool)), _ball(8)) for _ in range(2)]
    trained = lb.train_unet_models(pairs, fold=0, recipe=FAST)
    if trained["status"] != "ok":
        pytest.skip("tiny synthetic U-Net produced an empty output; unavailable path covered elsewhere")
    assert trained["seeds"] == [20261102, 20262102, 20263102]
    domain = np.ones(SHAPE, bool)
    domain[:, :, 6:] = False
    pred = lb.predict_unet(trained, pairs[0][0], domain=domain)
    assert len(pred["seed_probabilities"]) == 3
    assert np.allclose(pred["probability"], np.mean(pred["seed_probabilities"], axis=0))
    assert np.array_equal(pred["mask"], (pred["probability"] >= 0.5) & domain)


def test_training_is_seed_deterministic():
    pairs = [(lb.unet_inputs(_ball(4), _ball(6), np.ones(SHAPE, bool)), _ball(8))]
    one = dataclasses.replace(FAST, n_seeds=1)
    a = lb.train_unet_models(pairs, fold=1, recipe=one)
    b = lb.train_unet_models(pairs, fold=1, recipe=one)
    if a["status"] == "ok":
        pa = lb.predict_unet(a, pairs[0][0], domain=np.ones(SHAPE, bool))["probability"]
        pb = lb.predict_unet(b, pairs[0][0], domain=np.ones(SHAPE, bool))["probability"]
        assert np.array_equal(pa, pb)
    else:
        assert a["reason"] == b["reason"]


def test_unavailable_without_training_pairs():
    out = lb.train_unet_models([], fold=0, recipe=FAST)
    assert out["status"] == "unavailable"
    with pytest.raises(ValueError, match="unavailable"):
        lb.predict_unet(out, None, domain=np.ones(SHAPE, bool))


def test_new_label_pieces_definition():
    prev = _ball(4)
    cur = prev.copy()
    cur[0:2, 0:2, 0:2] = True          # 8-voxel far piece
    cur[11, 11, 7] = True              # 1-voxel piece, below minimum
    cur[0:2, 10:12, 6:8] = True        # second far piece
    assert lb.new_label_pieces(prev, cur, spacing_mm=(2, 2, 2)) == 2
    grown_edge = _ball(5.5)            # within 2 mm of prev -> not new
    assert lb.new_label_pieces(prev, grown_edge, spacing_mm=(2, 2, 2)) == 0
    dom = np.ones(SHAPE, bool)
    dom[0:2] = False
    assert lb.new_label_pieces(prev, cur, spacing_mm=(2, 2, 2), domain=dom) == 0
    assert lb.new_label_pieces(np.zeros(SHAPE, bool), _ball(4), spacing_mm=(2, 2, 2)) == 1


def test_count_model_recovers_rate_and_falls_back():
    rng = np.random.default_rng(0)
    lv = rng.uniform(6, 9, 400)
    t = rng.uniform(76, 104, 400)
    y = rng.poisson(np.exp(-6.0 + 0.4 * lv) * t)
    m = lb.fit_count_model(y, lv, t)
    assert m["status"] == "glm"
    assert abs(m["coef"][0] + 6.0) < 0.6 and abs(m["coef"][1] - 0.4) < 0.08
    fb = lb.fit_count_model([1, 0], [7.0, 8.0], [90, 90])
    assert fb["status"] == "mean_rate_fallback" and fb["mean_rate_per_day"] == pytest.approx(1 / 180)
    assert lb.count_rate_per_day(m, np.nan) == m["mean_rate_per_day"]


def test_count_metrics_values():
    m = lb.count_metrics(2.0, 0)
    assert m["poisson_deviance"] == pytest.approx(4.0) and m["abs_error"] == 2.0
    assert m["brier_any_new"] == pytest.approx((1 - np.exp(-2.0)) ** 2)
    m = lb.count_metrics(1.0, 1)
    assert m["poisson_deviance"] == pytest.approx(0.0)
    assert m["brier_any_new"] == pytest.approx(np.exp(-1.0) ** 2)
    with pytest.raises(ValueError):
        lb.count_metrics(0.0, 1)


def _cohort(n=10, empty=()):
    records, folds, domains = [], {}, {}
    for i in range(n):
        pid = f"s{i:02d}"
        v2 = np.zeros(SHAPE, bool) if pid in empty else _ball(6)
        v3 = _ball(8)
        if i % 2:
            v3 = v3.copy()
            v3[0:2, 0:2, 0:2] = True
        records.append(PatientSeries(pid, [0.0, 30.0, 30.0 + 76 + i], [_ball(4), v2, v3]))
        folds[pid] = i % 5
        domains[pid] = np.ones(SHAPE, bool)
    return records, folds, domains


def _run(records, folds, domains, recipe=FAST, seen=None):
    def fit(training):
        if seen is not None:
            seen.append(sorted(r.patient_id for r in training))
        return {}

    def fc(h, _n):
        last = np.asarray(h.masks[-1], bool)
        if not last.any():
            raise ForecastFailure("empty")
        return {"fallback_masks": {"hybrid": last, "no_pde": last, "persistence": last}, "history_mask": last}

    def fb(h, _n):
        last = np.asarray(h.masks[-1], bool)
        return {"fallback_masks": {"hybrid": last, "no_pde": last, "persistence": last}, "history_mask": last}

    arms = lb.add_learned_baselines_arm(fit, fc, fb, make_first_third_score_callback(2.0, valid_support_by_patient=domains),
                                        domain_by_patient=domains, spacing_mm=(2, 2, 2),
                                        valid_support_by_patient=domains, recipe=recipe)
    f, c, b, s = arms
    return run_first_third_90day_patient_cv(records, folds, f, c, s, fallback_forecast=b, n_folds=5)


def test_runner_fold_isolation_and_count_scoring():
    seen = []
    records, folds, domains = _cohort(empty=("s03",))
    cv = _run(records, folds, domains, seen=seen)
    for fold, train_ids in enumerate(seen):
        assert not set(train_ids) & {p for p, f in folds.items() if f == fold}
    for row in cv["rows"]:
        sc = row["scores"]
        assert sc["count_interval_days"] == row["target_day"] - 30.0
        expected = 1 if int(row["patient_id"][1:]) % 2 else 0
        if row["patient_id"] == "s03":
            expected += 1  # empty visit-2 core: the whole visit-3 core is a new piece
        assert sc["count_observed_new_pieces"] == expected
        assert "count_poisson_deviance" in sc and "count_mean_rate_brier_any_new" in sc
        assert not any(k.startswith("dice_count") for k in sc)
    s03 = next(r for r in cv["rows"] if r["patient_id"] == "s03")["scores"]
    assert s03["unet_status"] == "empty_history_copy_last" and s03["dice_unet"] == 0.0


def test_changing_test_target_does_not_change_forecasts():
    records, folds, domains = _cohort()
    base = _run(records, folds, domains)
    changed = [PatientSeries(r.patient_id, r.days, [r.masks[0], r.masks[1], _ball(2)])
               if folds[r.patient_id] == 0 else r for r in records]
    alt = _run(changed, folds, domains)
    for a, b in zip(base["rows"], alt["rows"]):
        if a["fold"] == 0:
            assert a["scores"]["count_predicted"] == b["scores"]["count_predicted"]
            assert a["scores"]["unet_status"] == b["scores"]["unet_status"]


def test_forecast_mask_identical_when_test_target_changes():
    records, folds, domains = _cohort()
    captured = {}

    def capture(records_):
        def fit(training):
            return {}

        def fc(h, _n):
            last = np.asarray(h.masks[-1], bool)
            return {"fallback_masks": {"hybrid": last, "no_pde": last, "persistence": last}, "history_mask": last}

        f, c, b, s = lb.add_learned_baselines_arm(fit, fc, fc, lambda p, t: {"dice_hybrid": 0.0},
                                                  domain_by_patient=domains, spacing_mm=(2, 2, 2), recipe=FAST)

        def c2(h, n):
            out = c(h, n)
            captured.setdefault(id(records_), {})[h.patient_id] = out["_unet_mask"]
            return out
        run_first_third_90day_patient_cv(records_, folds, f, c2, s, fallback_forecast=b, n_folds=5)
        return captured[id(records_)]

    a = capture(records)
    changed = [PatientSeries(r.patient_id, r.days, [r.masks[0], r.masks[1], _ball(2)])
               if folds[r.patient_id] == 0 else r for r in records]
    b = capture(changed)
    for pid, f in folds.items():
        if f == 0:
            assert (a[pid] is None and b[pid] is None) or np.array_equal(a[pid], b[pid])


# ---- Proposed (not adopted) history-centred crop, and U-Net reuse check ----
CROP = dataclasses.replace(FAST, crop_voxels=16)
BIG = (32, 32, 20)


def test_frozen_recipe_has_no_crop():
    assert lb.UNET_RECIPE.crop_voxels is None


def test_crop_box_uses_history_only_and_stays_inside_grid():
    v2 = _ball(4, c=(16, 16, 10), shape=BIG)
    box = lb.history_crop_box(v2, BIG, 16)
    assert [s.stop - s.start for s in box] == [16, 16, 16] and box[0] == slice(8, 24)
    edge = _ball(2, c=(1, 30, 19), shape=BIG)
    box = lb.history_crop_box(edge, BIG, 16)
    assert box[0] == slice(0, 16) and box[1] == slice(16, 32) and box[2] == slice(4, 20)
    assert lb.history_crop_box(v2, (32, 32, 5), 16)[2] == slice(0, 5)
    with pytest.raises(ValueError):
        lb.history_crop_box(np.zeros(BIG, bool), BIG, 16)


def test_uncrop_maps_back_to_full_grid():
    box = (slice(2, 6), slice(3, 7), slice(1, 5))
    small = np.zeros((4, 4, 4), bool)
    small[0, 0, 0] = small[3, 3, 3] = True
    full = lb.uncrop_mask(small, box, BIG)
    assert full.shape == BIG and full.sum() == 2 and full[2, 3, 1] and full[5, 6, 4]


def _big_cohort(n=10, far_piece=True):
    records, folds, domains = [], {}, {}
    for i in range(n):
        pid = f"b{i:02d}"
        v3 = _ball(8, c=(16, 16, 10), shape=BIG)
        if far_piece:
            v3 = v3.copy()
            v3[0:2, 0:2, 0:2] = True  # outside any 16-voxel box centred at the core
        records.append(PatientSeries(pid, [0.0, 30.0, 30.0 + 80 + i],
                                     [_ball(4, c=(16, 16, 10), shape=BIG), _ball(6, c=(16, 16, 10), shape=BIG), v3]))
        folds[pid] = i % 5
        domains[pid] = np.ones(BIG, bool)
    return records, folds, domains


def _capture_run(records, folds, domains, recipe, base_tag):
    captured = {}

    def fit(training):
        return {"tag": base_tag}  # stands in for the prior-bound-dependent PDE nuisance

    def fc(h, _n):
        last = np.asarray(h.masks[-1], bool)
        return {"fallback_masks": {"hybrid": last, "no_pde": last, "persistence": last}, "history_mask": last}

    f, c, b, s = lb.add_learned_baselines_arm(fit, fc, fc, lambda p, t: {"dice_hybrid": 0.0},
                                              domain_by_patient=domains, spacing_mm=(2, 2, 2), recipe=recipe)

    def c2(h, n):
        out = c(h, n)
        captured[h.patient_id] = (out["_unet_mask"], out["_unet_fingerprint"], out["_unet_crop_box"])
        return out
    cv = run_first_third_90day_patient_cv(records, folds, f, c2, s, fallback_forecast=b, n_folds=5)
    return captured, cv


def test_cropped_prediction_is_zero_outside_box_and_outside_target_is_a_miss():
    records, folds, domains = _big_cohort()
    captured, cv = _capture_run(records, folds, domains, CROP, "primary")
    for row in cv["rows"]:
        assert row["scores"]["unet_target_voxels_outside_crop"] == 8
    for mask, _fp, box in captured.values():
        if mask is not None:
            outside = mask.copy()
            outside[box] = False
            assert not outside.any() and mask.shape == BIG


def test_unet_identical_across_prior_runs_but_not_gap_subset():
    records, folds, domains = _big_cohort(far_piece=False)
    a, _ = _capture_run(records, folds, domains, FAST, "primary_bounds")
    b, _ = _capture_run(records, folds, domains, FAST, "narrower_bounds")
    for pid in a:
        assert a[pid][1] == b[pid][1]
        assert (a[pid][0] is None and b[pid][0] is None) or np.array_equal(a[pid][0], b[pid][0])
    subset = [r for r in records if r.patient_id not in ("b00",)]
    c, _ = _capture_run(subset, {r.patient_id: folds[r.patient_id] for r in subset}, domains, FAST, "gap")
    assert any(c[pid][1] != a[pid][1] for pid in c)


def test_fingerprint_covers_runtime_environment():
    x = np.zeros((3, 4, 4, 4), np.float32)
    y = np.zeros((4, 4, 4), bool)
    env = lb.unet_runtime_environment()
    for key in ("torch_version", "device", "num_threads", "deterministic_algorithms_during_training"):
        assert key in env
    base = lb.unet_training_fingerprint(["p"], [(x, y)], fold=0, recipe=FAST, environment=env)
    assert base == lb.unet_training_fingerprint(["p"], [(x, y)], fold=0, recipe=FAST, environment=dict(env))
    for key, value in (("torch_version", "0.0.0"), ("device", "cuda"), ("num_threads", env["num_threads"] + 1),
                       ("deterministic_algorithms_during_training", False)):
        changed = dict(env, **{key: value})
        assert lb.unet_training_fingerprint(["p"], [(x, y)], fold=0, recipe=FAST, environment=changed) != base


_REAL_TRAIN = lb.train_unet_models


def _reuse_run(records, folds, domains, cache, monkeypatch, calls, env=None):
    real = _REAL_TRAIN

    def counting(pairs, *, fold, recipe):
        calls.append(fold)
        return real(pairs, fold=fold, recipe=recipe)
    monkeypatch.setattr(lb, "train_unet_models", counting)
    if env is not None:
        monkeypatch.setattr(lb, "unet_runtime_environment", lambda: env)
    fc = lambda h, _n: {"fallback_masks": {"hybrid": np.asarray(h.masks[-1], bool)},
                        "history_mask": np.asarray(h.masks[-1], bool)}
    f, c, b, s = lb.add_learned_baselines_arm(lambda t: {}, fc, fc, lambda p, t: {"dice_hybrid": 0.0},
                                              domain_by_patient=domains, spacing_mm=(2, 2, 2), recipe=FAST,
                                              unet_reuse_cache=cache)
    masks = {}

    def c2(h, n):
        out = c(h, n)
        masks[h.patient_id] = out["_unet_mask"]
        return out
    run_first_third_90day_patient_cv(records, folds, f, c2, s, fallback_forecast=b, n_folds=5)
    return masks


def test_unet_reuse_only_on_full_fingerprint_match(monkeypatch):
    records, folds, domains = _big_cohort(far_piece=False)
    cache, calls = {}, []
    a = _reuse_run(records, folds, domains, cache, monkeypatch, calls)
    n_first = len(calls)
    b = _reuse_run(records, folds, domains, cache, monkeypatch, calls)
    assert n_first > 0 and len(calls) == n_first  # second prior run trained nothing
    for pid in a:
        assert (a[pid] is None and b[pid] is None) or np.array_equal(a[pid], b[pid])
    env = dict(lb.unet_runtime_environment(), num_threads=999)
    _reuse_run(records, folds, domains, cache, monkeypatch, calls, env=env)
    assert len(calls) == 2 * n_first  # environment mismatch -> retrain every fold


def test_gap_subset_and_no_cache_always_train(monkeypatch):
    records, folds, domains = _big_cohort(far_piece=False)
    cache, calls = {}, []
    _reuse_run(records, folds, domains, cache, monkeypatch, calls)
    n = len(calls)
    subset = [r for r in records if r.patient_id != "b00"]
    _reuse_run(subset, {r.patient_id: folds[r.patient_id] for r in subset}, domains, cache, monkeypatch, calls)
    assert len(calls) > n
    m = len(calls)
    _reuse_run(records, folds, domains, None, monkeypatch, calls)
    assert len(calls) == m + n


def test_protocol_wrapper_shares_cache_for_prior_runs_only(monkeypatch):
    from src import fk_mu_patient_runner as runner
    seen = []

    def fake_run(cohort, **kw):
        seen.append((kw.get("gap_sensitivity_only", False), kw["unet_reuse_cache"]))
        if len(seen) == 4:
            raise StopIteration
        return {"rows": []}
    monkeypatch.setattr(runner, "run_mu_first_two_to_third_development_cv", fake_run)
    with pytest.raises(StopIteration):
        runner.run_mu_first_two_to_third_with_prior_sensitivity(object())
    prior = [c for gap, c in seen if not gap]
    gap = [c for g, c in seen if g]
    assert len(prior) == 3 and prior[0] is prior[1] is prior[2]
    assert len(gap) == 1 and gap[0] is not prior[0] and gap[0] == {}
    with pytest.raises(ValueError, match="owns U-Net reuse"):
        runner.run_mu_first_two_to_third_with_prior_sensitivity(object(), unet_reuse_cache={})
