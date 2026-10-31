"""Synthetic tests for the Track B fixed-growth comparator. No patient data."""
import numpy as np
import pytest

from src.fk_fixed_growth import (
    RADIUS_GRID_MM,
    FixedGrowthCase,
    fixed_growth_forecast,
    run_outer_fixed_growth,
    select_radius,
)

SHAPE = (31, 31, 31)
CENTER = (15, 15, 15)


def _ball(radius_mm, spacing=(2.0, 2.0, 2.0), shape=SHAPE):
    grid = np.indices(shape).astype(float)
    d2 = sum(((grid[i] - CENTER[i]) * spacing[i]) ** 2 for i in range(3))
    return d2 <= radius_mm ** 2


def _case(pid, grow_mm, *, spacing=(2.0, 2.0, 2.0), target=True):
    last = _ball(6.0, spacing)
    support = np.ones(SHAPE, bool)
    return FixedGrowthCase(pid, last, support, spacing, _ball(6.0 + grow_mm, spacing) if target else None)


def test_radius_zero_is_copy_last_clipped_to_support():
    last = _ball(6.0)
    support = np.ones(SHAPE, bool)
    support[:, :, 16:] = False
    out = fixed_growth_forecast(last, support, spacing_mm=(2, 2, 2), radius_mm=0)
    assert np.array_equal(out, last & support)


def test_support_clipping_removes_voxels_outside_support():
    support = np.ones(SHAPE, bool)
    support[20:] = False
    out = fixed_growth_forecast(_ball(6.0), support, spacing_mm=(2, 2, 2), radius_mm=10)
    assert out.any() and not out[20:].any()


def test_physical_spacing_sets_growth_per_axis():
    last = np.zeros(SHAPE, bool)
    last[CENTER] = True
    out = fixed_growth_forecast(last, np.ones(SHAPE, bool), spacing_mm=(1.0, 2.0, 4.0), radius_mm=4.0)
    x, y, z = CENTER
    assert out[x + 4, y, z] and not out[x + 5, y, z]
    assert out[x, y + 2, z] and not out[x, y + 3, z]
    assert out[x, y, z + 1] and not out[x, y, z + 2]


def test_empty_history_and_negative_radius_rejected():
    with pytest.raises(ValueError, match="empty last history"):
        fixed_growth_forecast(np.zeros(SHAPE, bool), np.ones(SHAPE, bool), spacing_mm=(2, 2, 2), radius_mm=2)
    with pytest.raises(ValueError, match="nonnegative"):
        fixed_growth_forecast(_ball(6.0), np.ones(SHAPE, bool), spacing_mm=(2, 2, 2), radius_mm=-2)


def test_empty_target_scores_with_empty_rule():
    case = FixedGrowthCase("p", _ball(6.0), np.ones(SHAPE, bool), (2, 2, 2), np.zeros(SHAPE, bool))
    sel = select_radius([case])
    assert sel["radius_mm"] == 0.0 and sel["training_mean_dice"][0] == 0.0


def test_selection_recovers_true_growth():
    sel = select_radius([_case(f"p{i}", 4.0) for i in range(4)])
    assert sel["radius_mm"] == 4.0 and not sel["boundary_hit"]


def test_tie_goes_to_smaller_radius():
    target = _ball(6.0)
    case = FixedGrowthCase("p", _ball(6.0), target.copy(), (2, 2, 2), target)
    # Support equals the last core, so every radius gives the same forecast and the same Dice.
    sel = select_radius([case])
    assert len(set(np.round(sel["training_mean_dice"], 12))) == 1 and sel["radius_mm"] == 0.0


def test_boundary_hit_is_flagged():
    sel = select_radius([_case("p", 30.0)])
    assert sel["radius_mm"] == RADIUS_GRID_MM[-1] and sel["boundary_hit"]


def test_outer_test_targets_do_not_affect_selection_or_forecast():
    cases = [_case(f"p{i}", 4.0) for i in range(6)]
    folds = {f"p{i}": i % 3 for i in range(6)}
    base = run_outer_fixed_growth(cases, folds)
    # Replace fold-0 test targets with very different masks.
    changed = [FixedGrowthCase(c.patient_id, c.last_mask, c.support, c.spacing_mm,
                               _ball(20.0) if folds[c.patient_id] == 0 else c.target_mask) for c in cases]
    alt = run_outer_fixed_growth(changed, folds)
    assert base["selections"][0]["radius_mm"] == alt["selections"][0]["radius_mm"]
    for a, b in zip(base["rows"], alt["rows"]):
        if a["fold"] == 0:
            assert np.array_equal(a["forecast"], b["forecast"])
    # Other folds train on fold 0, so their selection may change; fold 0 must not.
    assert base["selections"][0]["n_training"] == 4


def test_training_case_without_target_is_rejected():
    cases = [_case(f"p{i}", 4.0) for i in range(4)]
    cases[0] = FixedGrowthCase("p0", cases[0].last_mask, cases[0].support, (2, 2, 2), None)
    with pytest.raises(ValueError, match="no target mask"):
        run_outer_fixed_growth(cases, {f"p{i}": int(i == 0) for i in range(4)})


def test_forecast_needs_no_target():
    sel = select_radius([_case(f"p{i}", 4.0) for i in range(3)])
    out = fixed_growth_forecast(_ball(6.0), np.ones(SHAPE, bool), spacing_mm=(2, 2, 2),
                                radius_mm=sel["radius_mm"])
    target = _ball(10.0)
    assert sel["radius_mm"] == 4.0
    assert 2 * (out & target).sum() / (out.sum() + target.sum()) > 0.95


def test_target_outside_support_counts_as_miss():
    support = np.ones(SHAPE, bool)
    support[CENTER[0] + 4:] = False
    target = _ball(10.0)
    case = FixedGrowthCase("p", _ball(6.0), support, (2, 2, 2), target)
    out = run_outer_fixed_growth([case, _case("q", 4.0)], {"p": 0, "q": 1})
    row = next(r for r in out["rows"] if r["patient_id"] == "p")
    assert row["target_voxels_outside_support"] == int((target & ~support).sum()) > 0
    assert row["dice"] < 1.0


# ---- End-to-end wiring through the patient-fold CV wrapper (synthetic only) ----
from src.fk_fixed_growth import add_fixed_growth_arm  # noqa: E402
from src.fk_mu_scoring import make_first_third_score_callback  # noqa: E402
from src.fk_patient_cv import ForecastFailure, PatientSeries, run_first_third_90day_patient_cv  # noqa: E402


def _cohort(n=10, grow_mm=4.0, empty_history_ids=()):
    records, folds, domains = [], {}, {}
    for i in range(n):
        pid = f"s{i:02d}"
        last = np.zeros(SHAPE, bool) if pid in empty_history_ids else _ball(6.0)
        records.append(PatientSeries(pid, [0.0, 30.0, 30.0 + 76 + (i % 29)],
                                     [_ball(4.0), last, _ball(6.0 + grow_mm)]))
        folds[pid] = i % 5
        domains[pid] = np.ones(SHAPE, bool)
    return records, folds, domains


def _run(records, folds, domains, seen_forecasts=None):
    def fit(training_records):
        return {"n": len(training_records)}

    def forecast(history, _nuisance):
        if seen_forecasts is not None:
            seen_forecasts.append((history.patient_id, history.forecast_day, len(history.masks)))
        last = np.asarray(history.masks[-1], bool)
        if not last.any():
            raise ForecastFailure("empty history")
        return {"fallback_masks": {"hybrid": last, "no_pde": last, "persistence": last}, "history_mask": last}

    def fallback(history, _nuisance):
        last = np.asarray(history.masks[-1], bool)
        return {"fallback_masks": {"hybrid": last, "no_pde": last, "persistence": last}, "history_mask": last}

    base_score = make_first_third_score_callback(spacing_mm=2.0, valid_support_by_patient=domains)
    fit_w, fc_w, fb_w, sc_w = add_fixed_growth_arm(fit, forecast, fallback, base_score, domain_by_patient=domains,
                                                   spacing_mm=(2.0, 2.0, 2.0), valid_support_by_patient=domains)
    return run_first_third_90day_patient_cv(records, folds, fit_w, fc_w, sc_w, fallback_forecast=fb_w, n_folds=5)


def test_runner_wiring_scores_fixed_growth_per_patient():
    cv = _run(*_cohort())
    assert len(cv["rows"]) == 10
    for row in cv["rows"]:
        assert row["scores"]["fixed_growth_radius_mm"] == 4.0
        assert row["scores"]["dice_fixed_growth"] > row["scores"]["dice_persistence"]
        assert "paired_delta_dice_hybrid_minus_fixed_growth" in row["scores"]


def test_runner_horizon_is_nominal_day_90_and_target_day_only_in_scoring():
    seen = []
    records, folds, domains = _cohort()
    cv = _run(records, folds, domains, seen_forecasts=seen)
    for pid, forecast_day, n_masks in seen:
        assert forecast_day == 30.0 + 90.0 and n_masks == 2
    gaps = sorted(row["target_day"] - 30.0 for row in cv["rows"])
    assert gaps[0] >= 76 and gaps[-1] <= 104


def test_runner_outer_test_target_does_not_change_radius_or_forecast():
    records, folds, domains = _cohort()
    base = _run(records, folds, domains)
    # Change only fold-0 test targets to a very different mask.
    changed = [PatientSeries(r.patient_id, r.days, [r.masks[0], r.masks[1], _ball(20.0)])
               if folds[r.patient_id] == 0 else r for r in records]
    alt = _run(changed, folds, domains)
    for a, b in zip(base["rows"], alt["rows"]):
        if a["fold"] == 0:
            assert a["scores"]["fixed_growth_radius_mm"] == b["scores"]["fixed_growth_radius_mm"]
            assert a["scores"]["fixed_growth_status"] == b["scores"]["fixed_growth_status"]


def test_runner_patient_grouped_folds_and_empty_history():
    records, folds, domains = _cohort(empty_history_ids=("s03",))
    cv = _run(records, folds, domains)
    row = next(r for r in cv["rows"] if r["patient_id"] == "s03")
    assert row["forecast_status"] == "fallback"
    assert row["scores"]["fixed_growth_status"] == "empty_history_copy_last"
    assert row["scores"]["dice_fixed_growth"] == 0.0
    assert {r["patient_id"] for r in cv["rows"]} == set(folds)
    assert all(a["n_test_patients"] == 2 and a["n_training_patients"] == 8 for a in cv["fold_audit"])


def test_runner_support_clipping_and_out_of_support_rule():
    records, folds, domains = _cohort()
    clipped = {pid: d.copy() for pid, d in domains.items()}
    for d in clipped.values():
        d[CENTER[0] + 6:] = False  # removes part of the 10 mm target too
    trimmed = [PatientSeries(r.patient_id, r.days, [m & clipped[r.patient_id] for m in r.masks]) for r in records]
    cv = _run(trimmed, folds, clipped)
    for row in cv["rows"]:
        assert row["scores"]["dice_fixed_growth"] > 0.9
    # A support that excludes forecast voxels but not the target must fail closed in scoring.
    tight = {pid: (_ball(10.0) | _ball(6.0)) for pid in domains}
    wide_domains = domains
    from src.fk_fixed_growth import add_fixed_growth_arm as arm
    fit_w, fc_w, fb_w, sc_w = arm(lambda t: {}, lambda h, n: {"fallback_masks": {"hybrid": h.masks[-1],
                                  "no_pde": h.masks[-1]}, "history_mask": h.masks[-1]},
                                  lambda h, n: h.masks[-1], lambda p, t: {"dice_hybrid": 1.0},
                                  domain_by_patient=wide_domains, spacing_mm=(2, 2, 2),
                                  valid_support_by_patient={pid: _ball(6.0) for pid in tight})
    with pytest.raises(ValueError, match="outside validated target support"):
        run_first_third_90day_patient_cv(records, folds, fit_w, fc_w, sc_w, fallback_forecast=fb_w, n_folds=5)
