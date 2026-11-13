"""Synthetic tests for the frozen v2 analysis. Made-up numbers; no policy outcome."""
import numpy as np

import trackc_v2_analysis as an


def diffs(mean_s, mean_m, sd=20.0, n=400, seed=0):
    rng = np.random.default_rng(seed)
    return {"single": rng.normal(mean_s, sd, n).tolist(), "mixture": rng.normal(mean_m, sd, n).tolist()}


def test_clear_gain_in_both_families_is_positive():
    assert an.analyse(diffs(60, 60), primary=True)["status"] == "positive"


def test_gain_in_one_family_only_is_not_positive():
    assert an.analyse(diffs(60, 0), primary=True)["status"] == "not_positive"


def test_clear_loss_is_negative():
    assert an.analyse(diffs(-60, -60), primary=True)["status"] == "negative_both_below_zero"


def test_too_few_worlds_fail_precision():
    out = an.analyse(diffs(60, 60, n=40), primary=True)
    assert out["status"] == "inconclusive_precision_not_met"


def test_sensitivity_cell_is_never_positive():
    assert an.analyse(diffs(60, 60), primary=False)["status"] == "descriptive_both_above_zero"


def test_missing_family_or_bad_value_is_inconclusive():
    d = diffs(60, 60)
    d["mixture"] = [1.0]
    assert an.analyse(d, primary=True)["status"] == "inconclusive"
    d = diffs(60, 60)
    d["single"][0] = 900.0
    assert an.analyse(d, primary=True)["status"] == "inconclusive"


def test_order_is_preserved_not_sorted(monkeypatch):
    seen = []
    monkeypatch.setattr(an, "ge_betting_interval", lambda x, **k: seen.append(list(x)) or (1.0, 2.0))
    d = {"single": [5.0, -3.0, 9.0], "mixture": [2.0, 1.0, 0.5]}
    an.analyse(d, primary=True)
    assert seen[0] == [5.0, -3.0, 9.0]


def _st(s):
    return {"status": s}


def test_integration_positive_only_if_all_three_positive():
    c = an.INTEGRATION_CONTRASTS
    assert an.integration_status({k: _st("positive") for k in c}) == "positive"
    mixed = {k: _st("positive") for k in c}
    mixed[c[1]] = _st("not_positive")
    assert an.integration_status(mixed) == "not_positive"
    mixed[c[2]] = _st("inconclusive_precision_not_met")
    assert an.integration_status(mixed) == "inconclusive"


def test_benchmark_names_integration_contrasts():
    from pathlib import Path
    src = (Path(an.__file__).parent / "211_trackc_policy_benchmark.py").read_text(encoding="utf-8")
    for k in an.INTEGRATION_CONTRASTS:
        assert f'"{k}"' in src
