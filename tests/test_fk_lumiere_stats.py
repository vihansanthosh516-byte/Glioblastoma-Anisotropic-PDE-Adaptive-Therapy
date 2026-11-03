"""Synthetic tests for the LUMIERE paired bootstrap and max-T pass test. No patient data."""
import numpy as np

from src.fk_lumiere_stats import max_t_pass_test

COMPS = ("a", "b", "c", "d")


def _dice(n, gains, noise=0.02, seed=0):
    rng = np.random.default_rng(seed)
    base = rng.uniform(0.4, 0.6, n)
    ids = [f"p{i}" for i in range(n)]
    out = {"candidate": {p: float(b) for p, b in zip(ids, base)}}
    for k, g in zip(COMPS, gains):
        out[k] = {p: float(b - g + rng.normal(0, noise)) for p, b in zip(ids, base)}
    return out


def test_clear_gain_on_every_contrast_passes():
    r = max_t_pass_test(_dice(28, (0.1, 0.1, 0.1, 0.1)), COMPS, n_bootstrap=2000)
    assert r["status"] == "pass" and all(v["adjusted_interval"][0] > 0 for v in r["contrasts"].values())


def test_one_null_contrast_fails_the_pass():
    r = max_t_pass_test(_dice(28, (0.1, 0.1, 0.1, 0.0)), COMPS, n_bootstrap=2000)
    assert r["status"] == "not passed" and r["contrasts"]["d"]["adjusted_interval"][0] <= 0


def test_missing_arm_is_incomplete_not_pass():
    d = _dice(28, (0.1, 0.1, 0.1, 0.1))
    d["d"]["p3"] = None
    assert max_t_pass_test(d, COMPS)["status"] == "incomplete"
    d2 = _dice(28, (0.1, 0.1, 0.1, 0.1))
    del d2["c"]["p0"]
    assert max_t_pass_test(d2, COMPS)["status"] == "incomplete"


def test_non_finite_scores_fail_closed():
    for bad in (float("nan"), float("inf")):
        d = _dice(28, (0.1, 0.1, 0.1, 0.1))
        d["b"]["p5"] = bad
        r = max_t_pass_test(d, COMPS)
        assert r["status"] == "incomplete" and r["reason"] == "non-finite patient score"


def test_max_t_is_wider_than_unadjusted_and_seeded():
    d = _dice(28, (0.05, 0.05, 0.05, 0.05), noise=0.05)
    r1, r2 = max_t_pass_test(d, COMPS, n_bootstrap=2000), max_t_pass_test(d, COMPS, n_bootstrap=2000)
    assert r1 == r2 and r1["critical_value"] > 1.96 and r1["seed"] == 20261110
