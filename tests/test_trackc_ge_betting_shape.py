"""Empirical check that the GE-betting accept set is an interval matching bisection.

Review finding C-g: the bisection assumes monotone accept sets, which the code does not
prove (bet sizes depend on the wealth path). This synthetic grid check is evidence, not
a proof. No policy outcome is involved.
"""
import numpy as np

import trackc_ge_betting as g


def test_accept_set_is_interval_and_inside_bisection_bounds():
    rng = np.random.default_rng(20261110)
    grid = np.linspace(0.0, 1.0, 801)
    checked = 0
    for n in (7, 30, 207):
        for rep in range(8):
            kind = rep % 4
            if kind == 0:
                x = rng.beta(2, 5, n)
            elif kind == 1:
                x = np.clip(0.5 + 0.02 * rng.standard_normal(n), 0, 1)
            elif kind == 2:
                x = (rng.random(n) < 0.3).astype(float)
            else:
                x = rng.uniform(0.45, 0.55, n)
            interval = g.ge_betting_interval_unit(x)
            if interval is None:
                continue
            xv = g._validated_unit_values(x)
            vp = g._variance_path(xv)
            accepted = np.array([g._accepts_mean(xv, m, 0.025, vp) for m in grid])
            idx = np.flatnonzero(accepted)
            lo, hi = interval
            assert not (accepted & ~((grid >= lo) & (grid <= hi))).any()
            assert idx.size == 0 or not (np.diff(idx) > 1).any()
            checked += 1
    assert checked >= 20
