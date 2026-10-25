"""Track C known-parameter oracle (src/trackc_oracle.py): the grid DP schedule must match an
exhaustive search on short horizons and must never lose to periodic schedules. Synthetic
parameters only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_oracle as to  # noqa: E402
import trackc_sim as tc  # noqa: E402

# fast growth so short horizons are not all capped
FAST = dict(tc.RANGES, r=(0.02, 0.05, "log"))


@pytest.mark.parametrize("n_cycles,budget", [(3, 3), (6, 2), (8, 3), (10, 4)])
def test_dp_matches_brute_force(n_cycles, budget):
    scen = tc.draw_scenarios(6, seed=101 + n_cycles, ranges=FAST)
    for i in range(6):
        p = to._scalar(scen, i)
        bf = to.brute_force(p, n_cycles, budget)
        out = to.best_known(p, n_cycles, budget, G=80)
        assert out["ttp_dp"] >= bf - 7, f"scenario {i}: raw DP {out['ttp_dp']} vs brute force {bf}"
        assert out["ttp"] <= bf + 1e-9
        assert out["ttp"] >= bf - 1e-9, f"scenario {i}: polished {out['ttp']} vs brute force {bf}"
        assert out["schedule"].sum() <= budget


def test_brute_force_is_not_trivial():
    scen = tc.draw_scenarios(6, seed=109, ranges=FAST)
    spreads = []
    for i in range(6):
        p = to._scalar(scen, i)
        rows = [np.isin(np.arange(8), c) for m in range(4) for c in __import__("itertools").combinations(range(8), m)]
        t = to.rollout(p, np.array(rows), 8)
        spreads.append(t.max() - t.min())
    assert max(spreads) > 28, "schedules must differ, or the brute-force test says nothing"


def test_rollout_matches_simulator():
    scen = tc.draw_scenarios(5, seed=110)
    out = tc.simulate(scen, tc.Periodic(4), seed=1)
    for i in range(5):
        p = to._scalar(scen, i)
        t = to.rollout(p, to.periodic_schedules(qs=[4]))[0]
        assert t == pytest.approx(out["ttp"][i])


def test_full_horizon_oracle_not_below_periodic():
    scen = tc.draw_scenarios(4, seed=111)
    for i in range(4):
        p = to._scalar(scen, i)
        out = to.best_known(p, G=60)
        per = to.rollout(p, to.periodic_schedules())
        assert out["ttp_dp"] >= per.max() - 14, f"scenario {i}: DP {out['ttp_dp']} vs periodic {per.max()}"
        assert out["ttp"] >= per.max()
        assert out["schedule"].sum() <= tc.BUDGET
        assert to.rollout(p, out["schedule"][None])[0] == pytest.approx(out["ttp"])
