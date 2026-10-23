"""Track C two-population simulator (script 76): equation, limiting cases, budget rule. Synthetic parameters only;
these check the assumed model behaves as written, not that it describes any patient or drug."""
from __future__ import annotations

import math
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest
import torch

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
_spec = spec_from_file_location("s76_t", REPO / "src" / "76_resistance_adaptive_equal_budget.py")
s76 = module_from_spec(_spec)
_spec.loader.exec_module(s76)

RHO = [0.02, 0.05, 0.1]
KILL = torch.tensor([[0.0, 0.03, 0.05, 0.08]] * 3, dtype=torch.float64)
NO_KILL = torch.zeros(3, 4, dtype=torch.float64)


def run(sim, actions, days=60):
    sim.reset()
    return torch.stack([sim.step(torch.full((sim.B,), actions, dtype=torch.long)) for _ in range(days)], 1)


def test_zero_resistance_matches_single_population():
    for occ, a in [("seed59", 0), ("seed59", 3)]:
        two = s76.TwoPopSim(RHO, KILL, 0.0, 0.25, 0.0, occ)
        one = s76.s67.KillReactionSim(RHO, KILL)
        assert torch.allclose(run(two, a), run(one, a), rtol=1e-12, atol=0)


@pytest.mark.parametrize("f_r0", [0.01, 0.1, 0.5])
def test_no_cost_no_drug_total_follows_single_population(f_r0):
    """With c = 0 and no kill, S and R are the same logistic population, so the total ignores the split."""
    two = s76.TwoPopSim(RHO, NO_KILL, f_r0, 0.0, 0.0, "seed59")
    one = s76.s67.KillReactionSim(RHO, NO_KILL)
    assert torch.allclose(run(two, 0), run(one, 0), rtol=1e-10)


def test_eps_one_no_cost_resistance_is_irrelevant():
    two = s76.TwoPopSim(RHO, KILL, 0.3, 0.0, 1.0, "seed59")
    one = s76.s67.KillReactionSim(RHO, KILL)
    assert torch.allclose(run(two, 3), run(one, 3), rtol=1e-10)


def test_fully_resistant_tumour_ignores_drug():
    on = s76.TwoPopSim(RHO, KILL, 1.0, 0.25, 0.0, "seed59")
    off = s76.TwoPopSim(RHO, KILL, 1.0, 0.25, 0.0, "seed59")
    assert torch.allclose(run(on, 3), run(off, 0), rtol=0, atol=0)


def test_drug_selects_for_resistance_and_cost_slows_regrowth():
    sim = s76.TwoPopSim(RHO, KILL, 0.01, 0.25, 0.0, "seed59")
    run(sim, 3, days=30)
    assert (sim.resistant_fraction() > 0.01).all()
    sim0 = s76.TwoPopSim(RHO, NO_KILL, 0.01, 0.25, 0.0, "seed59")
    run(sim0, 0, days=30)
    assert (sim0.resistant_fraction() < 0.01).all()   # with a cost and no drug, R loses ground


def test_densities_stay_in_bounds():
    sim = s76.TwoPopSim([0.5, 1.0, 2.0], KILL, 0.1, 0.0, 0.2, "near_capacity")
    for a in (0, 3, 0, 3):
        run(sim, a, days=20)
        assert (sim.us >= 0).all() and (sim.ur >= 0).all() and (sim.us + sim.ur <= s76.K * (1 + 1e-12)).all()


def _logistic_rel_error(monkeypatch, n_sub, dt_sub, rho=0.05, days=100):
    monkeypatch.setattr(s76, "N_SUB", n_sub)
    monkeypatch.setattr(s76, "DT_SUB", dt_sub)
    sim = s76.TwoPopSim([rho], NO_KILL[:1], 0.0, 0.0, 0.0, "seed59")
    u0 = sim.u0[0].clone()
    run(sim, 0, days=days)
    exact = s76.K / (1 + (s76.K / u0 - 1) * math.exp(-rho * days))
    return ((sim.u[0] - exact).abs() / exact).max().item()


def test_untreated_growth_converges_to_logistic_solution(monkeypatch):
    """Explicit Euler vs the closed-form logistic curve. One step() is one day. The default step (0.2 d) has a
    first-order error of about rho^2 * dt * t / 2 (2.5% here); halving dt must roughly halve it."""
    n, dt = s76.N_SUB, s76.DT_SUB
    assert math.isclose(n * dt, 1.0)
    e1 = _logistic_rel_error(monkeypatch, n, dt)
    e2 = _logistic_rel_error(monkeypatch, 2 * n, dt / 2)
    assert e1 < 0.05 ** 2 * dt * 100 / 2 * 1.1
    assert 1.8 < e1 / e2 < 2.2


def test_rollout_respects_budget_and_window():
    sim = s76.TwoPopSim(RHO, KILL, 0.01, 0.25, 0.0, "seed59")
    always = lambda obs: torch.full((obs["vol"].shape[0],), 3, dtype=torch.long)  # noqa: E731
    r = s76.rollout(sim, always, window=365)
    assert torch.allclose(r["auc"], torch.full((3,), s76.BUDGET))
    acts = r["actions"]
    assert (acts[:, int(s76.BUDGET):] == 0).all()          # budget spent on day 35, holidays after
    r = s76.rollout(sim, always, window=10)
    assert (r["actions"][:, 10:] == 0).all() and (r["auc"] <= s76.BUDGET).all()


@pytest.mark.parametrize("window", [90, 365])
def test_paced_policy_spends_budget_evenly_in_window(window):
    sim = s76.TwoPopSim(RHO, KILL, 0.01, 0.25, 0.0, "seed59")
    r = s76.rollout(sim, s76.paced_policy(window), window)
    on = (r["actions"][0] == 3).nonzero().flatten()
    assert len(on) == int(s76.BUDGET) and on.max() < window
    assert torch.allclose(r["auc"], torch.full((3,), s76.BUDGET))


def test_adaptive_rule_switches_on_thresholds():
    pol = s76.at_policy(0.5)
    v0 = torch.tensor([10.0, 10.0])
    seq = [(0, [10.0, 10.0], [3, 3]), (1, [6.0, 4.0], [3, 0]), (2, [5.0, 9.0], [0, 0]), (3, [7.0, 10.0], [0, 3])]
    for day, vol, want in seq:
        got = pol({"day": day, "vol": torch.tensor(vol), "v0": v0})
        assert got.tolist() == want


def test_rollout_is_deterministic():
    a = s76.rollout(s76.TwoPopSim(RHO, KILL, 0.01, 0.25, 0.0, "seed59"), s76.at_policy(0.5), 365)
    b = s76.rollout(s76.TwoPopSim(RHO, KILL, 0.01, 0.25, 0.0, "seed59"), s76.at_policy(0.5), 365)
    assert torch.equal(a["actions"], b["actions"]) and torch.equal(a["ttp"], b["ttp"])
