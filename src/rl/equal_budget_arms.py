"""
Equal-drug-budget treatment arms shared by scripts 60, 62, 64 and 65.

Four arms, all defined once here (policies and budget rule come from script 66):
  stupp               fixed Stupp schedule (the budget reference)
  heuristic_budgeted  script 59 volume-threshold heuristic, capped at the budget
  ppo                 PPO policy from script 66 (learned: combo from day 1)
  dagger_oracle       DAgger policy from script 66 (oracle: combo on the last 35 days)

Each arm is called as arm(solver, budget) and returns per-patient lists of
final_volume_mm3, drug_auc and action_history (plus mean burden and daily
volumes). `solver` is a batched solver: BatchedSolver wraps any script's own
numpy FastPDESolver instances (one per patient) so each script keeps its own
physics, kill rates and ablation flags.
"""
from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Callable, Dict, List, Sequence

import numpy as np
import torch

SRC = Path(__file__).resolve().parent.parent
_spec = spec_from_file_location("s66", SRC / "66_rl_equal_budget_study.py")
s66 = module_from_spec(_spec)
_spec.loader.exec_module(s66)

STUPP_BUDGET: float = s66.BUDGET
STUPP_SCHEDULE: List[int] = s66.STUPP_SCHEDULE
N_DAYS: int = s66.N_DAYS
WEIGHTS = {"ppo": "ppo_final_s0", "dagger_oracle": "bc_dagger_final"}
FACE_FIELDS = ["Dxx_xf", "Dxy_xf", "Dxz_xf", "Dyy_yf", "Dxy_yf", "Dyz_yf", "Dzz_zf", "Dxz_zf", "Dyz_zf"]


class BatchedSolver(s66.BatchedPDE):
    """Batched torch copy of a list of numpy FastPDESolver objects.

    kill_table[b][a] is the kill rate of action a for patient b, exactly as that
    script's rl_step computes it. The seed is each solver's reset() state times
    0.1, the initial condition every protocol in scripts 59-66 uses."""

    def __init__(self, solvers: Sequence, kill_table):
        s0 = solvers[0]
        for s in solvers:
            assert abs(s.dt * s66.N_SUB - s66.s59.DT_RL_DAYS) < 1e-12, "expected 5 sub-steps per day"
            assert (s.nx, s.ny, s.nz, s.dx) == (s0.nx, s0.ny, s0.nz, s0.dx)
        self.n, self.dx, self.B, self.dt = s0.nx, s0.dx, len(solvers), s0.dt
        self.rho = torch.tensor([float(s.rho) for s in solvers]).view(-1, 1, 1, 1)
        self.face = {f: torch.from_numpy(np.stack([getattr(s, f) for s in solvers])) for f in FACE_FIELDS}
        seeds = []
        for s in solvers:
            s.reset()
            seeds.append(s.u * 0.1)
        self.u0 = torch.from_numpy(np.stack(seeds))
        self.kill = torch.as_tensor(np.asarray(kill_table, dtype=float)).view(self.B, 4)
        self.reset()

    def step(self, actions: torch.Tensor) -> torch.Tensor:
        kill = self.kill.gather(1, actions.view(-1, 1)).view(-1, 1, 1, 1)
        K = s66.s59.K_CARRY
        for _ in range(s66.N_SUB):
            u = self.u
            self.u = torch.clamp(u + self.dt * (self._div(u) + self.rho * u * (1.0 - u / K) - kill * u), 0.0, K)
        return self.volume()


def validate_against_numpy(solver, kill_row: Sequence[float]) -> float:
    """Max relative daily-volume error of BatchedSolver vs the script's own numpy
    solver on the Stupp schedule. Catches a kill_table that disagrees with rl_step."""
    batched = s66.run_schedules(BatchedSolver([solver], [kill_row]), torch.tensor([STUPP_SCHEDULE]))[0]
    solver.reset()
    solver.u *= 0.1
    ref = []
    for a in STUPP_SCHEDULE:
        solver.rl_step(a)
        ref.append(float(solver.u.sum() * solver.dx ** 3))
    ref = torch.tensor(ref)
    return float(((batched - ref).abs() / ref).max())


# --------------------------------------------------------------------------- #
# Arms
# --------------------------------------------------------------------------- #
_nets: Dict[str, "s66.ActorCritic"] = {}


def _policy(arm: str) -> Callable:
    if arm not in _nets:
        _nets[arm] = s66.load_net(WEIGHTS[arm])
    return s66.greedy_policy(_nets[arm])


def _run(solver, policy: Callable, budget: float) -> Dict[str, list]:
    res = s66.rollout(solver, policy, budget)
    return {"final_volume_mm3": res["final"].tolist(), "drug_auc": res["auc"].tolist(),
            "action_history": res["actions"].tolist(), "mean_burden_mm3": res["mean_burden"].tolist(),
            "volume_history_mm3": res["vols"].tolist(), "initial_volume_mm3": res["v0"].tolist()}


def stupp(solver, budget: float = STUPP_BUDGET) -> Dict[str, list]:
    return _run(solver, s66.schedule_policy(STUPP_SCHEDULE), budget)


def heuristic_budgeted(solver, budget: float = STUPP_BUDGET) -> Dict[str, list]:
    return _run(solver, s66.heuristic59_policy, budget)


def ppo(solver, budget: float = STUPP_BUDGET) -> Dict[str, list]:
    return _run(solver, _policy("ppo"), budget)


def dagger_oracle(solver, budget: float = STUPP_BUDGET) -> Dict[str, list]:
    return _run(solver, _policy("dagger_oracle"), budget)


ARMS: Dict[str, Callable] = {"stupp": stupp, "heuristic_budgeted": heuristic_budgeted,
                             "ppo": ppo, "dagger_oracle": dagger_oracle}


def evaluate_arms(make_solvers: Sequence[Callable], kill_table, budget: float = STUPP_BUDGET,
                  arms: Dict[str, Callable] = ARMS, chunk: int = 30) -> Dict[str, Dict[str, list]]:
    """Run every arm on every patient. make_solvers[i]() builds patient i's numpy
    solver; patients are processed in chunks to bound memory."""
    out = {name: {} for name in arms}
    for i in range(0, len(make_solvers), chunk):
        solver = BatchedSolver([f() for f in make_solvers[i:i + chunk]], list(kill_table)[i:i + chunk])
        for name, arm in arms.items():
            for k, v in arm(solver, budget).items():
                out[name].setdefault(k, []).extend(v)
        del solver
    return out
