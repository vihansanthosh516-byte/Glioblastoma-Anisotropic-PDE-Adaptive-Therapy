#!/usr/bin/env python3
"""
Script 75: Equal-budget policies that do NOT receive kill rates
===============================================================
Two follow-ups to the Track C negatives, at the script-59/66 drug budget
(Stupp AUC = 35), on the script-67 evaluation sets, scored on day-90 volume
(scripts 59-67) and on RANO time to progression to day 365 (script 68).

C1  Unconditioned PPO loses to Stupp on the script-64 cohort (25% win), and
    script 67's "blind" policies (not told the kill rates) stay near 40-45%.
    Their only response signal is the day-to-day volume change, but the
    day-90 optimum is a LATE block, so nothing informative is observed before
    the action has to be chosen. Here the policy pays for information instead:

      probe   days 1..L off, L+1..2L chemo, 2L+1..3L RT; each day's log volume
              change is recorded (optionally with multiplicative imaging noise)
      infer   k_chemo = g_off - g_chemo, k_rt = g_off - g_rt, combo = sum
              (combo is additive in every script); no kill rate is given
      commit  either the script-67 efficiency rule on the inferred kills (late
              block of the most kill-per-drug action with the remaining budget),
              or the script-67/68 kill-conditioned networks fed the INFERRED kills

    The probe spends part of the budget (L * (0.5 + 0.75)), all inside the cap.

C2  The script-59 volume-threshold heuristic front-loads combo and runs out of
    budget by day 35. "Paced" caps cumulative AUC at budget * day/90 and
    downgrades the heuristic's action to the largest one under the pace line.

Primary comparison (fixed before running): cohort64 (the script-64 negative),
probe L=1 noise-free + efficiency rule vs Stupp, day-90 win rate.
Everything else is reported in full, including noise levels where it fails.

Sim: s67.KillReactionSim (exact voxel-class compression of the 64^3 model
with diffusion off; script 68 checked TTP against the full PDE on cohort64).
Output: output/probe_paced_policies/evaluate.json
"""
from __future__ import annotations

import json
import math
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Callable, Dict

import numpy as np
import torch
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUTPUT_DIR = PROJECT_ROOT / "output" / "probe_paced_policies"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

_spec = spec_from_file_location("s68", PROJECT_ROOT / "src" / "68_time_to_progression_equal_budget.py")
s68 = module_from_spec(_spec)
_spec.loader.exec_module(s68)
s67, s66 = s68.s67, s68.s66

N_DAYS, BUDGET, STUPP_SCHEDULE = s66.N_DAYS, s66.BUDGET, s66.STUPP_SCHEDULE
INTENSITY = s66.INTENSITY
PROBE_LENGTHS = (1, 3)
NOISE_SIGMAS = (0.0, 0.02, 0.05, 0.10)
NOISE_SEEDS = (0, 1, 2)


# --------------------------------------------------------------------------- #
# C1: probe -> infer -> commit
# --------------------------------------------------------------------------- #
def late_block_sched(best_a: torch.Tensor, remaining: torch.Tensor, start_day: int) -> torch.Tensor:
    """(B, 90) schedule: last n days = best action with n = floor(remaining/I),
    leftover spent on one day of the largest action that fits, placed just before."""
    B = len(best_a)
    S = torch.zeros(B, N_DAYS, dtype=torch.long)
    for b in range(B):
        a, R = int(best_a[b]), float(remaining[b])
        I = float(INTENSITY[a])
        n = min(int(math.floor(R / I + 1e-9)), N_DAYS - start_day)
        if n > 0:
            S[b, N_DAYS - n:] = a
        left = R - n * I
        fill = max([c for c in (1, 2, 3) if float(INTENSITY[c]) <= left + 1e-9], default=0,
                   key=lambda c: float(INTENSITY[c]))
        if fill and N_DAYS - n - 1 >= start_day:
            S[b, N_DAYS - n - 1] = fill
    return S


class ProbePolicy:
    """Stateful batched policy. commit in {"rule", "net"}; for "net", net is a
    kill-conditioned ActorCritic (script 67/68) fed the inferred kill table."""

    def __init__(self, L: int, sigma: float, seed: int, commit: str, net=None):
        self.L, self.sigma, self.commit, self.net = L, sigma, commit, net
        self.g = torch.Generator().manual_seed(seed)
        self.probe = [0] * L + [1] * L + [2] * L

    def _measure(self, v: torch.Tensor) -> torch.Tensor:
        if self.sigma == 0:
            return v
        return v * torch.exp(self.sigma * torch.randn(v.shape, generator=self.g, dtype=v.dtype))

    def __call__(self, obs: Dict) -> torch.Tensor:
        d, B = obs["day"], obs["vol"].shape[0]
        if d == 0:
            self.m_prev, self.r = self._measure(obs["vol"]), []
        else:
            m = self._measure(obs["vol"])
            self.r.append(torch.log(m.clamp_min(1e-12) / self.m_prev.clamp_min(1e-12)))
            self.m_prev = m
        if d < 3 * self.L:
            return torch.full((B,), self.probe[d], dtype=torch.long)
        if d == 3 * self.L:
            R = torch.stack(self.r, 1)
            L = self.L
            g_off, g_c, g_r = R[:, :L].mean(1), R[:, L:2 * L].mean(1), R[:, 2 * L:3 * L].mean(1)
            kc, kr = (g_off - g_c).clamp_min(1e-4), (g_off - g_r).clamp_min(1e-4)
            self.kill_hat = s67.kill_table(kc.numpy(), kr.numpy())
            if self.commit == "rule":
                best_a = s67.efficiency(self.kill_hat).argmax(1) + 1
                self.S = late_block_sched(best_a, obs["budget"] - obs["auc"], d)
            else:
                self.pol = s67.greedy(self.net, self.kill_hat, True)
        if self.commit == "rule":
            return self.S[:, d]
        return self.pol(obs)


# --------------------------------------------------------------------------- #
# C2: paced heuristic
# --------------------------------------------------------------------------- #
def paced_heuristic59(obs: Dict) -> torch.Tensor:
    a = s66.heuristic59_policy(obs)
    allowed = obs["budget"] * (obs["day"] + 1) / N_DAYS - obs["auc"]
    out = torch.zeros_like(a)
    for c in (1, 2, 3):  # ascending intensity; keep the largest that is <= a and fits the pace
        ok = (INTENSITY[c] <= INTENSITY[a]) & (INTENSITY[c] * s66.s59.DT_RL_DAYS <= allowed + 1e-9)
        out = torch.where(ok, torch.full_like(a, c), out)
    return out


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #
def stats_vs(r: Dict, ref: Dict, rng: np.random.Generator) -> Dict:
    f, fs = r["final"].numpy(), ref["final"].numpy()
    lr = np.log(f / fs)
    t, ts = r["ttp"].numpy(), ref["ttp"].numpy()
    dt = t - ts
    boot_lr = np.array([lr[rng.integers(0, len(lr), len(lr))].mean() for _ in range(2000)])
    boot_dt = np.array([dt[rng.integers(0, len(dt), len(dt))].mean() for _ in range(2000)])
    return {
        "day90_win_rate_vs_stupp_pct": float((f < fs).mean() * 100),
        "day90_mean_log_ratio_vs_stupp": float(lr.mean()),
        "day90_mean_log_ratio_ci95": [float(np.percentile(boot_lr, 2.5)), float(np.percentile(boot_lr, 97.5))],
        "day90_wilcoxon_p": float(stats.wilcoxon(np.log(f), np.log(fs)).pvalue) if np.any(lr != 0) else None,
        "ttp_rmst_days": float(t.mean()),
        "ttp_minus_stupp_days_mean": float(dt.mean()),
        "ttp_minus_stupp_days_ci95": [float(np.percentile(boot_dt, 2.5)), float(np.percentile(boot_dt, 97.5))],
        "ttp_n_longer_equal_shorter": [int((dt > 0).sum()), int((dt == 0).sum()), int((dt < 0).sum())],
        "ttp_wilcoxon_p": float(stats.wilcoxon(t, ts).pvalue) if np.any(dt != 0) else None,
        "mean_drug_auc": float(r["auc"].mean()), "max_drug_auc": float(r["auc"].max()),
    }


def arms(kill: torch.Tensor) -> Dict[str, Callable[[], Callable]]:
    """Factories (stateful policies must be rebuilt per rollout)."""
    dag_c, ppo_c = s67.load_run("dagger_cond"), s67.load_run("ppo_cond_ft")
    ttp_net = s67.make_net(True)
    ttp_net.load_state_dict(torch.load(s68.OUTPUT_DIR / "dagger_ttp.pt"))
    out = {
        "stupp": lambda: s66.schedule_policy(STUPP_SCHEDULE),
        # reference arms: blind (no kill rates), and privileged (true kill rates)
        "ppo66_blind": lambda: s66.greedy_policy(s66.load_net("ppo_final_s0")),
        "dagger_blind67": lambda: s67.greedy(s67.load_run("dagger_blind"), kill, False),
        "ppo_blind_ft67": lambda: s67.greedy(s67.load_run("ppo_blind_ft"), kill, False),
        "TRUEKILL_efficiency_rule": lambda: s66.sched_tensor_policy(s67.efficiency_rule_sched(kill)),
        "TRUEKILL_ppo_cond_ft": lambda: s67.greedy(ppo_c, kill, True),
        "TRUEKILL_dagger_ttp": lambda: s67.greedy(ttp_net, kill, True),
        # C2
        "heuristic59": lambda: s66.heuristic59_policy,
        "paced_heuristic59": lambda: paced_heuristic59,
    }
    for L in PROBE_LENGTHS:
        for sig in NOISE_SIGMAS:
            for seed in (NOISE_SEEDS if sig > 0 else (0,)):
                tag = f"L{L}_sig{sig:.2f}_s{seed}"
                out[f"probe_rule_{tag}"] = (lambda L=L, sig=sig, seed=seed: ProbePolicy(L, sig, seed, "rule"))
                out[f"probe_dagger_cond_{tag}"] = (lambda L=L, sig=sig, seed=seed:
                                                   ProbePolicy(L, sig, seed, "net", dag_c))
                out[f"probe_ppo_cond_ft_{tag}"] = (lambda L=L, sig=sig, seed=seed:
                                                   ProbePolicy(L, sig, seed, "net", ppo_c))
                out[f"probe_dagger_ttp_{tag}"] = (lambda L=L, sig=sig, seed=seed:
                                                  ProbePolicy(L, sig, seed, "net", ttp_net))
    return out


def collapse_seeds(block: Dict) -> Dict:
    """Mean over noise seeds of each probe arm's headline numbers."""
    groups: Dict[str, list] = {}
    for name, s in block.items():
        if name.startswith("probe_"):
            groups.setdefault(name.rsplit("_s", 1)[0], []).append(s)
    keys = ["day90_win_rate_vs_stupp_pct", "day90_mean_log_ratio_vs_stupp", "ttp_minus_stupp_days_mean",
            "mean_drug_auc"]
    return {g: {**{k: float(np.mean([r[k] for r in rows])) for k in keys},
                **{k + "_seed_range": [float(min(r[k] for r in rows)), float(max(r[k] for r in rows))]
                   for k in keys[:3]}, "n_seeds": len(rows)}
            for g, rows in groups.items()}


def main():
    rng = np.random.default_rng(75)
    out = {"budget": BUDGET, "horizon_day": s68.HORIZON, "progression_factor": s68.PROG_FACTOR,
           "probe_lengths": PROBE_LENGTHS, "noise_sigmas_log_volume": NOISE_SIGMAS,
           "noise_seeds": NOISE_SEEDS, "sim": "s67.KillReactionSim (reaction-only, exact voxel classes)",
           "primary": "cohort64, probe_rule_L1_sig0.00_s0 vs stupp, day-90 win rate", "sets": {}}
    for set_name, st in s68.reaction_sets().items():
        t0 = time.time()
        rho, kill = st["rho"], st["kill"]
        res = {}
        for name, make in arms(kill).items():
            res[name] = s68.rollout_ttp(s67.KillReactionSim(rho.tolist(), kill), make())
        assert all(float(r["auc"].max()) <= BUDGET + 1e-9 for r in res.values())
        block = {n: stats_vs(r, res["stupp"], rng) for n, r in res.items() if n != "stupp"}
        out["sets"][set_name] = {"source": st["source"], "n": len(rho), "arms": block,
                                 "probe_seed_means": collapse_seeds(block),
                                 "kill_true": kill.tolist(),
                                 "paced_heuristic_actions": res["paced_heuristic59"]["actions"].tolist()}
        print(f"\n[{set_name}] n={len(rho)} ({time.time() - t0:.0f}s)")
        for n, s in block.items():
            if n.startswith("probe_") and not n.endswith("_s0"):
                continue
            print(f"  {n:34s} d90 win {s['day90_win_rate_vs_stupp_pct']:5.1f}%  lr {s['day90_mean_log_ratio_vs_stupp']:+.3f}"
                  f"  TTP {s['ttp_minus_stupp_days_mean']:+6.1f} d  AUC {s['mean_drug_auc']:.1f}")
    path = OUTPUT_DIR / "evaluate.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"[saved] {path}")


if __name__ == "__main__":
    main()
