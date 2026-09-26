#!/usr/bin/env python3
"""
Script 67: Equal-budget policies that know the patient's drug sensitivity
=========================================================================
Follow-up to script 64, where the script-66 policies lost to Stupp at equal
drug (PPO win 25%, mean log-ratio +0.077). Those policies were trained with one
fixed kill-rate table (chemo 0.05, RT 0.08, combo 0.13), so they always chose
combo. In the virtual cohort 15/20 patients kill more tumour per unit of drug
with chemo alone, and Stupp's chemo-only days beat an all-combo policy.

Here training randomises each patient's kill rates (domain randomisation,
Tobin et al. 2017) and the policy is either
  cond   told the patient's kill rate per unit drug for each action, or
  blind  not told, so it can only infer sensitivity from the volume response.
Clinically, "cond" corresponds to a pre-treatment sensitivity biomarker
(e.g. MGMT promoter methylation for TMZ, Hegi et al. 2005); "blind" to pure
response-adaptive dosing.

Physics, seed tumour, intensity scale, budget rule and budget (= Stupp AUC)
are those of scripts 59/66, unchanged. Every arm is evaluated at that budget.

Stages (run all with no flags, or one with --stage):
  oracle    privileged per-patient optimum (structured search, CEM check);
            how often it equals "late block of the most drug-efficient action"
  train     DAgger (cond, blind) from the oracle; PPO from scratch (cond);
            PPO fine-tuned from DAgger (cond, blind; demonstration-initialised
            RL as in Rajeswaran et al. 2018)
  evaluate  full 64^3 PDE on: script-64 cohort (its own solver and kills),
            script-60 LHS scenarios (its own solver and kills), script-66 real
            test patients, and held-out randomised patients
  report    figure

Outputs: output/rl_kill_conditioned/{oracle,train,evaluate}.json,
         output/rl_kill_conditioned/kill_conditioned_policy.png
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import torch
import torch.nn.functional as F
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from rl.equal_budget_arms import BatchedSolver, s66, validate_against_numpy  # noqa: E402

OUTPUT_DIR = PROJECT_ROOT / "output" / "rl_kill_conditioned"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
s59 = s66.s59


def _load(name: str, file: str):
    spec = spec_from_file_location(name, PROJECT_ROOT / "src" / file)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


N_DAYS, N_SUB, DT_SUB = s66.N_DAYS, s66.N_SUB, s66.DT_SUB
INTENSITY, BUDGET, STUPP_SCHEDULE = s66.INTENSITY, s66.BUDGET, s66.STUPP_SCHEDULE
ACTION_NAMES = ["off", "chemo", "rad", "combo"]

# Training distribution. Covers the fixed 59/66 table (0.05, 0.08), script 60's
# alpha_sens-scaled table (chemo 0.025-0.075, RT 0.04-0.12) and script 64's
# ranges (chemo 0.02-0.08, RT 0.015-0.045). Combo is additive in every script.
KILL_RANGES = {"chemo": (0.01, 0.10), "rad": (0.01, 0.13)}
RHO_LOGUNIFORM = (0.002, 0.12)      # half the draws; the other half resample real training rho
N_TRAIN_POOL, N_VAL_POOL, N_SYNTH_TEST = 2000, 300, 40
SEED_TRAIN, SEED_VAL, SEED_TEST = 670, 671, 672


# --------------------------------------------------------------------------- #
# Simulators with a per-patient kill table
# --------------------------------------------------------------------------- #
class KillReactionSim(s66.ReactionSim):
    """s66.ReactionSim (exact voxel-class compression, no diffusion) with kill[b, a]."""

    def __init__(self, rho: Sequence[float], kill: torch.Tensor):
        super().__init__(list(rho))
        self.kill = torch.as_tensor(kill, dtype=torch.float64).view(self.B, 4)

    def step(self, actions: torch.Tensor) -> torch.Tensor:
        kill = self.kill.gather(1, actions.view(-1, 1))
        K = s59.K_CARRY
        for _ in range(N_SUB):
            u = self.u
            self.u = torch.clamp(u + DT_SUB * (self.rho * u * (1.0 - u / K) - kill * u), 0.0, K)
        return self.volume()


class KillPDE(s66.BatchedPDE):
    """s66.BatchedPDE (validated port of the script-59 PDE, diffusion on) with kill[b, a]."""

    def __init__(self, rho: Sequence[float], D_w: Sequence[float], kill: torch.Tensor, grid: int = 64):
        super().__init__(grid, list(rho), list(D_w))
        self.kill = torch.as_tensor(kill, dtype=torch.float64).view(self.B, 4)

    def step(self, actions: torch.Tensor) -> torch.Tensor:
        kill = self.kill.gather(1, actions.view(-1, 1)).view(-1, 1, 1, 1)
        K = s59.K_CARRY
        for _ in range(N_SUB):
            u = self.u
            self.u = torch.clamp(u + DT_SUB * (self._div(u) + self.rho * u * (1.0 - u / K) - kill * u), 0.0, K)
        return self.volume()


def kill_table(chemo: np.ndarray, rad: np.ndarray) -> torch.Tensor:
    chemo, rad = np.asarray(chemo, float), np.asarray(rad, float)
    return torch.tensor(np.stack([np.zeros_like(chemo), chemo, rad, chemo + rad], 1))


def sample_patients(n: int, seed: int, real_rho: np.ndarray) -> Dict[str, torch.Tensor]:
    rng = np.random.default_rng(seed)
    use_real = rng.random(n) < 0.5
    lo, hi = np.log(RHO_LOGUNIFORM[0]), np.log(RHO_LOGUNIFORM[1])
    rho = np.where(use_real, rng.choice(real_rho, n), np.exp(rng.uniform(lo, hi, n)))
    chemo = rng.uniform(*KILL_RANGES["chemo"], n)
    rad = rng.uniform(*KILL_RANGES["rad"], n)
    return {"rho": torch.tensor(rho), "kill": kill_table(chemo, rad)}


def efficiency(kill: torch.Tensor) -> torch.Tensor:
    """Kill rate per unit of drug AUC for chemo, rad, combo -> (B, 3)."""
    return kill[:, 1:] / INTENSITY[1:]


# --------------------------------------------------------------------------- #
# Privileged oracle and a non-learned rule
# --------------------------------------------------------------------------- #
def mixed_candidates() -> List[Dict]:
    """Two-action late blocks that pack the budget exactly: m days of action b, then n
    days of action a ending on day 90, for every split (n, m); the leftover budget goes
    to one extra day of the largest action that fits. When two actions have nearly the
    same kill per unit drug, a mix spends the budget more completely than one block
    (found by the unstructured CEM check)."""
    out = []
    for a in (1, 2, 3):
        for b in (1, 2, 3):
            if a == b:
                continue
            Ia, Ib = float(INTENSITY[a]), float(INTENSITY[b])
            for n in range(1, int(math.floor(BUDGET / Ia + 1e-9)) + 1):
                m = int(math.floor((BUDGET - n * Ia) / Ib + 1e-9))
                left = BUDGET - n * Ia - m * Ib
                fill = max([c for c in (1, 2, 3) if float(INTENSITY[c]) <= left + 1e-9], default=0,
                           key=lambda c: float(INTENSITY[c]))
                if m == 0 or n + m + (fill > 0) > N_DAYS:
                    continue
                sched = [0] * (N_DAYS - n - m) + [b] * m + [a] * n
                if fill:
                    sched[N_DAYS - n - m - 1] = fill
                out.append({"kind": "mixed", "action": a, "n": n, "second": b, "n_second": m, "fill": fill,
                            "start": N_DAYS - n - m + 1, "sched": sched})
    return out


CANDS = [{"kind": "block", **c} for c in s66.candidate_schedules(INTENSITY, BUDGET)] + mixed_candidates()
CAND_SCHED = torch.tensor([c["sched"] for c in CANDS])
assert all(sum(float(INTENSITY[a]) for a in c["sched"]) <= BUDGET + 1e-9 for c in CANDS)


def kill_oracle(rho: torch.Tensor, kill: torch.Tensor, chunk: int = 6000) -> Dict:
    """Best of the structured candidates for each patient: every single-action block
    that spends the budget at every start day (plus a one-day fill), and every
    two-action late block that packs the budget (mixed_candidates)."""
    P, K = len(rho), len(CANDS)
    ip, ik = torch.arange(P).repeat_interleave(K), torch.arange(K).repeat(P)
    final, burden = torch.empty(P * K), torch.empty(P * K)
    for i in range(0, P * K, chunk):
        sl = slice(i, i + chunk)
        res = s66.rollout(KillReactionSim(rho[ip[sl]].tolist(), kill[ip[sl]]),
                          s66.sched_tensor_policy(CAND_SCHED[ik[sl]]))
        final[sl], burden[sl] = res["final"], res["mean_burden"]
    final, burden = final.view(P, K), burden.view(P, K)
    best = final.argmin(1)
    return {"best": best, "sched": CAND_SCHED[best], "final": final, "burden": burden}


def efficiency_rule_sched(kill: torch.Tensor) -> torch.Tensor:
    """Non-learned rule: the most kill-per-unit-drug action, as a late block
    (latest start that spends the budget, with the candidate fill day)."""
    best_a = efficiency(kill).argmax(1) + 1
    idx = []
    for a in best_a.tolist():
        c = [k for k, c in enumerate(CANDS) if c["kind"] == "block" and c["action"] == a]
        idx.append(max(c, key=lambda k: CANDS[k]["start"]))
    return CAND_SCHED[idx]


def cem_kill(rho: Sequence[float], kill: torch.Tensor, n: int = 256, iters: int = 40, elite: float = 0.1,
             alpha: float = 0.7, seed: int = 0, init: Optional[torch.Tensor] = None) -> Dict:
    """s66.cem_search (unstructured, all 4^90 schedules) with per-patient kills.
    `init` (P, 90) warm-starts the sampling distribution at 50% on that schedule."""
    g = torch.Generator().manual_seed(seed)
    P = len(rho)
    probs = torch.full((P, N_DAYS, 4), 0.25)
    if init is not None:
        probs = 0.5 * F.one_hot(init.long(), 4).double() + 0.5 * probs
    best_val = torch.full((P,), float("inf"))
    best_sched = torch.zeros(P, N_DAYS, dtype=torch.long)
    n_el = max(1, int(n * elite))
    sim = KillReactionSim([r for r in rho for _ in range(n)], kill.repeat_interleave(n, 0))
    for _ in range(iters):
        S = torch.multinomial(probs.repeat_interleave(n, 0).view(-1, 4), 1, generator=g).view(P * n, N_DAYS)
        res = s66.rollout(sim, s66.sched_tensor_policy(S))
        val, exe = res["final"].view(P, n), res["actions"].view(P, n, N_DAYS)
        order = val.argsort(1)[:, :n_el]
        el = torch.gather(exe, 1, order.unsqueeze(2).expand(-1, -1, N_DAYS))
        probs = alpha * F.one_hot(el, 4).double().mean(1) + (1 - alpha) * probs
        cur = val.min(1)
        better = cur.values < best_val
        best_val = torch.where(better, cur.values, best_val)
        best_sched[better] = exe[better, cur.indices[better]]
    return {"best_value": best_val, "best_sched": best_sched}


def pools() -> Dict[str, Dict[str, torch.Tensor]]:
    real_rho = np.array([p["rho"] for p in s66.real_cohort()["train"]])
    return {"train": sample_patients(N_TRAIN_POOL, SEED_TRAIN, real_rho),
            "val": sample_patients(N_VAL_POOL, SEED_VAL, real_rho),
            "synth_test": sample_patients(N_SYNTH_TEST, SEED_TEST, real_rho)}


def save_json(obj, name: str):
    path = OUTPUT_DIR / name
    path.write_text(json.dumps(obj, indent=2))
    print(f"[saved] {path}")


def stage_oracle():
    out = {"budget": BUDGET, "kill_ranges": KILL_RANGES, "rho_loguniform": RHO_LOGUNIFORM,
           "n_candidates": len(CANDS), "sim": "reaction-only 64^3 (exact voxel classes)"}
    oracle_store = {}
    for name, pool in pools().items():
        t0 = time.time()
        o = kill_oracle(pool["rho"], pool["kill"])
        sim = lambda: KillReactionSim(pool["rho"].tolist(), pool["kill"])  # noqa: E731
        stupp = s66.rollout(sim(), s66.schedule_policy(STUPP_SCHEDULE))
        rule = s66.rollout(sim(), s66.sched_tensor_policy(efficiency_rule_sched(pool["kill"])))
        best_desc = [CANDS[b] for b in o["best"].tolist()]
        eff_a = (efficiency(pool["kill"]).argmax(1) + 1).tolist()
        latest = {a: max(c["start"] for c in CANDS if c["kind"] == "block" and c["action"] == a) for a in (1, 2, 3)}
        of = o["final"].gather(1, o["best"].view(-1, 1)).squeeze(1)
        blocks = [d for d in best_desc if d["kind"] == "block"]
        single = torch.tensor([c["kind"] == "block" for c in CANDS])
        of_single = o["final"][:, single].min(1).values
        out[name] = {
            "n": len(pool["rho"]), "seconds": time.time() - t0,
            "oracle_kind_counts": {k: sum(d["kind"] == k for d in best_desc) for k in ("block", "mixed")},
            "oracle_main_action_counts": {ACTION_NAMES[a]: sum(d["action"] == a for d in best_desc)
                                          for a in (1, 2, 3)},
            "oracle_main_action_equals_most_efficient_pct": 100 * float(np.mean([d["action"] == a for d, a
                                                                                in zip(best_desc, eff_a)])),
            "oracle_block_start_is_latest_pct": 100 * float(np.mean([d["start"] == latest[d["action"]]
                                                                    for d in blocks])) if blocks else None,
            "mixed_gain_over_single_block_log": float(torch.log(of / of_single).mean()),
            "mixed_gain_over_single_block_log_min": float(torch.log(of / of_single).min()),
            "oracle_final_log_ratio_vs_stupp": float(torch.log(of / stupp["final"]).mean()),
            "oracle_win_vs_stupp_pct": float((of < stupp["final"]).double().mean() * 100),
            "rule_final_log_ratio_vs_stupp": float(torch.log(rule["final"] / stupp["final"]).mean()),
            "rule_final_log_ratio_vs_oracle": float(torch.log(rule["final"] / of).mean()),
            "rule_final_log_ratio_vs_oracle_max": float(torch.log(rule["final"] / of).max()),
        }
        oracle_store[name] = {"rho": pool["rho"], "kill": pool["kill"], "sched": o["sched"], "final": of,
                              "stupp_final": stupp["final"]}
        print(f"  [{name}] {out[name]}")
    torch.save(oracle_store, OUTPUT_DIR / "oracle_pools.pt")

    # CEM check: does unstructured search beat the structured optimum? Patients chosen
    # to span the most-efficient action and rho.
    v = oracle_store["val"]
    eff_a = efficiency(v["kill"]).argmax(1) + 1
    pick = []
    for a in (1, 2, 3):
        ids = torch.nonzero(eff_a == a).flatten()
        ids = ids[v["rho"][ids].argsort()]
        pick += ids[torch.linspace(0, len(ids) - 1, 3).long()].tolist()
    out["cem_check"] = {"val_indices": pick, "rho": v["rho"][pick].tolist(),
                        "most_efficient_action": [ACTION_NAMES[a] for a in eff_a[pick].tolist()]}
    for mode, init in (("cold", None), ("warm_from_oracle", v["sched"][pick])):
        t0 = time.time()
        c = cem_kill(v["rho"][pick].tolist(), v["kill"][pick], init=init)
        rel = ((c["best_value"] - v["final"][pick]) / v["final"][pick]).tolist()
        out["cem_check"][mode] = {"cem_minus_structured_rel": rel, "cem_best_sched": c["best_sched"].tolist(),
                                  "seconds": time.time() - t0}
        print(f"  CEM {mode} ({time.time() - t0:.0f}s): rel diff {['%+.4f' % x for x in rel]}")
    save_json(out, "oracle.json")


# --------------------------------------------------------------------------- #
# Policies
# --------------------------------------------------------------------------- #
def feats(obs: Dict, kill: torch.Tensor, cond: bool) -> torch.Tensor:
    f = s66.features(obs)
    return torch.cat([f, efficiency(kill) * 10], 1) if cond else f


def make_net(cond: bool) -> s66.ActorCritic:
    return s66.ActorCritic(obs_dim=8 if cond else 5)


def greedy(net: s66.ActorCritic, kill: torch.Tensor, cond: bool) -> Callable:
    def pol(obs):
        with torch.no_grad():
            return net.pi(feats(obs, kill, cond)).masked_fill(~s66.action_mask(obs), -1e9).argmax(1)
    return pol


def eval_pool(pol_factory: Callable, pool: Dict) -> Dict:
    """Reaction-sim evaluation on a pool with stored oracle and Stupp finals."""
    res = s66.rollout(KillReactionSim(pool["rho"].tolist(), pool["kill"]), pol_factory(pool["kill"]))
    ls = torch.log(res["final"] / pool["stupp_final"])
    lo = torch.log(res["final"] / pool["final"])
    return {"final_log_ratio_vs_stupp": float(ls.mean()), "win_vs_stupp_pct": float((ls < 0).double().mean() * 100),
            "final_log_ratio_vs_oracle": float(lo.mean()), "mean_auc": float(res["auc"].mean())}


def val_score(row: Dict) -> float:
    """Checkpoint-selection score on the validation pool (lower is better). Script 68
    validation rows carry 'rmst' (days, higher is better) instead."""
    return -row["rmst"] if "rmst" in row else row["final_log_ratio_vs_stupp"]


def train_dagger(cond: bool, train: Dict, val: Dict, seed: int = 0, rounds: int = 6, steps: int = 3000,
                 mb: int = 2048, val_eval: Optional[Callable] = None) -> (s66.ActorCritic, List[Dict]):
    """DAgger with the privileged oracle as teacher (Ross et al. 2011). The teacher's
    label in any state is the patient's oracle-schedule action for that day, relabelled
    to a holiday if the budget no longer allows it. A blind student sees only volumes
    ("learning by cheating", Chen et al. 2019). val_eval(net) overrides the default
    day-90 validation metrics (script 68 passes a time-to-progression one).
    Returns the round with the best validation score (early stopping on the
    validation pool; the evaluation sets are never used for selection)."""
    torch.manual_seed(seed)
    net = make_net(cond)
    best_state, best = None, float("inf")
    opt = torch.optim.Adam(net.pi.parameters(), lr=1e-3)
    F_all, M_all, Y_all, log = [], [], [], []
    S = train["sched"]
    for r in range(rounds):
        student = greedy(net, train["kill"], cond)

        def pol(obs, r=r):
            m = s66.action_mask(obs)
            y = S[:, obs["day"]]
            y = torch.where(m.gather(1, y.view(-1, 1)).squeeze(1), y, torch.zeros_like(y))
            F_all.append(feats(obs, train["kill"], cond)); M_all.append(m); Y_all.append(y)
            return y if r == 0 else student(obs)

        s66.rollout(KillReactionSim(train["rho"].tolist(), train["kill"]), pol)
        f, m, y = torch.cat(F_all), torch.cat(M_all), torch.cat(Y_all)
        for _ in range(steps):
            i = torch.randint(len(y), (mb,))
            loss = F.cross_entropy(net.pi(f[i]).masked_fill(~m[i], -1e9), y[i])
            opt.zero_grad(); loss.backward(); opt.step()
        with torch.no_grad():
            acc = float((net.pi(f).masked_fill(~m, -1e9).argmax(1) == y).double().mean())
        row = {"round": r, "n_labels": len(y), "train_acc": acc,
               "val": val_eval(net) if val_eval else eval_pool(lambda k: greedy(net, k, cond), val)}
        log.append(row)
        print(f"    DAgger {'cond' if cond else 'blind'} round {r}: acc {acc:.4f}  val {row['val']}")
        if val_score(row["val"]) < best:
            best, best_state = val_score(row["val"]), copy.deepcopy(net.state_dict())
            row["selected"] = True
    net.load_state_dict(best_state)
    return net, log


def train_ppo_k(cond: bool, train: Dict, val: Dict, seed: int, iters: int, batch: int = 128,
                init: Optional[s66.ActorCritic] = None, lr: float = 3e-4, ent_coef: float = 0.01,
                critic_warmup: int = 0, eval_every: int = 10) -> (s66.ActorCritic, Dict):
    """s66.train_ppo (terminal reward -log(V_90 / V_90,Stupp), GAE 0.95, clip 0.2) with
    patients drawn from the randomised pool. With `init`, the actor starts from a DAgger
    policy and the critic is fitted alone for `critic_warmup` iterations first. Returns
    the evaluated checkpoint with the best validation score (early stopping)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    net = make_net(cond)
    if init is not None:
        net.pi.load_state_dict(init.pi.state_dict())
    opt_pi = torch.optim.Adam(net.pi.parameters(), lr=lr)
    opt_v = torch.optim.Adam(net.v.parameters(), lr=1e-3)
    curve, stats_ = [], []
    best_state, best, best_iter = None, float("inf"), None
    for it in range(iters + 1):
        if it % eval_every == 0 or it == iters:
            curve.append({"iter": it, "val": eval_pool(lambda k: greedy(net, k, cond), val)})
            if val_score(curve[-1]["val"]) < best:
                best, best_state, best_iter = val_score(curve[-1]["val"]), copy.deepcopy(net.state_dict()), it
        if it == iters:
            break
        idx = rng.integers(len(train["rho"]), size=batch)
        kill = train["kill"][idx]
        sim = KillReactionSim(train["rho"][idx].tolist(), kill)
        ref = s66.rollout(sim, s66.schedule_policy(STUPP_SCHEDULE))
        buf = {"f": [], "m": [], "a": [], "logp": [], "v": []}

        def pol(obs):
            f, m = feats(obs, kill, cond), s66.action_mask(obs)
            with torch.no_grad():
                dist = torch.distributions.Categorical(logits=net.pi(f).masked_fill(~m, -1e9))
                a = dist.sample()
                buf["f"].append(f); buf["m"].append(m); buf["a"].append(a)
                buf["logp"].append(dist.log_prob(a)); buf["v"].append(net.v(f).squeeze(1))
            return a

        res = s66.rollout(sim, pol)
        assert torch.equal(res["actions"], torch.stack(buf["a"], 1)), "mask/budget rule mismatch"
        R = s66.compute_rewards({"reward": "final_rel"}, res, ref)
        Vt = torch.stack(buf["v"], 1)
        adv = torch.zeros_like(R)
        last = torch.zeros(R.shape[0])
        for t in reversed(range(N_DAYS)):
            nxt = Vt[:, t + 1] if t + 1 < N_DAYS else torch.zeros_like(last)
            last = R[:, t] + nxt - Vt[:, t] + 0.95 * last
            adv[:, t] = last
        ret = (adv + Vt).flatten()
        f = torch.stack(buf["f"], 1).flatten(0, 1)
        m = torch.stack(buf["m"], 1).flatten(0, 1)
        a = torch.stack(buf["a"], 1).flatten()
        old = torch.stack(buf["logp"], 1).flatten()
        advf = adv.flatten()
        advf = (advf - advf.mean()) / (advf.std() + 1e-8)
        kls, ents = [], []
        for _ in range(4):
            for b in torch.randperm(len(a)).split(1024):
                if it >= critic_warmup:
                    dist = torch.distributions.Categorical(logits=net.pi(f[b]).masked_fill(~m[b], -1e9))
                    lp = dist.log_prob(a[b])
                    ratio = torch.exp(lp - old[b])
                    pl = -torch.min(ratio * advf[b], ratio.clamp(0.8, 1.2) * advf[b]).mean()
                    ent = dist.entropy().mean()
                    opt_pi.zero_grad()
                    (pl - ent_coef * ent).backward()
                    torch.nn.utils.clip_grad_norm_(net.pi.parameters(), 0.5)
                    opt_pi.step()
                    kls.append(float((old[b] - lp).mean().detach()))
                    ents.append(float(ent.detach()))
                vl = F.mse_loss(net.v(f[b]).squeeze(1), ret[b])
                opt_v.zero_grad()
                vl.backward()
                torch.nn.utils.clip_grad_norm_(net.v.parameters(), 0.5)
                opt_v.step()
        stats_.append({"iter": it, "train_final_log_ratio": float(torch.log(res["final"] / ref["final"]).mean()),
                       "entropy": float(np.mean(ents)) if ents else None,
                       "approx_kl": float(np.mean(kls)) if kls else None, "mean_auc": float(res["auc"].mean())})
        if it % eval_every == 0:
            print(f"    iter {it:4d}  val {curve[-1]['val']}")
    net.load_state_dict(best_state)
    return net, {"curve": curve, "stats": stats_, "selected_iter": best_iter}


RUNS = ["dagger_cond", "dagger_blind", "ppo_cond", "ppo_cond_ft", "ppo_blind_ft"]


def load_run(name: str) -> s66.ActorCritic:
    net = make_net("blind" not in name)
    net.load_state_dict(torch.load(OUTPUT_DIR / f"{name}.pt"))
    return net


def stage_train(runs: Optional[List[str]] = None, ppo_iters: int = 300, ft_iters: int = 200):
    store = torch.load(OUTPUT_DIR / "oracle_pools.pt")
    train, val = store["train"], store["val"]
    out_path = OUTPUT_DIR / "train.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    out.setdefault("runs", {})
    ref_rows = {}
    for name, factory in (("stupp", lambda k: s66.schedule_policy(STUPP_SCHEDULE)),
                          ("efficiency_rule", lambda k: s66.sched_tensor_policy(efficiency_rule_sched(k))),
                          ("ppo66", lambda k: s66.greedy_policy(s66.load_net("ppo_final_s0"))),
                          ("dagger66", lambda k: s66.greedy_policy(s66.load_net("bc_dagger_final")))):
        ref_rows[name] = eval_pool(factory, val)
    out["val_reference"] = ref_rows
    print(f"  val reference: {ref_rows}")

    def done(name, net, log):
        torch.save(net.state_dict(), OUTPUT_DIR / f"{name}.pt")
        out["runs"][name] = log
        out_path.write_text(json.dumps(out, indent=2))

    for name in RUNS:
        if runs and name not in runs:
            continue
        t0 = time.time()
        cond = "blind" not in name
        print(f"  [{name}]")
        if name.startswith("dagger"):
            net, log = train_dagger(cond, train, val)
            log = {"dagger": log}
        elif name == "ppo_cond":
            net, log = train_ppo_k(cond, train, val, seed=0, iters=ppo_iters)
        else:
            init = load_run(name.replace("ppo_", "dagger_").replace("_ft", ""))
            net, log = train_ppo_k(cond, train, val, seed=0, iters=ft_iters, init=init, lr=1e-4,
                                   ent_coef=0.001, critic_warmup=20)
        log["seconds"] = time.time() - t0
        done(name, net, log)
        print(f"  [{name}] {log['seconds']:.0f}s")
    print(f"[saved] {out_path}")


# --------------------------------------------------------------------------- #
# Evaluation on the full PDE
# --------------------------------------------------------------------------- #
def eval_sets() -> Dict[str, Dict]:
    """Each set: solver factory (batched full PDE), kill table, rho, and ids."""
    s64 = _load("s64", "64_virtual_cohort_simulation.py")
    s60 = _load("s60", "60_baselines_and_ablation.py")
    sets = {}

    cohort = s64.generate_virtual_cohort(s64.N_PATIENTS, seed=s64.COHORT_SEED)
    k64 = torch.tensor([s64.kill_row(p) for p in cohort])
    sets["cohort64"] = {"ids": [int(p["patient_id"]) for p in cohort], "rho": [p["rho"] for p in cohort],
                        "kill": k64, "source": "script 64 solver and kill rates (seed 42)",
                        "solver": lambda: BatchedSolver([s64.make_solver(p) for p in cohort], k64.tolist()),
                        "validate": lambda: validate_against_numpy(s64.make_solver(cohort[0]), k64[0].tolist())}

    sc = s66.lhs_scenarios()
    k60 = torch.tensor([s60.kill_row(p["alpha_sens"]) for p in sc])

    def mk60(p):
        return s60.FastPDESolver(grid_size=s60.EVAL_GRID, dt_pde=s60.DT_PDE_EVAL, rho=p["rho"], D_white=p["D_w"],
                                 alpha_sens=p["alpha_sens"], **s60.ABLATIONS[s60.FULL])
    sets["lhs60"] = {"ids": list(range(len(sc))), "rho": [p["rho"] for p in sc], "kill": k60,
                     "source": "script 60 full-model solver and alpha_sens kill rates",
                     "solver": lambda: BatchedSolver([mk60(p) for p in sc], k60.tolist()),
                     "validate": lambda: validate_against_numpy(mk60(sc[0]), k60[0].tolist())}

    rt = s66.real_cohort()["test"]
    krt = s66.KILL.view(1, 4).repeat(len(rt), 1)
    sets["real_test"] = {"ids": [p["id"] for p in rt], "rho": [p["rho"] for p in rt], "kill": krt,
                         "source": "script 66 real test patients, fixed kill table",
                         "solver": lambda: KillPDE([p["rho"] for p in rt], [p["D_w"] for p in rt], krt),
                         "validate": None}

    st = pools()["synth_test"]
    sets["synth_test"] = {"ids": list(range(len(st["rho"]))), "rho": st["rho"].tolist(), "kill": st["kill"],
                          "source": f"held-out randomised patients (seed {SEED_TEST}), script-59 physics",
                          "solver": lambda: KillPDE(st["rho"].tolist(), [s59.D_WHITE_BASE] * len(st["rho"]),
                                                    st["kill"]),
                          "validate": None}
    return sets


def arm_stats(res: Dict, ref: Dict, oracle: Dict, rng: np.random.Generator) -> Dict:
    f, fs, fo = res["final"].numpy(), ref["final"].numpy(), oracle["final"].numpy()
    lr = np.log(f / fs)
    boot = np.array([lr[rng.integers(0, len(lr), len(lr))].mean() for _ in range(2000)])
    out = {"mean_final_volume_mm3": float(f.mean()), "mean_drug_auc": float(res["auc"].mean()),
           "max_drug_auc": float(res["auc"].max()),
           "win_rate_vs_stupp_pct": float((f < fs).mean() * 100),
           "mean_log_ratio_final_vs_stupp": float(lr.mean()),
           "mean_log_ratio_final_vs_stupp_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
           "mean_log_ratio_final_vs_oracle": float(np.log(f / fo).mean()),
           "mean_log_ratio_burden_vs_stupp": float(torch.log(res["mean_burden"] / ref["mean_burden"]).mean()),
           "wilcoxon_p_vs_stupp": None, "final_volume_mm3": f.tolist()}
    if np.any(lr != 0):
        out["wilcoxon_p_vs_stupp"] = float(stats.wilcoxon(np.log(f), np.log(fs)).pvalue)
    return out


def stage_evaluate():
    nets = {name: load_run(name) for name in RUNS}
    out = {"budget": BUDGET, "physics": "full 64^3 PDE incl. diffusion; each set uses its own solver",
           "kill_ranges_train": KILL_RANGES, "sets": {}}
    rng = np.random.default_rng(67)
    for set_name, st in eval_sets().items():
        t0 = time.time()
        val_err = st["validate"]() if st["validate"] else None
        if val_err is not None:
            assert val_err < 1e-9, f"{set_name}: batched solver disagrees with numpy solver ({val_err})"
        kill = st["kill"]
        oracle_sched = kill_oracle(torch.tensor(st["rho"]), kill)["sched"]
        arms = {"stupp": s66.schedule_policy(STUPP_SCHEDULE), "heuristic59": s66.heuristic59_policy,
                "ppo66": s66.greedy_policy(s66.load_net("ppo_final_s0")),
                "dagger66": s66.greedy_policy(s66.load_net("bc_dagger_final")),
                "efficiency_rule": s66.sched_tensor_policy(efficiency_rule_sched(kill)),
                "oracle": s66.sched_tensor_policy(oracle_sched)}
        for name, net in nets.items():
            arms[name] = greedy(net, kill, "blind" not in name)
        solver = st["solver"]()
        res = {name: s66.rollout(solver, pol) for name, pol in arms.items()}
        block = {name: arm_stats(r, res["stupp"], res["oracle"], rng) for name, r in res.items()}
        eff_a = (efficiency(kill).argmax(1) + 1).tolist()
        out["sets"][set_name] = {
            "source": st["source"], "n": len(st["rho"]), "ids": st["ids"], "rho": st["rho"],
            "kill": kill.tolist(), "most_efficient_action": [ACTION_NAMES[a] for a in eff_a],
            "validation_max_rel_err_vs_numpy_solver": val_err, "arms": block,
            "actions": {name: r["actions"].tolist() for name, r in res.items()}}
        print(f"\n  [{set_name}] n={len(st['rho'])}  ({time.time() - t0:.0f}s, validation {val_err})")
        for name, s in block.items():
            print(f"    {name:16s} mean {s['mean_final_volume_mm3']:8.2f}  win {s['win_rate_vs_stupp_pct']:5.1f}%  "
                  f"vs Stupp {s['mean_log_ratio_final_vs_stupp']:+.4f} "
                  f"[{s['mean_log_ratio_final_vs_stupp_ci95'][0]:+.3f}, {s['mean_log_ratio_final_vs_stupp_ci95'][1]:+.3f}]"
                  f"  vs oracle {s['mean_log_ratio_final_vs_oracle']:+.4f}  AUC max {s['max_drug_auc']:.2f}  "
                  f"p {s['wilcoxon_p_vs_stupp']}")
        save_json(out, "evaluate.json")


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def stage_report():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ev = json.loads((OUTPUT_DIR / "evaluate.json").read_text())
    tr = json.loads((OUTPUT_DIR / "train.json").read_text())
    arms = ["heuristic59", "ppo66", "dagger66", "efficiency_rule", "dagger_blind", "ppo_blind_ft",
            "ppo_cond", "dagger_cond", "ppo_cond_ft", "oracle"]
    sets = list(ev["sets"])
    fig, axes = plt.subplots(2, 2, figsize=(17, 12))

    ax = axes[0, 0]
    w = 0.8 / len(sets)
    for j, s in enumerate(sets):
        m = [ev["sets"][s]["arms"][a]["mean_log_ratio_final_vs_stupp"] for a in arms]
        ci = np.array([ev["sets"][s]["arms"][a]["mean_log_ratio_final_vs_stupp_ci95"] for a in arms])
        x = np.arange(len(arms)) + (j - (len(sets) - 1) / 2) * w
        ax.bar(x, m, w, label=f"{s} (n={ev['sets'][s]['n']})")
        ax.errorbar(x, m, yerr=[np.array(m) - ci[:, 0], ci[:, 1] - np.array(m)], fmt="none", ecolor="k", lw=0.8)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels(arms, rotation=35, ha="right")
    ax.set_ylabel("mean log(V90 / V90,Stupp)   (< 0 beats Stupp)")
    ax.set_title("A. Equal drug budget: effect vs Stupp (bootstrap 95% CI)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[0, 1]
    for j, s in enumerate(sets):
        g = [ev["sets"][s]["arms"][a]["mean_log_ratio_final_vs_oracle"] for a in arms[:-1]]
        ax.bar(np.arange(len(arms) - 1) + (j - (len(sets) - 1) / 2) * w, g, w, label=s)
    ax.set_xticks(range(len(arms) - 1))
    ax.set_xticklabels(arms[:-1], rotation=35, ha="right")
    ax.set_ylabel("mean log(V90 / V90,oracle)   (0 = oracle)")
    ax.set_title("B. Gap to the privileged per-patient oracle")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1, 0]
    s = ev["sets"]["cohort64"]
    st = np.array(s["arms"]["stupp"]["final_volume_mm3"])
    cols = {"chemo": "#fdae61", "rad": "#abd9e9", "combo": "#2c7bb6"}
    for a, mk in (("ppo66", "x"), ("dagger_cond", "o"), ("oracle", "s")):
        lr = np.log(np.array(s["arms"][a]["final_volume_mm3"]) / st)
        ax.scatter(np.array(s["rho"]), lr, marker=mk, s=60,
                   c=[cols[e] for e in s["most_efficient_action"]], edgecolors="k", label=a)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("rho (1/day)")
    ax.set_ylabel("log(V90 / V90,Stupp)")
    ax.set_title("C. Script-64 cohort per patient (colour = most drug-efficient action)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    for name, run in tr["runs"].items():
        if "curve" in run:
            c = run["curve"]
            ax.plot([r["iter"] for r in c], [r["val"]["final_log_ratio_vs_oracle"] for r in c], label=name)
        else:
            d = run["dagger"]
            ax.plot([r["round"] * 10 for r in d], [r["val"]["final_log_ratio_vs_oracle"] for r in d], "--o",
                    label=f"{name} (x = round*10)")
    for name, r in tr["val_reference"].items():
        ax.axhline(r["final_log_ratio_vs_oracle"], ls=":", lw=1, color="gray")
        ax.text(0, r["final_log_ratio_vs_oracle"], name, fontsize=7, va="bottom", color="gray")
    ax.set_xlabel("PPO iteration")
    ax.set_ylabel("validation mean log(V90 / V90,oracle)")
    ax.set_title("D. Training curves (reaction-only sim, 300 held-out patients)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    plt.suptitle(f"Script 67: policies that know patient drug sensitivity, equal drug budget (AUC {BUDGET:g})",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    path = OUTPUT_DIR / "kill_conditioned_policy.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[saved] {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["oracle", "train", "evaluate", "report"])
    ap.add_argument("--runs", nargs="*")
    args = ap.parse_args()
    print("=" * 70)
    print("SCRIPT 67: KILL-RATE-CONDITIONED POLICIES AT EQUAL DRUG BUDGET")
    print("=" * 70)
    stages = [args.stage] if args.stage else ["oracle", "train", "evaluate", "report"]
    for st in stages:
        t0 = time.time()
        print(f"\n--- stage {st} ---")
        {"oracle": stage_oracle, "train": lambda: stage_train(args.runs), "evaluate": stage_evaluate,
         "report": stage_report}[st]()
        print(f"--- stage {st} done in {time.time() - t0:.0f}s ---")


if __name__ == "__main__":
    main()
