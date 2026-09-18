#!/usr/bin/env python3
"""Trust-gated MPC-RL hybrid: 4-arm virtual trial on the real-patient cohort.

Arms (all inside the same 3D anisotropic PDE environment):
  0. Standard Stupp
  1. MPC-alone (corrected robust MPC as a policy)
  2. RL-alone (trust-conditioned policy WITHOUT the trust feature -> obs_dim 5)
  3. Hybrid-A (hard gate: |benefit| > tau -> MPC, else RL)
  4. Hybrid-B (trust-conditioned RL: obs_dim 6, trust injected as a feature)

Patient parameters (rho, D) come from output/real_patient_validation.csv
(joint treatment-aware inverse estimation on the MU-Glioma-Post cohort).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import hybrid_controller as hc  # noqa: E402
import trust_signal  # noqa: E402
from hybrid_controller import (
    MPCAsPolicy,
    HybridGateA,
    HybridGateB,
    TrustConditionedPolicy,
    train_reinforce,
)

_spec = importlib.util.spec_from_file_location("vcs", ROOT / "src" / "64_virtual_cohort_simulation.py")
vcs = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(vcs)

_spec_r = importlib.util.spec_from_file_location("resenv", ROOT / "src" / "resistance_env.py")
resenv = importlib.util.module_from_spec(_spec_r)
assert _spec_r.loader is not None
_spec_r.loader.exec_module(resenv)

FastPDESolver = vcs.FastPDESolver
GbmTherapyEnv = vcs.GbmTherapyEnv
RL_REWARD_WEIGHTS = vcs.RL_REWARD_WEIGHTS
# Reward weights that also penalize resistance selection.
# Balanced: volume control dominates, resistance penalty moderate (avoid both
# over-dosing that breeds resistance and under-dosing that lets tumor escape).
# Endpoint penalty strongly penalizes finishing with a resistant-heavy tumor.
RL_REWARD_WEIGHTS_RES = dict(RL_REWARD_WEIGHTS, lambda_res=8.0, lambda_res_end=30.0)


def make_env(rho: float, D_white: float, grid_size=(16, 16, 16), resistance: bool = True, max_steps: int = 90):
    if resistance:
        solver = resenv.ResistancePDESolver(
            grid_size=grid_size,
            dt_pde=0.2,
            rho=rho,
            D_white=D_white,
            alpha_sens=1.0,
            gamma_chemo=0.08,
            alpha_rt=0.045,
            is_training=False,
        )
        env = resenv.ResistanceGbmEnv(solver, RL_REWARD_WEIGHTS_RES, max_steps=max_steps)
    else:
        solver = FastPDESolver(
            grid_size=grid_size,
            dt_pde=0.2,
            rho=rho,
            D_white=D_white,
            alpha_sens=1.0,
            gamma_chemo=0.08,
            alpha_rt=0.045,
            is_training=False,
        )
        env = GbmTherapyEnv(solver, RL_REWARD_WEIGHTS)
    env.reset()  # initialize solver.initial_volume
    return env


def run_arm(env, action_fn, reset=True) -> Dict[str, Any]:
    """Run one episode. action_fn(env, obs) -> (action, diag)."""
    obs, _ = env.reset()
    trajectory = []
    total_reward = 0.0
    for step in range(env.max_steps):
        action, diag = action_fn(env, obs)
        obs, reward, terminated, truncated, _ = env.step(action)
        total_reward += reward
        trajectory.append({
            "day": step + 1,
            "action": int(action),
            "volume_mm3": float(env.solver.u.sum() * env.solver.dx**3),
        })
        if terminated or truncated:
            break
    return {
        "final_volume_mm3": float(env.solver.u.sum() * env.solver.dx**3),
        "total_reward": float(total_reward),
        "trajectory": trajectory,
        "dose_days": sum(1 for t in trajectory if t["action"] > 0),
        "peak_u_max": float(env.solver.u.max()),
        "resistant_fraction": float(
            getattr(env.solver, "u_r", 0.0).sum() / max(float(env.solver.u.sum()), 1e-12)
        ) if hasattr(env.solver, "u_r") else 0.0,
    }


def stupp_fn(env, obs):
    step = env.solver.step_count
    day = step + 1
    total = env.max_steps
    # RT: first 30 days (days 20-50) -> Combo
    if 20 <= day < 50:
        action = 3
    elif day >= 50:
        # TMZ 5-on/28-off cycles for the remainder
        action = 1 if (int(day) % 28) < 5 else 0
    else:
        action = 0
    return action, {"arm": "stupp"}


def run_trial(
    patients: List[Dict[str, float]],
    policy_b: Optional[TrustConditionedPolicy] = None,
    policy_rl: Optional[TrustConditionedPolicy] = None,
    tau: float = 0.05,
    grid_size=(16, 16, 16),
    trust_mode: str = "margin",
    days: int = 90,
) -> List[Dict[str, Any]]:
    results = []
    for pat in patients:
        rho = pat["rho"]
        D = pat["D"]
        env = make_env(rho, D, grid_size, max_steps=days)
        target_vol = env.solver.initial_volume * 0.12

        # Arm 0: Stupp
        r_stupp = run_arm(env, stupp_fn)

        # Arm 1: MPC-alone
        mpc = MPCAsPolicy(rho=rho, D=D, target_volume_mm3=target_vol)
        mpc.reset()
        r_mpc = run_arm(env, lambda e, o, m=mpc: m.choose_action(e, o))

        # Arm 2: RL-alone (obs_dim 5)
        def rl_fn(env, obs):
            a = hc.rl_choose_action(policy_rl, np.asarray(obs, dtype=np.float32))
            return a, {"arm": "rl"}
        r_rl = run_arm(env, rl_fn)

        # Arm 3: Hybrid-A (hard gate delegates to 5-dim RL policy)
        gate_a = HybridGateA(rho=rho, D=D, target_volume_mm3=target_vol, policy=policy_rl, tau=tau, trust_mode=trust_mode)
        gate_a.reset()
        r_a = run_arm(env, lambda e, o, g=gate_a: g.choose_action(e, o))

        # Arm 4: Hybrid-B (trust-conditioned 6-dim policy)
        gate_b = HybridGateB(rho=rho, D=D, target_volume_mm3=target_vol, policy=policy_b, trust_mode=trust_mode)
        gate_b.reset()
        r_b = run_arm(env, lambda e, o, g=gate_b: g.choose_action(e, o))

        results.append({
            "patient_id": pat["patient_id"],
            "rho": rho,
            "D": D,
            "stupp": r_stupp,
            "mpc": r_mpc,
            "rl": r_rl,
            "hybrid_a": r_a,
            "hybrid_b": r_b,
            "target_volume_mm3": target_vol,
        })
    return results


def summarize(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    arms = ["stupp", "mpc", "rl", "hybrid_a", "hybrid_b"]
    out: Dict[str, Any] = {"n_patients": len(results)}
    for arm in arms:
        vols = np.array([r[arm]["final_volume_mm3"] for r in results])
        dose_days = np.array([r[arm]["dose_days"] for r in results])
        res_frac = np.array([r[arm]["resistant_fraction"] for r in results])
        out[arm] = {
            "mean_final_volume_mm3": float(np.mean(vols)),
            "median_final_volume_mm3": float(np.median(vols)),
            "std_final_volume_mm3": float(np.std(vols)),
            "mean_dose_days": float(np.mean(dose_days)),
            "mean_resistant_fraction": float(np.mean(res_frac)),
        }
    stupp_dose = out["stupp"]["mean_dose_days"]
    for arm in arms:
        out[arm]["dose_sparing_vs_stupp"] = (
            float(1.0 - out[arm]["mean_dose_days"] / stupp_dose) if stupp_dose > 0 else 0.0
        )
    # Paired statistics (hybrid_b vs each)
    from scipy import stats
    for base in ["stupp", "mpc", "rl", "hybrid_a"]:
        b_vols = np.array([r["hybrid_b"]["final_volume_mm3"] for r in results])
        base_vols = np.array([r[base]["final_volume_mm3"] for r in results])
        if len(b_vols) > 1 and np.std(b_vols - base_vols) > 0:
            t, p = stats.ttest_rel(b_vols, base_vols)
            w, wp = stats.wilcoxon(b_vols, base_vols)
            out[f"hybrid_b_vs_{base}"] = {
                "mean_diff_mm3": float(np.mean(b_vols - base_vols)),
                "paired_t_p": float(p),
                "wilcoxon_p": float(wp),
            }
    return out


def train_policy_curriculum(
    obs_dim: int,
    rho_values=(0.005, 0.015, 0.025, 0.035, 0.05),
    episodes_per_rho: int = 20,
    grid_size=(16, 16, 16),
    trust_mode: str = "margin",
    days: int = 90,
) -> Optional[TrustConditionedPolicy]:
    """Train a policy across a rho curriculum so it generalizes to the cohort."""
    policy = None
    for rho in rho_values:
        env = make_env(rho, 0.0012, grid_size, max_steps=days)
        if obs_dim == 6:
            from hybrid_controller import TrustObsWrapper, make_trust_fn, RobustMPCController, CalibrationTrust
            target = env.solver.initial_volume * 0.12
            ctrl = RobustMPCController(seed=0)
            if trust_mode == "calibration":
                cal = CalibrationTrust(window=7, lambda_cal=3.0)
                mpc_pred = MPCAsPolicy(rho=rho, D=0.0012, target_volume_mm3=target, controller=ctrl)
                mpc_pred.reset()
                pending: list = []

                def tr_fn(e=env, cal=cal, mpc=mpc_pred):
                    vol = float(e.solver.u.sum() * e.solver.dx ** 3)
                    # MPC predicts the horizon from the current volume
                    pred = trust_signal._predict_volume_corrected(
                        vol, rho, 0.0012, mpc.controller.horizon, mpc.step,
                        drug_on=True, target_volume=target,
                        w_tumor=mpc.controller.w_tumor, w_drug=mpc.controller.w_drug,
                    )
                    pending.append((pred, vol))
                    if len(pending) > 1:
                        old_pred, old_vol = pending.pop(0)
                        cal.record(old_pred, vol)  # compare prior prediction vs current reality
                    return cal.trust()
            else:
                tr_fn = make_trust_fn(
                    ctrl, rho, 0.0012, target,
                    lambda e=env: float(e.solver.u.sum() * e.solver.dx ** 3),
                )
            env = TrustObsWrapper(env, tr_fn)
        policy = train_reinforce(env, obs_dim=obs_dim, episodes=episodes_per_rho)
    return policy


def main() -> int:
    parser = argparse.ArgumentParser(description="Trust-gated MPC-RL hybrid virtual trial")
    parser.add_argument("--params", type=Path, default=ROOT / "output/real_patient_validation.csv")
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--tau", type=float, default=0.05)
    parser.add_argument("--grid", type=int, default=16)
    parser.add_argument("--train-episodes", type=int, default=40)
    parser.add_argument("--trust-mode", choices=["margin", "calibration"], default="margin")
    parser.add_argument("--no-resistance", action="store_true")
    parser.add_argument("--days", type=int, default=90, help="Episode horizon in days (RL steps)")
    parser.add_argument("--output", type=Path, default=ROOT / "output/hybrid_cohort_metrics.json")
    args = parser.parse_args()

    df = pd.read_csv(args.params)
    df = df.head(args.limit)
    patients = [
        {"patient_id": str(r["patient_id"]), "rho": float(r["rho"]), "D": float(r["D"])}
        for _, r in df.iterrows()
    ]

    grid = (args.grid, args.grid, args.grid)

    # Train policies across a rho curriculum so the arms share the same
    # learned policy and generalize to the real-parameter cohort.
    print("[Trial] Training trust-conditioned policy (obs_dim=6) ...")
    policy_b = train_policy_curriculum(
        obs_dim=6, episodes_per_rho=args.train_episodes, grid_size=grid,
        trust_mode=args.trust_mode, days=args.days,
    )

    print("[Trial] Training RL-alone policy (obs_dim=5) ...")
    policy_rl = train_policy_curriculum(
        obs_dim=5, episodes_per_rho=args.train_episodes, grid_size=grid,
        trust_mode=args.trust_mode, days=args.days,
    )

    print(f"[Trial] Running {len(patients)} patients x 5 arms on {grid} grid "
          f"(trust_mode={args.trust_mode}, resistance={not args.no_resistance}, days={args.days}) ...")
    results = run_trial(
        patients, policy_b=policy_b, policy_rl=policy_rl, tau=args.tau,
        grid_size=grid, trust_mode=args.trust_mode, days=args.days,
    )

    summary = summarize(results)

    # Lightweight per-patient summary (omit trajectories for size)
    light = []
    for r in results:
        light.append({
            "patient_id": r["patient_id"], "rho": r["rho"], "D": r["D"],
            **{a: {k: v for k, v in r[a].items() if k != "trajectory"} for a in
               ["stupp", "mpc", "rl", "hybrid_a", "hybrid_b"]},
        })

    payload = {"summary": summary, "patients": light}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"\nSaved -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())