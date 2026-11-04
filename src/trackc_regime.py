"""Helpers for the Stage 2 regime map (test E). SIMULATION only; every rule and range is an assumption.

    candidate_schedules   a small family of open-loop schedules (at most 12 courses), used for the best fixed
                          schedule and for the informed upper bound
    rmst, paired_stats    restricted mean time and paired differences between two arms on the same patients
    evaluate_schedules    outcome of every candidate schedule on every patient (open loop, so no noise)
    policy_arms           the named rules compared in test E
"""
from __future__ import annotations

import numpy as np

import trackc_policies as pol
import trackc_sim as tc

TAUS = (365.0, 730.0)
SCAN_EVERY = 3


def candidate_schedules():
    """Rows are 0/1 schedules over the 52 cycles. Front-loaded 0..12 courses; every q-th cycle from a start cycle
    (at most 12 courses); two blocks with a gap. Returned with a short name per row."""
    rows, names = [], []

    def add(name, cycles):
        mask = np.zeros(tc.N_CYCLES, bool)
        for c in cycles:
            if c < tc.N_CYCLES:
                mask[c] = True
        if mask.sum() <= tc.BUDGET and not any((mask == r).all() for r in rows):
            rows.append(mask)
            names.append(name)

    for k in range(0, tc.BUDGET + 1):
        add(f"front{k}", range(k))
    for q in (2, 3, 4, 6):
        for start in (0, 1, 2, 3):
            add(f"every{q}_from{start}", [start + q * i for i in range(tc.BUDGET)])
    for k1 in (2, 4, 6):
        for gap in (2, 4, 8):
            for k2 in (2, 4, 6):
                add(f"block{k1}_gap{gap}_block{k2}", list(range(k1)) + list(range(k1 + gap, k1 + gap + k2)))
    return np.array(rows), names


def rmst(ttp, tau):
    """Restricted mean of min(ttp, tau)."""
    return float(np.mean(np.minimum(np.asarray(ttp, float), tau)))


def paired_stats(ttp_arm, ttp_ref, tau, alpha=0.10):
    """Mean paired difference in restricted time (arm minus reference) with a normal 95% interval over patients,
    and the paired-regret tail: the mean of the lowest `alpha` share of the paired differences."""
    d = np.minimum(np.asarray(ttp_arm, float), tau) - np.minimum(np.asarray(ttp_ref, float), tau)
    se = d.std(ddof=1) / np.sqrt(d.size) if d.size > 1 else float("nan")
    k = max(1, int(np.ceil(alpha * d.size)))
    return {"diff": float(d.mean()), "ci_lo": float(d.mean() - 1.96 * se), "ci_hi": float(d.mean() + 1.96 * se),
            "regret_cvar": float(np.sort(d)[:k].mean()), "frac_worse": float(np.mean(d < 0))}


def evaluate_schedules(scen, schedules, endpoint, scan_every=SCAN_EVERY, seed=(0, 0)):
    """TTP of every schedule on every patient: array (n_schedules, n_patients). Open loop, so the seed does not matter."""
    n_s, n_p = len(schedules), len(scen["r"])
    tiled = {k: np.tile(v, n_s) for k, v in scen.items()}
    tiled["scan_every"] = np.full(n_s * n_p, int(scan_every))
    mask = np.repeat(schedules, n_p, axis=0)
    out = tc.simulate(tiled, pol.Schedule(mask), list(seed), endpoint=endpoint)
    return out["ttp"].reshape(n_s, n_p), out["courses"].reshape(n_s, n_p)


def policy_arms():
    """name -> (policy factory, course budget). All read scans every 84 days where they read at all."""
    return {
        "soc6": (lambda: tc.Periodic(1), 6),
        "front3": (lambda: tc.Periodic(1), 3),
        "front9": (lambda: tc.Periodic(1), 9),
        "front12": (lambda: tc.Periodic(1), 12),
        "paced2": (lambda: tc.Periodic(2), 12),
        "paced3": (lambda: tc.Periodic(3), 12),
        "raw_zhang_1.0_0.5": (lambda: tc.RawThreshold(1.0, 0.5), 12),
        "raw_0.6_0.6": (lambda: tc.RawThreshold(0.6, 0.6), 12),
        "nadir_k3_on1.25": (lambda: pol.NadirContainment(3, 1.25, 1.0), 12),
        "nadir_k3_on1.15": (lambda: pol.NadirContainment(3, 1.15, 1.0), 12),
        "nadir_k6_on1.25": (lambda: pol.NadirContainment(6, 1.25, 1.0), 12),
    }


def run_policy(scen, factory, budget, endpoint, seed, scan_every=SCAN_EVERY):
    s = dict(scen)
    s["scan_every"] = np.full(len(s["r"]), int(scan_every))
    out = tc.simulate(s, factory(), list(seed), endpoint=endpoint, budget=budget)
    return out["ttp"], out["courses"]
