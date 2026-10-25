"""
Track C oracle with known parameters (Claude plan trackAC_new_results_plan_claude.md, result 1).
==============================================================================================
SIMULATION ONLY. The oracle knows the true scenario parameters and state, so it is an
upper reference, not a usable policy.

The model in trackc_sim is autonomous (no time dependence), so one cycle maps (S, R) to the
next (S, R) the same way at every cycle. Value iteration on a log-spaced (S, R) grid:

    V_k(x, b) = max over a in {0, 1 if b > 0} of
                hit time inside the cycle            if the cycle reaches 1.2 B0
                28 + V_{k-1}(x', b - a)              otherwise
    V_0 = 0,  k = cycles left,  b = courses left.

V_{k-1}(x') is bilinear in (log S, log R). The grid value is only used to choose actions:
`oracle_schedule` rolls the TRUE state forward and, at each cycle, picks the action with the
higher one-step exact cycle + interpolated value. The returned TTP is from an exact rollout of
that schedule, so it is achievable (a lower bound on the true optimum). `best_known` also takes
the max with every periodic schedule, so it is never below them.
"""
from __future__ import annotations

import itertools

import numpy as np

import trackc_sim as tc

LOG_MIN = np.log(1e-12)


def _scalar(scen, i):
    return {key: np.asarray(val[i : i + 1], float) for key, val in scen.items()}


def _grid(p, level, G, budget):
    """Per-axis log grids from just below the lowest reachable value (all courses kill,
    no regrowth) up to the progression level. A tight range keeps the cell size small."""
    k5 = tc.COURSE_DAYS * p["m"][0] * p["r"][0]
    S0, R0 = (1 - p["f0"][0]) * p["B0"][0], p["f0"][0] * p["B0"][0]
    lo_s = max(np.log(S0) - k5 * budget - 1.0, LOG_MIN)
    lo_r = max(np.log(R0) - p["eps"][0] * k5 * budget - 1.0, LOG_MIN)
    axs, axr = np.linspace(lo_s, np.log(level), G), np.linspace(lo_r, np.log(level), G)
    LS, LR = np.meshgrid(axs, axr, indexing="ij")
    return (axs, axr), np.exp(LS.ravel()), np.exp(LR.ravel())


def _interp(V, ax, S, R):
    """Bilinear interpolation of V[..., G, G] at (log S, log R), clamped to the grid."""
    def coord(a, x):
        u = (np.clip(np.log(np.maximum(x, 1e-300)), a[0], a[-1]) - a[0]) / (a[1] - a[0])
        i = np.minimum(np.floor(u).astype(int), a.size - 2)
        return i, u - i

    i, fi = coord(ax[0], S)
    j, fj = coord(ax[1], R)
    return ((1 - fi) * (1 - fj) * V[..., i, j] + fi * (1 - fj) * V[..., i + 1, j]
            + (1 - fi) * fj * V[..., i, j + 1] + fi * fj * V[..., i + 1, j + 1])


def _step(S, R, p, a, level):
    """One exact cycle from t = 0. Returns next S, R and hit time (NaN = no progression)."""
    ttp = np.full(np.shape(S), np.nan)
    S1, R1 = tc.run_cycle(S, R, p, np.full(np.shape(S), float(a)), 0.0,
                          prog_level=level, ttp=ttp)
    return S1, R1, ttp


def value_table(p, n_cycles=tc.N_CYCLES, budget=tc.BUDGET, G=80):
    """V[k, b, i, j] for k = 0..n_cycles cycles left and b = 0..budget courses left."""
    level = float(tc.PROG_FACTOR * p["B0"][0])
    ax, S, R = _grid(p, level, G, budget)
    alive = (S + R) < level
    nxt = [_step(S, R, p, a, level) for a in (0, 1)]
    V = np.zeros((n_cycles + 1, budget + 1, G, G))
    for k in range(1, n_cycles + 1):
        prev = V[k - 1]
        q0 = np.where(np.isnan(nxt[0][2]), tc.CYCLE_DAYS + _interp(prev, ax, nxt[0][0], nxt[0][1]),
                      nxt[0][2])
        q1 = np.where(np.isnan(nxt[1][2]), tc.CYCLE_DAYS + _interp(prev[:-1], ax, nxt[1][0], nxt[1][1]),
                      nxt[1][2])
        best = q0.copy()
        best[1:] = np.maximum(q0[1:], q1)
        V[k] = (best * alive).reshape(budget + 1, G, G)
    return V, ax


def rollout(p, schedules, n_cycles=tc.N_CYCLES):
    """Exact capped TTP for each 0/1 schedule row (no noise: TTP depends on the truth only)."""
    sched = np.asarray(schedules, bool)
    n = sched.shape[0]
    level = tc.PROG_FACTOR * p["B0"][0]
    S = np.full(n, (1 - p["f0"][0]) * p["B0"][0])
    R = np.full(n, p["f0"][0] * p["B0"][0])
    ttp = np.full(n, np.nan)
    for j in range(n_cycles):
        treat = sched[:, j] & np.isnan(ttp)
        S, R = tc.run_cycle(S, R, p, treat, j * tc.CYCLE_DAYS, prog_level=level, ttp=ttp)
        if not np.isnan(ttp).any():
            break
    return np.where(np.isnan(ttp), n_cycles * tc.CYCLE_DAYS, ttp)


def oracle_schedule(p, n_cycles=tc.N_CYCLES, budget=tc.BUDGET, G=80):
    """Greedy exact-state rollout against the value table. Returns (schedule, TTP, V estimate)."""
    V, ax = value_table(p, n_cycles, budget, G)
    level = tc.PROG_FACTOR * p["B0"][0]
    S = np.array([(1 - p["f0"][0]) * p["B0"][0]])
    R = np.array([p["f0"][0] * p["B0"][0]])
    v_est = float(_interp(V[n_cycles, budget], ax, S, R)[0])
    sched = np.zeros(n_cycles, bool)
    b = budget
    for j in range(n_cycles):
        k = n_cycles - j
        q = []
        for a in ((0, 1) if b > 0 else (0,)):
            S1, R1, hit = _step(S, R, p, a, level)
            q.append(hit[0] if not np.isnan(hit[0]) else tc.CYCLE_DAYS + _interp(V[k - 1, b - a], ax, S1, R1)[0])
        a = int(np.argmax(q))
        sched[j] = bool(a)
        S, R, hit = _step(S, R, p, a, level)
        b -= a
        if not np.isnan(hit[0]):
            break
    return sched, float(rollout(p, sched[None], n_cycles)[0]), v_est


def periodic_schedules(n_cycles=tc.N_CYCLES, budget=tc.BUDGET, qs=range(1, 7)):
    rows = []
    for q in qs:
        s = np.zeros(n_cycles, bool)
        s[np.arange(0, n_cycles, q)[:budget]] = True
        rows.append(s)
    return np.array(rows)


def _neighbours(s, budget):
    """Every schedule one exact edit away: move one course, add one, or drop one."""
    on, off = np.flatnonzero(s), np.flatnonzero(~s)
    rows = []
    for i in on:
        t = s.copy()
        t[i] = False
        rows.append(t)
        for j in off:
            u = t.copy()
            u[j] = True
            rows.append(u)
    if on.size < budget:
        for j in off:
            u = s.copy()
            u[j] = True
            rows.append(u)
    return np.array(rows)


def polish(p, sched, budget=tc.BUDGET, n_cycles=tc.N_CYCLES, max_iter=50):
    """Hill-climb on exact rollouts. The grid value has jumps (hit times are discrete), so
    interpolation alone can miss the optimum by a few days; this fixes the misses found in
    the brute-force tests. It can only raise the TTP."""
    best_s = sched.copy()
    best_t = float(rollout(p, best_s[None], n_cycles)[0])
    for _ in range(max_iter):
        nb = _neighbours(best_s, budget)
        t = rollout(p, nb, n_cycles)
        i = int(np.argmax(t))
        if t[i] <= best_t + 1e-9:
            break
        best_s, best_t = nb[i], float(t[i])
    return best_s, best_t


def best_known(p, n_cycles=tc.N_CYCLES, budget=tc.BUDGET, G=80):
    """Polished best of (DP schedule, periodic q = 1..6). Every value is an exact rollout,
    so this is an achievable lower bound on the true known-parameter optimum."""
    sched, ttp_dp, v_est = oracle_schedule(p, n_cycles, budget, G)
    per_rows = periodic_schedules(n_cycles, budget)
    per = rollout(p, per_rows, n_cycles)
    start = sched if ttp_dp >= per.max() else per_rows[int(np.argmax(per))]
    s, t = polish(p, start, budget, n_cycles)
    return {"ttp": t, "ttp_dp": ttp_dp, "v_est": v_est,
            "ttp_best_periodic": float(per.max()), "schedule": s}


def brute_force(p, n_cycles, budget):
    """Exact optimum over every schedule with at most `budget` courses (small cases only)."""
    rows = [np.isin(np.arange(n_cycles), c) for m in range(budget + 1)
            for c in itertools.combinations(range(n_cycles), m)]
    return float(rollout(p, np.array(rows), n_cycles).max())
