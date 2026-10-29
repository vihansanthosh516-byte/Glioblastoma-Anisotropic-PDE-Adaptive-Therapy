"""
Track C simulator, particle filter and dosing policies (Claude final plan, appendix Track C).
==========================================================================================
SIMULATION ONLY. Every parameter range is an assumption, not a patient or drug measurement.

Two-population model (sensitive S, resistant R, K = 1), with optional asymmetric
competition used only by the misspecified / positive-control suites:

    dS/dt = r S (1 - (S + a_rs R)) - k u S,   k = m r
    dR/dt = r (1 - c) R (1 - (a_sr S + R)) - eps k u R

Treatment: 28-day cycles; a treated cycle gets u = 1 on days 0-4. Budget 12 courses,
horizon 52 cycles = 1,456 days. One noisy reading per cycle start:
y_j = log B(t_j) + e_j, e_j ~ N(0, sigma^2), sigma known to the filter.
Simulated progression: first time the TRUE burden B = S + R reaches 1.2 B0 (not RANO).
With B0 <= 0.5 and K = 1 this level (<= 0.6) is reachable without treatment.

Particle filter: 500 particles over (log r, log k, log f0, logit c, logit eps, logit B0),
drawn from the generating priors; Gaussian likelihood; systematic resampling when
ESS < P/2, then one Liu-West move (h = 0.90, chosen on tuning seed 555 for interval coverage) and deterministic re-propagation of S, R
from t = 0 under the known past actions (no process noise).

Policies decide once per cycle: periodic schedules, raw-reading on/off thresholds
(Zhang et al. 2017 style), and the filtered on/off rule with a resistance mode.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import qmc

CYCLE_DAYS = 28
COURSE_DAYS = 5
N_CYCLES = 52
HORIZON = N_CYCLES * CYCLE_DAYS          # 1,456 days
BUDGET = 12
PROG_FACTOR = 1.2
DT = 0.25                                 # RK4 step; tests check it against 0.05 (rel. error < 1e-4)

PARAMS = ("r", "m", "f0", "c", "eps", "B0")      # kill rate k = m * r
# m = k / r >= 6 keeps every scenario treatable: one 5-day course outweighs 28 days of growth
# (pilot 1 showed k ~ r scenarios progress in ~13 days under every policy). Assumed, not measured.
RANGES = {"r": (0.001, 0.05, "log"), "m": (6.0, 20.0, "log"), "f0": (0.001, 0.1, "log"),
          "c": (0.0, 0.25, "lin"), "eps": (0.0, 0.2, "lin"), "B0": (0.05, 0.5, "lin")}
SIGMAS = (0.05, 0.10, 0.20)
RANGES_WIDE_M = dict(RANGES, m=(3.0, 20.0, "log"))     # filter prior that can hold a weak drug


# ---------------------------------------------------------------- scenarios
def draw_scenarios(n, seed, ranges=RANGES, sigmas=SIGMAS, a_sr=1.0, a_rs=1.0):
    """Latin-hypercube scenarios; each gets one sigma, balanced across the list."""
    rng = np.random.default_rng(seed)
    u = qmc.LatinHypercube(d=len(PARAMS), seed=rng).random(n)
    out = {}
    for i, name in enumerate(PARAMS):
        lo, hi, kind = ranges[name]
        if kind == "log":
            out[name] = np.exp(np.log(lo) + u[:, i] * (np.log(hi) - np.log(lo)))
        else:
            out[name] = lo + u[:, i] * (hi - lo)
    out["sigma"] = rng.permutation(np.resize(np.asarray(sigmas, float), n))
    out["a_sr"] = np.broadcast_to(np.asarray(a_sr, float), (n,)).copy()
    out["a_rs"] = np.broadcast_to(np.asarray(a_rs, float), (n,)).copy()
    return out


def subset(scen, idx):
    return {key: val[idx] for key, val in scen.items()}


# ---------------------------------------------------------------- dynamics
def _deriv(S, R, p, u):
    dS = p["r"] * S * (1.0 - (S + p["a_rs"] * R)) - p["m"] * p["r"] * u * S
    dR = p["r"] * (1.0 - p["c"]) * R * (1.0 - (p["a_sr"] * S + R)) - p["eps"] * p["m"] * p["r"] * u * R
    if "mu" in p:                     # script 169: drug-induced switching S -> R at rate mu (truth only)
        dS, dR = dS - p["mu"] * u * S, dR + p["mu"] * u * S
    return dS, dR


def _rk4(S, R, p, u, dt):
    k1s, k1r = _deriv(S, R, p, u)
    k2s, k2r = _deriv(S + 0.5 * dt * k1s, R + 0.5 * dt * k1r, p, u)
    k3s, k3r = _deriv(S + 0.5 * dt * k2s, R + 0.5 * dt * k2r, p, u)
    k4s, k4r = _deriv(S + dt * k3s, R + dt * k3r, p, u)
    return (S + dt / 6.0 * (k1s + 2 * k2s + 2 * k3s + k4s),
            R + dt / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r))


# Event-time mode (2026-10-26, audit C1). "refine": a step that crosses the progression level
# is re-solved by bisection on a single RK4 sub-step from the start of that step, so the
# reported time is the RK4 crossing time, not the next grid point. "grid": the old rule
# (first grid point at or above the level), kept only to reproduce outputs of scripts 150-158.
EVENT_MODE = "refine"
_BISECT_ITERS = 40


def _segments(dt):
    """Steps for the on (0-5 d) and off (5-28 d) parts. Each part gets a whole number of equal
    steps no longer than dt, so the treatment switch always falls on a step boundary. For
    dt = 0.25 this is 20 + 92 steps, identical to the earlier single-grid loop."""
    n_on = int(np.ceil(COURSE_DAYS / dt - 1e-9))
    n_off = int(np.ceil((CYCLE_DAYS - COURSE_DAYS) / dt - 1e-9))
    return [(n_on, COURSE_DAYS / n_on, True), (n_off, (CYCLE_DAYS - COURSE_DAYS) / n_off, False)]


def _crossing_time(S, R, p, u, h, level):
    """Bisection for tau in (0, h] with B(tau) = level, B from one RK4 step of size tau."""
    lo, hi = np.zeros_like(h), h.copy()
    for _ in range(_BISECT_ITERS):
        mid = 0.5 * (lo + hi)
        s, r = _rk4(S, R, p, u, mid)
        above = (s + r) >= level
        hi = np.where(above, mid, hi)
        lo = np.where(above, lo, mid)
    return hi


def _take(x, m):
    x = np.asarray(x)
    return x[m] if x.ndim and x.shape[0] == m.shape[0] else x


def run_cycle(S, R, p, treat, t0, dt=None, prog_level=None, ttp=None, snaps=None):
    """Integrate one 28-day cycle. `treat` broadcasts against S. If `prog_level` is given,
    record in `ttp` (NaN = not yet) the first time S + R reaches prog_level (see EVENT_MODE).
    dt=None uses the module value DT at call time. `snaps` {day offset: None} is filled with
    (S, R) at those offsets into the cycle (offsets must fall on a step boundary)."""
    dt = DT if dt is None else dt
    treat = np.asarray(treat, float)
    t = t0
    for n, h, on in _segments(dt):
        u = treat if on else 0.0
        for _ in range(n):
            S1, R1 = _rk4(S, R, p, u, h)
            if prog_level is not None:
                hit = np.isnan(ttp) & (S1 + R1 >= prog_level)
                if hit.any():
                    if EVENT_MODE == "grid":
                        ttp[hit] = t + h
                    else:
                        S0b, R0b = np.broadcast_arrays(S, R, ttp)[:2]
                        pp = {k: _take(np.broadcast_to(v, ttp.shape) if np.ndim(v) else v, hit)
                              for k, v in p.items()}
                        ub = _take(np.broadcast_to(u, ttp.shape), hit) if np.ndim(u) else u
                        lv = _take(np.broadcast_to(prog_level, ttp.shape), hit)
                        tau = _crossing_time(S0b[hit], R0b[hit], pp, ub, np.full(int(hit.sum()), h), lv)
                        ttp[hit] = t + tau
            S, R = S1, R1
            t += h
            if snaps is not None:
                for d in snaps:
                    if abs(t - t0 - d) < 1e-6:
                        snaps[d] = (S, R)
    if snaps is not None and any(v is None for v in snaps.values()):
        raise ValueError("snapshot offset not on a step boundary")
    return S, R


# ---------------------------------------------------------------- particle filter
def _to_z(theta, ranges=RANGES):
    z = np.empty(theta[PARAMS[0]].shape + (len(PARAMS),))
    for i, name in enumerate(PARAMS):
        lo, hi, kind = ranges[name]
        x = theta[name]
        if kind == "log":
            z[..., i] = np.log(x)
        else:
            q = np.clip((x - lo) / (hi - lo), 1e-9, 1 - 1e-9)
            z[..., i] = np.log(q / (1 - q))
    return z


def _from_z(z, ranges=RANGES):
    theta = {}
    for i, name in enumerate(PARAMS):
        lo, hi, kind = ranges[name]
        if kind == "log":
            a, b = np.log(lo), np.log(hi)
            w = b - a
            t = np.mod(z[..., i] - a, 2 * w)          # reflect at prior bounds
            z[..., i] = a + np.where(t > w, 2 * w - t, t)
            theta[name] = np.exp(z[..., i])
        else:
            z[..., i] = np.clip(z[..., i], -20, 20)
            theta[name] = lo + (hi - lo) / (1 + np.exp(-z[..., i]))
    return theta


class ParticleFilter:
    """Batch of independent filters, one per scenario (axis 0), P particles (axis 1)."""

    def __init__(self, sigma, n_particles, rng, h=0.90, dt=None, ranges=RANGES):
        self.N, self.P, self.rng, self.h = len(sigma), n_particles, rng, h
        self.dt = DT if dt is None else dt
        self.ranges = ranges
        self.sigma = np.asarray(sigma, float)[:, None]
        u = rng.random((self.N, self.P, len(PARAMS)))
        theta = {}
        for i, name in enumerate(PARAMS):
            lo, hi, kind = ranges[name]
            theta[name] = (np.exp(np.log(lo) + u[..., i] * (np.log(hi) - np.log(lo)))
                           if kind == "log" else lo + u[..., i] * (hi - lo))
        self.z = _to_z(theta, ranges)
        self._set_theta(theta)
        self.S = (1 - theta["f0"]) * theta["B0"]
        self.R = theta["f0"] * theta["B0"]
        self.logw = np.zeros((self.N, self.P))
        self.actions = np.zeros((self.N, N_CYCLES), bool)
        self.n_resample = 0

    def _set_theta(self, theta):
        self.theta = dict(theta, a_sr=np.ones_like(theta["r"]), a_rs=np.ones_like(theta["r"]))

    @property
    def w(self):
        lw = self.logw - self.logw.max(1, keepdims=True)
        w = np.exp(lw)
        return w / w.sum(1, keepdims=True)

    def update(self, j, y):
        B = np.maximum(self.S + self.R, 1e-300)
        seen = ~np.isnan(y)[:, None]            # NaN = known missing scan (script 170): no update
        self.logw += np.where(seen, -0.5 * ((np.where(seen, y[:, None], 0.0) - np.log(B)) / self.sigma) ** 2, 0.0)
        w = self.w
        ess = 1.0 / (w ** 2).sum(1)
        idx = np.flatnonzero(ess < self.P / 2)
        if idx.size:
            self._resample_move(idx, j, w[idx])

    def _resample_move(self, idx, j, w):
        n, P, d = idx.size, self.P, len(PARAMS)
        z = self.z[idx]
        zbar = (w[..., None] * z).sum(1)
        dz = z - zbar[:, None, :]
        V = np.einsum("np,npi,npj->nij", w, dz, dz) + 1e-10 * np.eye(d)
        # systematic resampling
        pos = (self.rng.random((n, 1)) + np.arange(P)[None, :]) / P
        cw = np.cumsum(w, 1)
        cw[:, -1] = 1.0
        pick = np.array([np.searchsorted(cw[i], pos[i]) for i in range(n)])
        z_res = np.take_along_axis(z, pick[..., None], 1)
        L = np.linalg.cholesky((1 - self.h ** 2) * V)
        noise = np.einsum("nij,npj->npi", L, self.rng.standard_normal((n, P, d)))
        z_new = self.h * z_res + (1 - self.h) * zbar[:, None, :] + noise
        theta_new = _from_z(z_new, self.ranges)
        self.z[idx] = z_new
        for name in PARAMS:
            self.theta[name][idx] = theta_new[name]
        # deterministic re-propagation from t = 0 under the known actions
        p = {key: val[idx] for key, val in self.theta.items()}
        S = (1 - p["f0"]) * p["B0"]
        R = p["f0"] * p["B0"]
        for c in range(j):
            S, R = run_cycle(S, R, p, self.actions[idx, c][:, None], c * CYCLE_DAYS, self.dt)
        self.S[idx], self.R[idx] = S, R
        self.logw[idx] = 0.0
        self.n_resample += n

    def propagate(self, j, treat):
        self.actions[:, j] = treat
        self.S, self.R = run_cycle(self.S, self.R, self.theta, treat[:, None], j * CYCLE_DAYS, self.dt)

    def burden_mean(self):
        return (self.w * (self.S + self.R)).sum(1)

    def fr(self):
        return self.R / np.maximum(self.S + self.R, 1e-300)

    def prob_fr_ge(self, b):
        return (self.w * (self.fr() >= b)).sum(1)

    def fr_interval(self, lo=0.05, hi=0.95):
        """Weighted equal-tailed interval of the resistant fraction, per scenario."""
        fr, w = self.fr(), self.w
        order = np.argsort(fr, 1)
        fr_s = np.take_along_axis(fr, order, 1)
        cw = np.cumsum(np.take_along_axis(w, order, 1), 1)
        rows = np.arange(self.N)
        return (fr_s[rows, np.argmax(cw >= lo, 1)], fr_s[rows, np.argmax(cw >= hi, 1)])


# ---------------------------------------------------------------- policies
class Periodic:
    """A course every q cycles from cycle `start` (q = 1 front-loaded fixed, q = 4 paced;
    start > 0 delays the first course, script 167)."""
    uses_filter = False

    def __init__(self, q, start=0):
        self.q, self.start = q, start

    def reset(self, scen, rng):
        pass

    def act(self, j, y):
        return np.full(y.shape, j >= self.start and (j - self.start) % self.q == 0)

    def advance(self, j, treat):
        pass


class RawThreshold:
    """On/off rule on raw readings m_j = exp(y_j - y_0): start when m >= a_on, stop when
    m <= a_off (Zhang et al. 2017 is a_on = 1.0, a_off = 0.5)."""
    uses_filter = False

    def __init__(self, a_on, a_off):
        self.a_on, self.a_off = a_on, a_off

    def reset(self, scen, rng):
        self.y0, self.on, self.last = None, None, None

    def act(self, j, y):
        if j == 0:
            self.y0, self.on = y.copy(), np.zeros(y.shape, bool)
        if self.last is not None:               # NaN = known missing scan: reuse the last reading
            y = np.where(np.isnan(y), self.last, y)
        self.last = y
        m = np.exp(y - self.y0)
        self.on = np.where(self.on, m > self.a_off, m >= self.a_on)
        return self.on.copy()

    def advance(self, j, treat):
        pass


class FilteredThreshold:
    """Hybrid: on/off rule on filtered burden / filtered baseline. The start threshold is
    a_R instead of a_on while P(f_R >= b) >= p_trig (resistance mode). a_R = a_on turns
    the resistance mode off (ablation); single=True drops the on/off memory and treats
    exactly when m >= a (no-hysteresis ablation)."""
    uses_filter = True

    def __init__(self, a_on, a_off, a_R, b, single=False, p_trig=0.8, n_particles=500, ranges=RANGES):
        self.a_on, self.a_off, self.a_R, self.b, self.single = a_on, a_off, a_R, b, single
        self.p_trig, self.n_particles, self.ranges = p_trig, n_particles, ranges

    def reset(self, scen, rng):
        self.pf = ParticleFilter(scen["sigma"], self.n_particles, rng, ranges=self.ranges)
        self.ref, self.on = None, None

    def act(self, j, y):
        self.pf.update(j, y)
        Bhat = self.pf.burden_mean()
        if j == 0:
            self.ref, self.on = Bhat.copy(), np.zeros(y.shape, bool)
        m = Bhat / self.ref
        a = np.where(self.pf.prob_fr_ge(self.b) >= self.p_trig, self.a_R, self.a_on)
        self.on = (m >= a) if self.single else np.where(self.on, m > self.a_off, m >= a)
        return self.on.copy()

    def advance(self, j, treat):
        self.pf.propagate(j, treat)


class RiskContainment:
    """Belief containment: treat when the filter's probability of reaching the
    progression level within `horizon` untreated cycles is >= p_star. Each particle uses
    its own 1.2 x B0 level (B0 is a filtered parameter). Motivated by the known-parameter
    oracle, which mostly doses near 1.1-1.19 x B0 (script 153 follow-up; training seeds)."""
    uses_filter = True

    def __init__(self, p_star, horizon=1, guard=0.0, wide=False, n_particles=500, ranges=None, front=0):
        """guard > 0: weak-drug guard (script 156). Treat every cycle while the filter gives
        P(m < guard) >= 0.5, since waiting near the line fails when one course cannot push
        burden back (script 155, M3). wide=True widens the filter prior on m to [3, 20] so
        the filter can represent a weak drug. front > 0: also treat in cycles 0..front-1
        (front-load then contain, script 171; aimed at the lower tail)."""
        self.p_star, self.horizon, self.guard, self.n_particles = p_star, horizon, guard, n_particles
        self.ranges = ranges or (RANGES_WIDE_M if wide else RANGES)
        self.front = front

    def reset(self, scen, rng):
        self.pf = ParticleFilter(scen["sigma"], self.n_particles, rng, ranges=self.ranges)

    def act(self, j, y):
        pf = self.pf
        pf.update(j, y)
        level = PROG_FACTOR * pf.theta["B0"]
        ttp = np.full(pf.S.shape, np.nan)
        S, R = pf.S, pf.R
        for h in range(self.horizon):
            S, R = run_cycle(S, R, pf.theta, 0.0, 0.0, pf.dt, prog_level=level, ttp=ttp)
        risk = (pf.w * ~np.isnan(ttp)).sum(1)
        treat = risk >= self.p_star
        if self.guard > 0:
            treat |= (pf.w * (pf.theta["m"] < self.guard)).sum(1) >= 0.5
        if j < self.front:
            treat[:] = True
        return treat

    def advance(self, j, treat):
        self.pf.propagate(j, treat)


# ---------------------------------------------------------------- episode
def simulate(scen, policy, seed, track_fr=False):
    """Run one noise replicate of every scenario under `policy` (shared course budget).
    Returns capped TTP (days), courses used, and optional filter coverage records."""
    rng = np.random.default_rng(seed)
    N = len(scen["r"])
    # wrong-model suites (script 155) may set: sigma_true (actual noise sd; policies still get
    # "sigma"), noise_df (Student-t readings scaled to the same sd), bias (log offset);
    # all are per-scenario arrays so subset() keeps working
    sd = scen.get("sigma_true", scen["sigma"])[:, None]
    if "noise_df" in scen:
        df = float(scen["noise_df"][0])
        noise = rng.standard_t(df, (N, N_CYCLES)) * np.sqrt((df - 2) / df) * sd
    else:
        noise = rng.standard_normal((N, N_CYCLES)) * sd
    if "bias" in scen:
        noise = noise + scen["bias"][:, None]
    # script 169 truth-only options (per-scenario arrays; policies are not told): mu (S -> R
    # switching under drug, in _deriv), drift (fractional loss of m per treated course),
    # miss_p (reading lost; the last reading is repeated, or NaN if miss_visible is set so the
    # policy knows the scan is missing, script 170), delay (days 0/7/14: the reading shows
    # burden that long before the decision), prog_factor (true endpoint level / B0).
    # miss_p_hi (script 171): state-linked loss, prob miss_p_hi while true burden >= B0, else miss_p
    p = scen
    if "drift" in scen:
        p = dict(scen, m=scen["m"].copy())
    u_miss = None
    if "miss_p" in scen:                    # own stream, so the noise draws above are unchanged
        u_miss = np.random.default_rng([*np.atleast_1d(seed).tolist(), 1]).random((N, N_CYCLES))
    delays = sorted(set(np.asarray(scen.get("delay", [0])).tolist()) - {0})
    S = (1 - scen["f0"]) * scen["B0"]
    R = scen["f0"] * scen["B0"]
    past = {d: S + R for d in delays}
    level = scen.get("prog_factor", PROG_FACTOR) * scen["B0"]
    ttp = np.full(N, np.nan)
    used = np.zeros(N, int)
    hits, total = 0, 0
    policy.reset(scen, rng)
    y_prev = None
    for j in range(N_CYCLES):
        B = S + R
        for d in delays:
            B = np.where(scen["delay"] == d, past[d], B)
        y = np.log(B) + noise[:, j]
        if u_miss is not None:
            if j > 0:
                q = scen["miss_p"]
                if "miss_p_hi" in scen:
                    q = np.where(S + R >= scen["B0"], scen["miss_p_hi"], q)
                y = np.where(u_miss[:, j] < q, np.nan if "miss_visible" in scen else y_prev, y)
            y_prev = y
        treat = policy.act(j, y) & (used < BUDGET)
        if track_fr and policy.uses_filter:
            lo, hi = policy.pf.fr_interval()
            fr_true = R / (S + R)
            live = np.isnan(ttp)
            hits += int(((fr_true >= lo) & (fr_true <= hi) & live).sum())
            total += int(live.sum())
        treat = treat & np.isnan(ttp)           # nothing is given after progression
        used += treat
        snaps = {CYCLE_DAYS - d: None for d in delays} if delays else None
        S, R = run_cycle(S, R, p, treat, j * CYCLE_DAYS, prog_level=level, ttp=ttp, snaps=snaps)
        for d in delays:
            past[d] = sum(snaps[CYCLE_DAYS - d])
        if "drift" in scen:
            p["m"] = np.where(treat, p["m"] * (1.0 - scen["drift"]), p["m"])
        policy.advance(j, treat)
        if not np.isnan(ttp).any():
            break
    capped = np.isnan(ttp)
    ttp[capped] = HORIZON
    out = {"ttp": ttp, "capped": capped, "courses": used}
    if track_fr and policy.uses_filter:
        out["fr_cover"] = (hits, total)
    return out


# ---------------------------------------------------------------- statistics
def cvar_lower(x, alpha=0.10):
    """Mean of the lowest alpha fraction, with fractional weight at the cutoff."""
    x = np.sort(np.asarray(x, float))
    m = alpha * x.size
    k = int(np.floor(m))
    tail = x[:k].sum() + (m - k) * (x[k] if k < x.size else 0.0)
    return tail / m


def paired_bootstrap_cvar(a, b, n_boot=2000, seed=0, alpha=0.10):
    """CI of CVaR(a) - CVaR(b), resampling scenarios (pairs kept together)."""
    rng = np.random.default_rng(seed)
    n = len(a)
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        diffs[i] = cvar_lower(a[idx], alpha) - cvar_lower(b[idx], alpha)
    return float(cvar_lower(a, alpha) - cvar_lower(b, alpha)), np.percentile(diffs, [2.5, 97.5]).tolist()


def evaluate(scen, policy, seeds):
    """Capped TTP per scenario averaged over noise replicates ("ttp"), plus the unaveraged
    replicate x scenario array ("ttp_reps", shape (len(seeds), N)) and max courses used."""
    runs = [simulate(scen, policy, s) for s in seeds]
    reps = np.stack([r["ttp"] for r in runs])
    return {"ttp": reps.mean(0), "ttp_reps": reps,
            "capped_frac": float(np.mean([r["capped"].mean() for r in runs])),
            "courses": np.mean([r["courses"] for r in runs], 0),
            "max_courses": int(max(r["courses"].max() for r in runs))}


# ---------------------------------------------------------------- estimands (audit C2, 2026-10-26)
# Two different lower-tail quantities; name the one you report.
#   cvar_scenario_mean: CVaR_alpha over scenarios of the noise-averaged TTP. This is what
#                       scripts 150-158 report as "cvar10".
#   cvar_joint:         CVaR_alpha over all scenario x noise-replicate outcomes.
def cvar_scenario_mean(reps, alpha=0.10):
    return cvar_lower(np.asarray(reps).mean(0), alpha)


def cvar_joint(reps, alpha=0.10):
    return cvar_lower(np.ravel(reps), alpha)


def hier_bootstrap_diff(a_reps, b_reps, stat, n_boot=2000, seed=0):
    """Paired two-stage bootstrap of stat(a) - stat(b) for (n_rep, N) outcome arrays.
    Stage 1 resamples scenarios; stage 2 resamples noise replicates within each drawn
    scenario. Both arms get the same indices, since replicate k uses the same noise seed in
    every arm. With few replicates stage 2 is crude; report the replicate count."""
    a, b = np.asarray(a_reps, float), np.asarray(b_reps, float)
    n_rep, n = a.shape
    rng = np.random.default_rng(seed)
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        sc = np.broadcast_to(rng.integers(0, n, n), (n_rep, n))
        rp = rng.integers(0, n_rep, (n_rep, n))
        diffs[i] = stat(a[rp, sc]) - stat(b[rp, sc])
    return float(stat(a) - stat(b)), np.percentile(diffs, [2.5, 97.5]).tolist()


def score_spec_reps(args):
    """Worker entry point keeping replicates: (spec, scen, seeds) -> (spec, (n_rep, N) TTP)."""
    spec, scen, seeds = args
    return spec, evaluate(scen, make_policy(spec), seeds)["ttp_reps"]


def make_policy(spec):
    kind = spec[0]
    if kind == "periodic":
        return Periodic(*spec[1:])
    if kind == "raw":
        return RawThreshold(*spec[1:])
    if kind == "filt":
        return FilteredThreshold(*spec[1:])
    if kind == "filtw":                     # script 169 ablation: filtered single threshold, wide prior
        a = spec[1]
        return FilteredThreshold(a, a, a, 0.5, single=True, ranges=RANGES_WIDE_M)
    if kind == "risk":
        return RiskContainment(*spec[1:])
    if kind == "riskf":                     # script 171 hybrid: front-load K courses, then wide containment
        p_star, horizon, k = spec[1:]
        return RiskContainment(p_star, horizon, 0.0, True, front=k)
    raise ValueError(kind)


def score_spec(args):
    """Worker entry point: (spec, scen, seeds) -> (spec, per-scenario mean TTP)."""
    spec, scen, seeds = args
    res = evaluate(scen, make_policy(spec), seeds)
    return spec, res["ttp"]
