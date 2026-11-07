"""Tabular Q-learning baseline for Track C (SIMULATION ONLY).

Adapted from Yauney & Shah (2018): learn dosing choices from a discretized recent
tumor-burden history, add an action-derived course penalty, and use an episode outcome
as the terminal reward. This is not a replication of their Double-DQN/TMZ-PCV study.
Reference: Yauney & Shah, PMLR 85 (2018), 161-226,
https://proceedings.mlr.press/v85/yauney18a.html.

The policy sees only noisy scan readings, elapsed cycle, and courses actually given.
Training uses only caller-supplied selection scenarios. Simulated TTP is used as the
selection-cohort reward; evaluation patients are never accepted by the trainer.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

import trackc_sim as tc


# Bins for log(burden / first observed burden). Edges correspond to burden ratios
# 0.5, 0.8, 1.0, 1.25 and 1.5. They are fixed method choices, not fitted on evaluation data.
BURDEN_LOG_EDGES = (-0.6931471805599453, -0.2231435513142097, 0.0,
                    0.22314355131420976, 0.4054651081081644)
TREND_LOG_THRESHOLD = 0.10
N_BURDEN_BINS = len(BURDEN_LOG_EDGES) + 1
N_TREND_BINS = 3
N_ACTIONS = 2                         # 0 = no course, 1 = give one course
DEFAULT_SCAN_EVERY = 3                # 3 x 28 days = 84 days
DEFAULT_TAU_DAYS = 730.0


def q_shape(course_cap: int = tc.BUDGET) -> tuple[int, ...]:
    """Shape of the Q table: cycle x burden x trend x courses-used x action."""
    cap = _check_course_cap(course_cap)
    return (tc.N_CYCLES, N_BURDEN_BINS, N_TREND_BINS, cap + 1, N_ACTIONS)


def _check_course_cap(course_cap: int) -> int:
    if isinstance(course_cap, bool) or int(course_cap) != course_cap:
        raise ValueError("course_cap must be an integer")
    cap = int(course_cap)
    if not 0 <= cap <= tc.BUDGET:
        raise ValueError(f"course_cap must be between 0 and {tc.BUDGET}")
    return cap


def _check_scan_every(scan_every: int) -> int:
    if isinstance(scan_every, bool) or int(scan_every) != scan_every or int(scan_every) < 1:
        raise ValueError("scan_every must be a positive integer number of 28-day cycles")
    return int(scan_every)


def bin_state(cycle: int, log_relative_burden: float, trend_delta: float,
              courses_used: int, *, course_cap: int = tc.BUDGET) -> tuple[int, int, int, int]:
    """Map scan history to (cycle, burden bin, trend bin, courses used).

    log_relative_burden is the latest noisy log burden minus the first noisy
    reading. trend_delta is the change between the latest two available scans;
    negative, near-zero, and positive changes map to falling, stable, and rising.
    """
    cap = _check_course_cap(course_cap)
    cycle_i = int(cycle)
    if isinstance(cycle, bool) or cycle_i != cycle or not 0 <= cycle_i < tc.N_CYCLES:
        raise ValueError(f"cycle must be between 0 and {tc.N_CYCLES - 1}")
    rel = float(log_relative_burden)
    change = float(trend_delta)
    if not np.isfinite(rel) or not np.isfinite(change):
        raise ValueError("burden and trend values must be finite")
    used = int(courses_used)
    if isinstance(courses_used, bool) or used != courses_used or used < 0:
        raise ValueError("courses_used must be a non-negative integer")
    burden_bin = int(np.searchsorted(BURDEN_LOG_EDGES, rel, side="right"))
    trend_bin = 0 if change < -TREND_LOG_THRESHOLD else (2 if change > TREND_LOG_THRESHOLD else 1)
    return cycle_i, burden_bin, trend_bin, min(used, cap)


def bellman_update(q_values: np.ndarray, state: tuple[int, ...], action: int,
                   reward: float, next_state: tuple[int, ...] | None, *,
                   learning_rate: float = 0.1, discount: float = 1.0,
                   terminal: bool = False) -> float:
    """Apply one tabular Q-learning Bellman update in place; return the new value."""
    q = np.asarray(q_values)
    if q.ndim != 5 or q.shape[-1] != N_ACTIONS:
        raise ValueError("q_values must have cycle, burden, trend, course, action axes")
    if q.dtype.kind not in "fc":
        raise ValueError("q_values must use a floating-point dtype")
    if not 0.0 < learning_rate <= 1.0:
        raise ValueError("learning_rate must be in (0, 1]")
    if not 0.0 <= discount <= 1.0:
        raise ValueError("discount must be in [0, 1]")
    if action not in (0, 1):
        raise ValueError("action must be 0 or 1")
    if not np.isfinite(reward):
        raise ValueError("reward must be finite")
    old = float(q[state + (int(action),)])
    target = float(reward)
    if not terminal:
        if next_state is None:
            raise ValueError("next_state is required for a non-terminal update")
        target += discount * float(np.max(q[next_state]))
    q[state + (int(action),)] = old + learning_rate * (target - old)
    return float(q[state + (int(action),)])


class QTablePolicy:
    """Deterministic Track C policy implementing reset/act/advance.

    Give trackc_sim.simulate scenarios with scan_every=policy.scan_every and
    budget=policy.course_cap. evaluate_q_policy applies both automatically.
    The policy's reset deliberately ignores the scenario dictionary; its only
    patient-specific inputs are the noisy y readings and delivered actions.
    """
    uses_filter = False

    def __init__(self, q_values: np.ndarray, *, course_cap: int = tc.BUDGET,
                 scan_every: int = DEFAULT_SCAN_EVERY):
        self.course_cap = _check_course_cap(course_cap)
        self.scan_every = _check_scan_every(scan_every)
        expected = q_shape(self.course_cap)
        q = np.asarray(q_values, dtype=float)
        if q.shape != expected or not np.all(np.isfinite(q)):
            raise ValueError(f"q_values must be finite and have shape {expected}")
        self.q_values = q.copy()
        self._record = False
        self._epsilon = 0.0
        self._initialize_state()

    def _initialize_state(self):
        self._rng = None
        self._baseline = None
        self._last_scan = None
        self._relative = None
        self._trend_sign = None
        self._used = None
        self._pending_states = None
        self._pending_actions = None
        self._state_history = []
        self._action_history = []
        self._delivered_history = []

    def reset(self, scen: Mapping, rng: np.random.Generator):
        """Reset per-patient memory without reading any scenario parameter."""
        del scen
        self._initialize_state()
        self._rng = rng

    def _state_rows(self, cycle: int) -> np.ndarray:
        return np.asarray([
            bin_state(cycle, self._relative[i], self._trend_sign[i] * 2 * TREND_LOG_THRESHOLD,
                      self._used[i], course_cap=self.course_cap)
            for i in range(len(self._used))
        ], dtype=np.int16)

    def act(self, j: int, y: np.ndarray) -> np.ndarray:
        """Choose 0/1 treatment using the latest visible noisy scan; carry it over gaps."""
        obs = np.asarray(y, dtype=float)
        if obs.ndim != 1:
            raise ValueError("y must be a one-dimensional patient vector")
        if obs.size == 0:
            raise ValueError("y must contain at least one patient")
        cycle = int(j)
        if cycle != j or not 0 <= cycle < tc.N_CYCLES:
            raise ValueError(f"cycle must be between 0 and {tc.N_CYCLES - 1}")
        if self._baseline is None:
            n = obs.size
            self._baseline = np.full(n, np.nan)
            self._last_scan = np.full(n, np.nan)
            self._relative = np.zeros(n, dtype=float)
            self._trend_sign = np.zeros(n, dtype=np.int8)
            self._used = np.zeros(n, dtype=np.int16)
        elif obs.size != self._baseline.size:
            raise ValueError("patient batch size changed within an episode")

        seen = np.isfinite(obs)
        old_scan = self._last_scan.copy()
        initial = seen & ~np.isfinite(self._baseline)
        self._baseline[initial] = obs[initial]
        self._last_scan[initial] = obs[initial]
        has_pair = seen & np.isfinite(old_scan) & ~initial
        delta = obs - old_scan
        self._trend_sign[has_pair & (delta < -TREND_LOG_THRESHOLD)] = -1
        self._trend_sign[has_pair & (delta > TREND_LOG_THRESHOLD)] = 1
        self._trend_sign[has_pair & (np.abs(delta) <= TREND_LOG_THRESHOLD)] = 0
        comparable = seen & np.isfinite(self._baseline)
        self._relative[comparable] = obs[comparable] - self._baseline[comparable]
        self._last_scan[seen] = obs[seen]

        states = self._state_rows(cycle)
        q_rows = np.stack([self.q_values[tuple(row)] for row in states])
        actions = np.argmax(q_rows, axis=1).astype(np.int8)  # ties favor no treatment
        if self._epsilon > 0.0:
            if self._rng is None:
                raise RuntimeError("training policy has not been reset with an RNG")
            explore = self._rng.random(obs.size) < self._epsilon
            actions[explore] = self._rng.integers(0, N_ACTIONS, int(explore.sum()), dtype=np.int8)
        actions[self._used >= self.course_cap] = 0
        self._pending_states = states.copy()
        self._pending_actions = actions.copy()
        return actions.astype(bool)

    def advance(self, j: int, treat: np.ndarray):
        """Track only courses actually delivered by the simulator's budget mask."""
        del j
        actual = np.asarray(treat, dtype=bool)
        if self._used is None or actual.shape != self._used.shape:
            raise ValueError("advance must follow act with one treatment flag per patient")
        if self._record:
            self._state_history.append(self._pending_states.copy())
            self._action_history.append(self._pending_actions.copy())
            self._delivered_history.append(actual.astype(np.int8))
        self._used = np.minimum(self.course_cap, self._used + actual.astype(np.int16))


class _ExploringPolicy(QTablePolicy):
    def __init__(self, q_values: np.ndarray, *, course_cap: int, scan_every: int, epsilon: float):
        super().__init__(q_values, course_cap=course_cap, scan_every=scan_every)
        self._epsilon = float(epsilon)
        self._record = True


def _seed_components(seed) -> list[int]:
    """Normalize an integer or structured seed to non-negative integer components."""
    raw = np.atleast_1d(seed).tolist()
    components = []
    for value in raw:
        if isinstance(value, (bool, np.bool_)):
            raise ValueError("seed components must be non-negative integers")
        try:
            component = int(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("seed components must be non-negative integers") from exc
        if component != value or component < 0:
            raise ValueError("seed components must be non-negative integers")
        components.append(component)
    if not components:
        raise ValueError("seed must contain at least one component")
    return components


def train_q_learning(selection_scenarios: Mapping, *, episodes: int = 50,
                     course_cap: int = tc.BUDGET, course_penalty: float = 10.0,
                     scan_every: int = DEFAULT_SCAN_EVERY, tau_days: float = DEFAULT_TAU_DAYS,
                     learning_rate: float = 0.1, discount: float = 1.0,
                     epsilon: float = 0.25, epsilon_decay: float = 0.98,
                     min_epsilon: float = 0.02, seed: int = 0,
                     endpoint: str = "rano") -> QTablePolicy:
    """Fit a shared tabular policy on selection scenarios only.

    Each cycle has a binary action (one course or none). Observations are supplied by
    the simulator every scan_every cycles, defaulting to 84 days. The episode return
    is restricted TTP in days minus course_penalty days per delivered course. A
    terminal Q reward uses simulated TTP; intermediate rewards are the action-derived
    course penalty. The undiscounted course penalty covers every delivered course through
    the simulator horizon, including courses after ``tau_days``. No internal patient/world
    draws or evaluation scenarios are used.
    """
    cap = _check_course_cap(course_cap)
    scans = _check_scan_every(scan_every)
    if isinstance(episodes, bool) or int(episodes) != episodes or int(episodes) < 1:
        raise ValueError("episodes must be a positive integer")
    if not np.isfinite(course_penalty) or course_penalty < 0:
        raise ValueError("course_penalty must be finite and non-negative")
    if not np.isfinite(tau_days) or tau_days <= 0:
        raise ValueError("tau_days must be finite and positive")
    if not 0.0 < learning_rate <= 1.0 or not 0.0 <= discount <= 1.0:
        raise ValueError("learning_rate must be in (0, 1] and discount in [0, 1]")
    if float(discount) != 1.0:
        raise ValueError("this baseline uses an undiscounted RMST-minus-course-cost objective; discount must be 1")
    if not 0.0 <= min_epsilon <= epsilon <= 1.0 or not 0.0 <= epsilon_decay <= 1.0:
        raise ValueError("epsilon values must satisfy 0 <= min_epsilon <= epsilon <= 1 and 0 <= epsilon_decay <= 1")
    seed_components = _seed_components(seed)

    q = np.zeros(q_shape(cap), dtype=float)
    training_scenarios = dict(selection_scenarios)
    training_scenarios["scan_every"] = scans
    current_epsilon = float(epsilon)
    for episode in range(int(episodes)):
        learner = _ExploringPolicy(q, course_cap=cap, scan_every=scans, epsilon=current_epsilon)
        sim_seed = [*seed_components, episode]
        outcome = tc.simulate(training_scenarios, learner, sim_seed, endpoint=endpoint, budget=cap)
        if not learner._state_history:
            raise RuntimeError("simulator produced no policy transitions")
        states = np.stack(learner._state_history)           # (cycle, patient, state-dimension)
        actions = np.stack(learner._action_history)
        delivered = np.stack(learner._delivered_history)
        n_steps, n_patients = actions.shape
        for patient in range(n_patients):
            ttp = float(outcome["ttp"][patient])
            terminal_step = int(np.ceil(ttp / tc.CYCLE_DAYS - 1e-12)) - 1
            terminal_step = int(np.clip(terminal_step, 0, n_steps - 1))
            for step in range(terminal_step, -1, -1):
                state = tuple(int(v) for v in states[step, patient])
                action = int(actions[step, patient])
                reward = -float(course_penalty) * int(delivered[step, patient])
                terminal = step == terminal_step
                if terminal:
                    reward += min(ttp, float(tau_days))
                    next_state = None
                else:
                    next_state = tuple(int(v) for v in states[step + 1, patient])
                bellman_update(q, state, action, reward, next_state,
                               learning_rate=learning_rate, discount=discount, terminal=terminal)
        current_epsilon = max(float(min_epsilon), current_epsilon * float(epsilon_decay))

    return QTablePolicy(q, course_cap=cap, scan_every=scans)


def evaluate_q_policy(evaluation_scenarios: Mapping, policy: QTablePolicy, seed,
                      *, endpoint: str = "rano") -> dict:
    """Evaluate a frozen deterministic policy with 84-day scans and its fitted course cap."""
    if not isinstance(policy, QTablePolicy):
        raise TypeError("policy must be a fitted QTablePolicy")
    eval_scenarios = dict(evaluation_scenarios)
    eval_scenarios["scan_every"] = policy.scan_every
    return tc.simulate(eval_scenarios, policy, seed, endpoint=endpoint, budget=policy.course_cap)


def train_q_learning_shared_endpoints(
    selection_scenarios: Mapping, *, episodes: int = 50,
    course_cap: int = tc.BUDGET, course_penalty: float = 10.0,
    scan_every: int = DEFAULT_SCAN_EVERY, tau_days: float = DEFAULT_TAU_DAYS,
    learning_rate: float = 0.1, discount: float = 1.0,
    epsilon: float = 0.25, epsilon_decay: float = 0.98,
    min_epsilon: float = 0.02, seed: int = 0,
    stop_endpoint: str = "either",
) -> QTablePolicy:
    """Fit Q-learning to both endpoints on the simulator's shared treatment path.

    Training episodes use the requested shared treatment-stop rule. For latent-stop
    modes, the terminal reward is assigned at the first endpoint event. For
    ``observed_two_scan``, hidden endpoint times do not truncate learning: treatment
    follows the observed two-scan stop rule, endpoint utility is scored through
    ``tau_days``, and the terminal reward is assigned at the observed stop or the
    restriction horizon. Selection scenarios only are accepted by this function;
    callers must keep evaluation scenarios separate.

    This is a simulation baseline, not a replication of Yauney & Shah (2018),
    a patient-specific model, or a clinical treatment policy.
    """
    cap = _check_course_cap(course_cap)
    scans = _check_scan_every(scan_every)
    if stop_endpoint not in ("baseline", "rano", "either", "observed_two_scan"):
        raise ValueError("stop_endpoint must be baseline, rano, either, or observed_two_scan")
    if isinstance(episodes, bool) or int(episodes) != episodes or int(episodes) < 1:
        raise ValueError("episodes must be a positive integer")
    if not np.isfinite(course_penalty) or course_penalty < 0:
        raise ValueError("course_penalty must be finite and non-negative")
    if not np.isfinite(tau_days) or tau_days <= 0:
        raise ValueError("tau_days must be finite and positive")
    if not 0.0 < learning_rate <= 1.0 or not 0.0 <= discount <= 1.0:
        raise ValueError("learning_rate must be in (0, 1] and discount in [0, 1]")
    if float(discount) != 1.0:
        raise ValueError("this baseline uses an undiscounted shared-endpoint objective; discount must be 1")
    if not 0.0 <= min_epsilon <= epsilon <= 1.0 or not 0.0 <= epsilon_decay <= 1.0:
        raise ValueError("epsilon values must satisfy 0 <= min_epsilon <= epsilon <= 1 and 0 <= epsilon_decay <= 1")
    seed_components = _seed_components(seed)

    q = np.zeros(q_shape(cap), dtype=float)
    training_scenarios = dict(selection_scenarios)
    training_scenarios["scan_every"] = scans
    current_epsilon = float(epsilon)
    for episode in range(int(episodes)):
        learner = _ExploringPolicy(q, course_cap=cap, scan_every=scans, epsilon=current_epsilon)
        sim_seed = [*seed_components, episode]
        outcome = tc.simulate_two_endpoints(
            training_scenarios, learner, sim_seed, stop_endpoint=stop_endpoint, budget=cap
        )
        if not learner._state_history:
            raise RuntimeError("simulator produced no policy transitions")
        required = ("ttp_baseline", "ttp_rano")
        if not all(key in outcome for key in required):
            raise ValueError("shared-endpoint simulator must return both progression times")
        ttp_baseline = np.asarray(outcome["ttp_baseline"], dtype=float)
        ttp_rano = np.asarray(outcome["ttp_rano"], dtype=float)
        if (ttp_baseline.ndim != 1 or ttp_rano.shape != ttp_baseline.shape
                or not np.isfinite(ttp_baseline).all() or not np.isfinite(ttp_rano).all()
                or np.any(ttp_baseline < 0) or np.any(ttp_rano < 0)):
            raise ValueError("both progression-time arrays must be matching, finite, and non-negative")
        states = np.stack(learner._state_history)
        actions = np.stack(learner._action_history)
        delivered = np.stack(learner._delivered_history)
        n_steps, n_patients = actions.shape
        if n_patients != ttp_baseline.size:
            raise ValueError("simulator outcomes and policy transitions have different patient counts")
        if "courses" in outcome:
            courses = np.asarray(outcome["courses"])
            if courses.shape != (n_patients,) or not np.isfinite(courses).all():
                raise ValueError("shared-endpoint course counts must match the training patients")

        effective_ttp = np.minimum(np.minimum(ttp_baseline, ttp_rano), float(tau_days))
        if stop_endpoint == "observed_two_scan":
            if "observed_stop_day" not in outcome:
                raise ValueError("observed-stop training requires observed_stop_day outcomes")
            observed_stop_day = np.asarray(outcome["observed_stop_day"], dtype=float)
            if observed_stop_day.shape != (n_patients,) or np.isinf(observed_stop_day).any():
                raise ValueError("observed_stop_day must be a matching vector of finite values or NaN")
            if np.any(observed_stop_day[np.isfinite(observed_stop_day)] < 0):
                raise ValueError("observed_stop_day values must be non-negative")
            last_tau_step = int(np.ceil(float(tau_days) / tc.CYCLE_DAYS - 1e-12)) - 1
            last_tau_step = int(np.clip(last_tau_step, 0, n_steps - 1))
            observed_steps = np.full(n_patients, last_tau_step, dtype=int)
            has_observed_stop = np.isfinite(observed_stop_day)
            observed_steps[has_observed_stop] = np.rint(
                observed_stop_day[has_observed_stop] / tc.CYCLE_DAYS
            ).astype(int)
            terminal_steps = np.clip(np.minimum(observed_steps, last_tau_step), 0, n_steps - 1)
        else:
            terminal_steps = np.ceil(effective_ttp / tc.CYCLE_DAYS - 1e-12).astype(int) - 1
            terminal_steps = np.clip(terminal_steps, 0, n_steps - 1)
        for patient in range(n_patients):
            terminal_step = int(terminal_steps[patient])
            for step in range(terminal_step, -1, -1):
                state = tuple(int(v) for v in states[step, patient])
                action = int(actions[step, patient])
                reward = -float(course_penalty) * int(delivered[step, patient])
                terminal = step == terminal_step
                if terminal:
                    reward += float(effective_ttp[patient])
                    next_state = None
                else:
                    next_state = tuple(int(v) for v in states[step + 1, patient])
                bellman_update(q, state, action, reward, next_state,
                               learning_rate=learning_rate, discount=discount, terminal=terminal)
        current_epsilon = max(float(min_epsilon), current_epsilon * float(epsilon_decay))

    return QTablePolicy(q, course_cap=cap, scan_every=scans)


def evaluate_q_policy_shared_endpoints(
    evaluation_scenarios: Mapping,
    policy: QTablePolicy,
    seed,
    *,
    stop_endpoint: str = "either",
):
    """Evaluate a frozen Q policy on one shared path for both progression endpoints."""
    if not isinstance(policy, QTablePolicy):
        raise TypeError("policy must be a fitted QTablePolicy")
    eval_scenarios = dict(evaluation_scenarios)
    eval_scenarios["scan_every"] = policy.scan_every
    return tc.simulate_two_endpoints(
        eval_scenarios, policy, seed, stop_endpoint=stop_endpoint, budget=policy.course_cap
    )

