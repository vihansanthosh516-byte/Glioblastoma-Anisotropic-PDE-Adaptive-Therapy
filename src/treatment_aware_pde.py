"""Treatment-aware Fisher--Kolmogorov volume surrogate.

The full spatial PDE can use the same sink terms pointwise. This module keeps
inverse estimation fast while preserving the treatment semantics:

  dV/dt = rho V (1 - V/K) + c_diff D V^(1/3) - alpha C(t)V - beta R(t)V
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np

try:
    from .radiation_model import RadiationSchedule
    from .tmz_pk import TMZPK
except ImportError:  # import-by-path execution
    import importlib.util
    _root = __import__("pathlib").Path(__file__).parent
    def _load(name: str):
        spec = importlib.util.spec_from_file_location(name, _root / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        import sys
        sys.modules[name] = module
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module
    RadiationSchedule = _load("radiation_model").RadiationSchedule
    TMZPK = _load("tmz_pk").TMZPK


@dataclass(frozen=True)
class TreatmentSchedule:
    tmz_bolus_days: tuple[float, ...] = ()
    radiation: RadiationSchedule | None = None

    def tmz_concentration(self, day: float, pk: TMZPK | None = None) -> float:
        return (pk or TMZPK()).concentration(day, self.tmz_bolus_days)

    def radiation_dose_rate(self, day: float) -> float:
        return self.radiation.dose_rate(day) if self.radiation else 0.0

    def active_days(self, start_day: float, end_day: float) -> float:
        """Approximate treatment-active duration for interval masking."""
        if end_day <= start_day:
            return 0.0
        points = np.linspace(start_day, end_day, max(2, int(np.ceil(end_day - start_day)) + 1))
        active = [self.tmz_concentration(float(day)) > 1e-8 or self.radiation_dose_rate(float(day)) > 0 for day in points]
        return float(np.trapz(np.asarray(active, dtype=float), points))


def treatment_aware_ode_model(
    rho: float,
    D: float,
    V0: float,
    delta_t: float,
    schedule: TreatmentSchedule | None = None,
    alpha: float = 0.08,
    beta: float = 0.03,
    K: float = 1.0e6,
    start_day: float = 0.0,
    dt: float = 0.1,
) -> float:
    if V0 <= 0:
        return 0.0
    if delta_t < 0 or dt <= 0:
        raise ValueError("delta_t must be non-negative and dt must be positive")
    schedule = schedule or TreatmentSchedule()
    n_steps = max(1, int(np.ceil(delta_t / dt))) if delta_t else 0
    step = delta_t / n_steps if n_steps else 0.0
    volume = float(V0)
    c_diff = (36.0 * np.pi) ** (1.0 / 3.0)
    for index in range(n_steps):
        day = start_day + index * step
        C = schedule.tmz_concentration(day)
        R = schedule.radiation_dose_rate(day)
        derivative = (
            rho * volume * (1.0 - volume / K)
            + c_diff * D * volume ** (1.0 / 3.0)
            - alpha * C * volume
            - beta * R * volume
        )
        volume = max(0.0, volume + step * derivative)
    return volume
