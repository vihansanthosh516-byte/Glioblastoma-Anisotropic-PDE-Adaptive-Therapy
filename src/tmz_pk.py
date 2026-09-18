"""One-compartment temozolomide pharmacokinetics.

Concentrations are normalized to the amount of a single bolus. Time is in days.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp, log
from typing import Iterable, Sequence


@dataclass(frozen=True)
class TMZPK:
    half_life_days: float = 1.8 / 24.0

    @property
    def elimination_rate(self) -> float:
        if self.half_life_days <= 0:
            raise ValueError("half_life_days must be positive")
        return log(2.0) / self.half_life_days

    def concentration(self, time_days: float, bolus_days: Iterable[float] = ()) -> float:
        """Return normalized concentration after boluses at the supplied times."""
        if time_days < 0:
            return 0.0
        k_el = self.elimination_rate
        return sum(exp(-k_el * (time_days - day)) for day in bolus_days if day <= time_days)

    def profile(self, times_days: Sequence[float], bolus_days: Iterable[float] = ()) -> list[float]:
        boluses = tuple(bolus_days)
        return [self.concentration(time, boluses) for time in times_days]


def tmz_concentration(time_days: float, bolus_days: Iterable[float] = (), half_life_days: float = 1.8 / 24.0) -> float:
    return TMZPK(half_life_days).concentration(time_days, bolus_days)
