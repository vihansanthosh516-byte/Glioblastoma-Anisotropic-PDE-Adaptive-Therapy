"""Fractionated external-beam radiation schedule utilities."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class RadiationSchedule:
    start_day: float
    end_day: float
    dose_per_fraction_gy: float = 2.0
    fractions_per_week: int = 5
    oer: float = 1.0
    let_correction: float = 1.0

    def __post_init__(self) -> None:
        if self.end_day < self.start_day:
            raise ValueError("end_day must be >= start_day")
        if self.dose_per_fraction_gy < 0 or self.fractions_per_week < 1:
            raise ValueError("invalid radiation dose or fraction count")
        if self.oer <= 0 or self.let_correction <= 0:
            raise ValueError("OER and LET correction must be positive")

    @property
    def corrected_dose_per_fraction(self) -> float:
        return self.dose_per_fraction_gy * self.oer * self.let_correction

    def dose_on_day(self, day: float) -> float:
        """Return Gy delivered on a day; weekdays are anchored at start_day."""
        if day < self.start_day or day > self.end_day:
            return 0.0
        treatment_day = int(day - self.start_day + 1e-9)
        return self.corrected_dose_per_fraction if treatment_day % 7 < self.fractions_per_week else 0.0

    def total_dose(self) -> float:
        first = int(self.start_day)
        last = int(self.end_day)
        return sum(self.dose_on_day(day) for day in range(first, last + 1))

    def dose_rate(self, day: float, fraction_duration_days: float = 1.0) -> float:
        if fraction_duration_days <= 0:
            raise ValueError("fraction_duration_days must be positive")
        return self.dose_on_day(day) / fraction_duration_days


def standard_stupp_schedule(start_day: float, end_day: float | None = None) -> RadiationSchedule:
    """Create the conventional 60 Gy / 30 weekday fraction schedule."""
    if end_day is None:
        end_day = start_day + 41
    return RadiationSchedule(start_day=start_day, end_day=end_day)
