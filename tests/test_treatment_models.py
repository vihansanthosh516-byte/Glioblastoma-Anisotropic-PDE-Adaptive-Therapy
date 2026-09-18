from src.radiation_model import RadiationSchedule
from src.tmz_pk import TMZPK
from src.treatment_aware_pde import TreatmentSchedule, treatment_aware_ode_model


def test_tmz_bolus_decays_with_half_life():
    pk = TMZPK()
    assert abs(pk.concentration(pk.half_life_days, [0.0]) - 0.5) < 1e-12


def test_radiation_has_thirty_fractions_and_sixty_gy():
    schedule = RadiationSchedule(0.0, 41.0)
    assert sum(schedule.dose_on_day(day) > 0 for day in range(42)) == 30
    assert schedule.total_dose() == 60.0


def test_treatment_reduces_predicted_volume():
    untreated = treatment_aware_ode_model(0.02, 0.01, 1000.0, 30.0)
    treated = treatment_aware_ode_model(
        0.02,
        0.01,
        1000.0,
        30.0,
        TreatmentSchedule(
            tmz_bolus_days=tuple(float(day) for day in range(30)),
            radiation=RadiationSchedule(0.0, 29.0),
        ),
    )
    assert treated < untreated
