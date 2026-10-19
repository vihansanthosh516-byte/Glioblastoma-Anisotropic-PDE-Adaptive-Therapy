"""Small synthetic edge cases from audit chunk 069; no patient data or output writes."""
import pandas as pd
import pytest

from src.patient_stats import bootstrap_ci, cluster_bootstrap
from src.radiation_model import RadiationSchedule, standard_stupp_schedule


def test_empty_bootstrap_rejects_missing_patient_values():
    with pytest.raises(ValueError, match="at least one patient value"):
        bootstrap_ci([], n_boot=10)


def test_cluster_bootstrap_rejects_missing_patient_groups():
    for frame in (
        pd.DataFrame(columns=["patient_id", "value"]),
        pd.DataFrame({"patient_id": [None], "value": [1.0]}),
    ):
        with pytest.raises(ValueError, match="at least one patient group"):
            cluster_bootstrap(frame, lambda sample: sample["value"].mean(), n_boot=10)


def test_total_radiation_dose_is_invariant_to_fractional_time_origin():
    for start in (0.0, 0.25, 0.5, -0.5):
        assert standard_stupp_schedule(start).total_dose() == 60.0
        assert RadiationSchedule(start, start).total_dose() == 2.0
        assert RadiationSchedule(start, start + 4.9).total_dose() == 10.0


def test_tensor_builder_rejects_nonfinite_inputs_and_regularizes_negative_eigenvalues():
    import numpy as np
    from src.load_ucsf_tensor import build_pd_tensor

    for values in ((float("nan"), 1.0, 0.0), (1.0, float("inf"), 0.0), (1.0, 1.0, float("nan"))):
        with pytest.raises(ValueError, match="must be finite"):
            build_pd_tensor(*values)
    tensor = build_pd_tensor(-2.0, 0.0, 0.6)
    assert np.allclose(tensor, tensor.T)
    assert np.linalg.eigvalsh(tensor).min() >= 1e-6 - 1e-12
