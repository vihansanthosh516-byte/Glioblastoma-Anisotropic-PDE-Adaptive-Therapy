import importlib.util
from pathlib import Path

import numpy as np

from src.treatment_aware_pde import treatment_aware_ode_model


_spec = importlib.util.spec_from_file_location("joint", Path("src/51_joint_inverse_estimation.py"))
_joint = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_joint)


def test_joint_estimation_uses_all_intervals():
    rho, diffusion = 0.025, 0.012
    days = [0.0, 10.0, 20.0, 30.0]
    volumes = [1000.0]
    for start, end in zip(days, days[1:]):
        volumes.append(treatment_aware_ode_model(rho, diffusion, volumes[-1], end - start, start_day=start))
    result = _joint.estimate_joint_parameters(volumes, days, regularization=0.0)
    assert abs(result["rho"] - rho) < 0.003
    assert abs(result["D"] - diffusion) < 0.004
    assert result["n_intervals"] == 3
