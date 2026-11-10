import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

import trackc_sim as tc
from trackc_integrated_sampler import (
    B0_CENTER_HALF_WIDTH,
    COHORT_ROLES,
    PATIENT_DIMENSION_ORDER,
    SIGMA_LEVELS,
    WORLD_AXIS_BOUNDS,
    WORLD_MU_NU_LEVELS,
    WORLD_Q_SIGNAL_LEVELS,
    assign_balanced_sigma,
    draw_patient_cohort,
    draw_world_parameters,
)


def test_world_draw_is_deterministic_and_uses_independent_world_support():
    first = draw_world_parameters("single", 7001)
    again = draw_world_parameters("single", 7001)
    other = draw_world_parameters("single", 7002)
    assert first == again
    assert first != other
    assert first["mu"] in WORLD_MU_NU_LEVELS
    assert first["nu"] in WORLD_MU_NU_LEVELS
    assert first["q_signal"] in WORLD_Q_SIGNAL_LEVELS
    assert "c" not in first


def test_mixture_world_keeps_response_fraction_distinct_from_signal_quality():
    world = draw_world_parameters("mixture", 7011)
    assert "pi_response" in world
    assert 0.1 <= world["pi_response"] <= 0.6
    assert "q" not in world
    assert world["q_signal"] in WORLD_Q_SIGNAL_LEVELS


def test_single_patient_draw_uses_declared_axis_order_and_patient_level_cost():
    seed, n = 7021, 12
    world = draw_world_parameters("single", 7020)
    actual = draw_patient_cohort(world, n, seed)

    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    unit = qmc.LatinHypercube(d=len(PATIENT_DIMENSION_ORDER), scramble=True, seed=rng).random(n)
    r_bounds = (tc.RANGES["r"][0] * np.exp(world["log_r_shift"]), tc.RANGES["r"][1] * np.exp(world["log_r_shift"]))
    m_scale = np.exp(world["log_m_scale"])
    m_bounds = (tc.RANGES["m"][0] * m_scale, tc.RANGES["m"][1] * m_scale)
    f0_scale = np.exp(world["log_f0_shift"])
    f0_bounds = (tc.RANGES["f0"][0] * f0_scale, min(tc.RANGES["f0"][1] * f0_scale, 0.9))
    b0_bounds = (world["B0_center"] - 0.1, world["B0_center"] + 0.1)

    expected = {
        "r": np.exp(np.log(r_bounds[0]) + unit[:, 0] * (np.log(r_bounds[1]) - np.log(r_bounds[0]))),
        "m": np.exp(np.log(m_bounds[0]) + unit[:, 1] * (np.log(m_bounds[1]) - np.log(m_bounds[0]))),
        "f0": np.exp(np.log(f0_bounds[0]) + unit[:, 2] * (np.log(f0_bounds[1]) - np.log(f0_bounds[0]))),
        "B0": b0_bounds[0] + unit[:, 3] * (b0_bounds[1] - b0_bounds[0]),
        "c": tc.RANGES["c"][0] + unit[:, 4] * (tc.RANGES["c"][1] - tc.RANGES["c"][0]),
        "eps": tc.RANGES["eps"][0] + unit[:, 5] * (tc.RANGES["eps"][1] - tc.RANGES["eps"][0]),
    }
    for name, values in expected.items():
        np.testing.assert_array_equal(actual[name], values)
    assert len(np.unique(actual["c"])) == n
    assert tuple(PATIENT_DIMENSION_ORDER) == ("r", "m", "f0", "B0", "c", "epsilon")


def test_mixture_patients_are_reproducible_and_within_declared_ranges():
    world = draw_world_parameters("mixture", 7030)
    first = draw_patient_cohort(world, 40, 7031)
    again = draw_patient_cohort(world, 40, 7031)
    for key in first:
        np.testing.assert_array_equal(first[key], again[key])
    assert np.all((first["m"] > 0.0) & (first["m"] <= 20.0))
    assert np.all((first["c"] >= 0.0) & (first["c"] <= 0.25))
    assert np.all((first["eps"] >= 0.0) & (first["eps"] <= 0.2))
    assert "response_class" not in first


def test_sigma_is_balanced_once_across_all_world_cohorts():
    world = draw_world_parameters("single", 7040)
    n = 2
    scenarios = {
        role: draw_patient_cohort(world, n, 7041 + i)
        for i, role in enumerate(COHORT_ROLES)
    }
    patient_ids = {
        role: tuple(f"world-{role}-{i}" for i in range(n))
        for role in COHORT_ROLES
    }
    assigned, record = assign_balanced_sigma(scenarios, patient_ids, world_seed=7040)
    repeated, repeated_record = assign_balanced_sigma(scenarios, patient_ids, world_seed=7040)
    assert record == repeated_record
    assert record["balanced_count_difference_max"] == 0
    assert record["counts"] == {"0.05": 2, "0.1": 2, "0.2": 2}
    for role in COHORT_ROLES:
        np.testing.assert_array_equal(assigned[role]["sigma"], repeated[role]["sigma"])
        assert len(assigned[role]["sigma"]) == n


def test_sigma_assignment_refuses_duplicate_patient_ids_across_cohorts():
    world = draw_world_parameters("single", 7050)
    scenarios = {role: draw_patient_cohort(world, 1, 7051 + i) for i, role in enumerate(COHORT_ROLES)}
    patient_ids = {role: ("duplicate",) for role in COHORT_ROLES}
    try:
        assign_balanced_sigma(scenarios, patient_ids, world_seed=7050)
    except ValueError as exc:
        assert "unique" in str(exc)
    else:
        raise AssertionError("duplicate patient IDs must be refused")


def test_integrated_config_matches_the_sampler_contract():
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "configs" / "trackABC_integrated_v2.json").read_text(encoding="utf-8"))
    sampling = config["tumor_model"]["sampling_algorithm"]
    world_axes = config["tumor_model"]["world_axes"]
    patient_axes = config["tumor_model"]["patient_axes"]
    randomness = config["randomness_specification_draft"]

    assert "independently and uniformly" in sampling["mu_nu_assignment"]
    assert sampling["patient_dimension_order"] == list(PATIENT_DIMENSION_ORDER)
    assert "all sorted patient IDs" in sampling["observation_noise_assignment"]
    assert config["signal_model"]["primary_q_signal_sampling"] == (
        "uniform draw from [0.6, 0.7, 0.8, 0.9] once per world; constant across that world's virtual patients"
    )
    assert config["signal_model"]["primary_q_signal_sampling_levels"] == list(WORLD_Q_SIGNAL_LEVELS)
    assert sampling["world_dimension_order_single_response"] == [name for name, _, _ in WORLD_AXIS_BOUNDS["single"]]
    assert sampling["world_dimension_order_mixture_response"] == [name for name, _, _ in WORLD_AXIS_BOUNDS["mixture"]]
    assert world_axes["log_r_shift"] == list(WORLD_AXIS_BOUNDS["single"][0][1:])
    assert world_axes["log_f0_shift"] == list(WORLD_AXIS_BOUNDS["single"][2][1:])
    assert world_axes["B0_center"] == list(WORLD_AXIS_BOUNDS["single"][3][1:])
    np.testing.assert_allclose(
        world_axes["m_scale_single_response_factors"],
        [float(np.exp(WORLD_AXIS_BOUNDS["single"][1][1])), float(np.exp(WORLD_AXIS_BOUNDS["single"][1][2]))],
        rtol=0.0,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        world_axes["m_scale_nonresponder_mixture_factors"],
        [float(np.exp(WORLD_AXIS_BOUNDS["mixture"][1][1])), float(np.exp(WORLD_AXIS_BOUNDS["mixture"][1][2]))],
        rtol=0.0,
        atol=1e-15,
    )
    assert world_axes["pi_response_mixture"] == list(WORLD_AXIS_BOUNDS["mixture"][4][1:])
    assert world_axes["mu_nu_day_minus1_levels"] == list(WORLD_MU_NU_LEVELS)
    assert patient_axes["r_day_minus1_log_uniform"] == list(tc.RANGES["r"][:2])
    assert patient_axes["m_log_uniform"] == list(tc.RANGES["m"][:2])
    assert patient_axes["f0_log_uniform"] == list(tc.RANGES["f0"][:2])
    assert patient_axes["f0_after_world_shift_cap"] == 0.9
    assert patient_axes["B0_center_half_width"] == B0_CENTER_HALF_WIDTH
    assert patient_axes["c_uniform"] == list(tc.RANGES["c"][:2])
    assert patient_axes["epsilon_uniform"] == list(tc.RANGES["eps"][:2])
    assert patient_axes["mixture_responder_m"] == list(tc.RANGES["m"][:2])
    assert patient_axes["sigma_balanced_levels"] == list(SIGMA_LEVELS)
    assert "pi_response" in sampling["world_dimension_order_mixture_response"]
    assert "PCG64DXSM" in randomness["generator"]
    assert "calibration_joint_signal" in randomness["streams"]
    assert "evaluation_joint_signal" in randomness["streams"]
    assert "sigma_assignment_child" in randomness["streams"]
    assert randomness["sigma_assignment_child"] == 0x5349474D
    assert "0x5349474D" in sampling["observation_noise_assignment"]
    assert "molecular_signal" not in randomness["streams"]
