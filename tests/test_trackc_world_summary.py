import numpy as np
import pytest

from trackc_world_summary import (
    paired_world_summary,
    simultaneous_max_t_world_intervals,
    summarize_benchmark_worlds,
)


def test_equal_world_summary_and_seeded_interval():
    candidate = np.array([[1.0, 3.0], [10.0, 14.0], [5.0, 7.0]])
    comparator = np.array([[0.0, 0.0], [8.0, 8.0], [5.0, 5.0]])
    args = (candidate, comparator, ["a", "b", "c"])
    first = paired_world_summary(*args, resamples=2000, seed=11)
    second = paired_world_summary(*args, resamples=2000, seed=11)
    np.testing.assert_allclose(first["world_mean_differences"], [2, 4, 1])
    assert first["equal_world_mean_difference"] == pytest.approx(7 / 3)
    np.testing.assert_array_equal(first["bootstrap_percentile_95"], second["bootstrap_percentile_95"])
    assert first["bootstrap_unit"] == "world"
    assert first["rng"] == "numpy.random.Generator(PCG64DXSM)"


@pytest.mark.parametrize("candidate,comparator,ids", [
    ([[1]], [[0]], ["a"]),
    ([[1], [2]], [[0]], ["a", "b"]),
    ([[1], [np.nan]], [[0], [0]], ["a", "b"]),
    ([[1], [2]], [[0], [0]], ["a", "a"]),
])
def test_invalid_pairs_are_rejected(candidate, comparator, ids):
    with pytest.raises(ValueError):
        paired_world_summary(candidate, comparator, ids, resamples=10, seed=1)


@pytest.mark.parametrize("seed", [-1, True, 1.5])
def test_invalid_world_bootstrap_seed_is_rejected(seed):
    with pytest.raises(ValueError, match="seed must be a nonnegative integer"):
        paired_world_summary([[1.0], [2.0]], [[0.0], [0.0]], ["a", "b"], resamples=10, seed=seed)


def test_benchmark_worlds_reject_duplicate_inputs_and_inconsistent_contrast():
    def world(i):
        return {
            "analysis_type": "synthetic_virtual_world_development", "world_id": str(i),
            "world_provenance": {"family": "single", "scan_loss": None, "world_generation_seed": 400 + i,
                                 "world_parameters_sha256": format(1000 + i, "064x")},
            "settings": {"stop": "observed"}, "inference_model_parameters": [{"d": 0}],
            "signal_generation": {"generator": "test", "parameters": {"d": 0.2}, "seed": 100 + i,
                                  "f0_source": "evaluation_scenarios['f0']", "mri_template_sha256": None},
            "calibration_signal_generation": {"generator": "test", "parameters": {"d": 0.2},
                                              "f0_source": "calibration_scenarios['f0']"},
            "input_fingerprints": {name: format(i * 3 + j, "064x") for j, name in enumerate(
                ("selection_scenarios_sha256", "calibration_scenarios_sha256", "evaluation_scenarios_sha256"))},
            "seeds": {"selection": [200 + i], "evaluation": [300 + i]},
            "arm_outcomes": {"loo_dependence_aware_mpc": {"robust_utility": [1, 2]},
                             "loo_conditional_independence_mpc": {"robust_utility": [0, 0]},
                             "finite_set_minimax_mpc": {"robust_utility": [3, 4]},
                             "finite_set_d0_minimax_mpc": {"robust_utility": [1, 1]}},
            "primary_contrast_name": "minimax",
            "paired_contrasts": {
                "loo": {"contrast": "loo_dependence_aware_mpc minus loo_conditional_independence_mpc",
                        "candidate_arm": "loo_dependence_aware_mpc", "control_arm": "loo_conditional_independence_mpc",
                        "robust_utility_difference": [1, 2]},
                "minimax": {"contrast": "finite_set_minimax_mpc minus finite_set_d0_minimax_mpc",
                            "candidate_arm": "finite_set_minimax_mpc", "control_arm": "finite_set_d0_minimax_mpc",
                            "robust_utility_difference": [2, 3]},
            },
            "primary_paired_contrast": {"contrast": "finite_set_minimax_mpc minus finite_set_d0_minimax_mpc",
                                        "candidate_arm": "finite_set_minimax_mpc", "control_arm": "finite_set_d0_minimax_mpc",
                                        "robust_utility_difference": [2, 3]},
            "n_evaluation_patients": 2,
        }

    first, second = world(1), world(2)
    result = summarize_benchmark_worlds([first, second], resamples=100, seed=3)
    assert result["equal_world_mean_difference"] == pytest.approx(2.5)
    assert result["primary_contrast_name"] == "minimax"
    assert result["weight_source"] == "fixed finite model sets"
    loo = summarize_benchmark_worlds([first, second], resamples=100, seed=3, contrast="loo")
    assert loo["equal_world_mean_difference"] == pytest.approx(1.5)
    assert loo["weight_source"] == "evaluation cohort peers"
    minimax = summarize_benchmark_worlds([first, second], resamples=100, seed=3, contrast="minimax")
    assert minimax["equal_world_mean_difference"] == pytest.approx(2.5)
    assert minimax["weight_source"] == "fixed finite model sets"
    assert minimax["primary_contrast_name"] == "minimax"
    mixed_primary = world(2)
    mixed_primary["primary_contrast_name"] = "loo"
    with pytest.raises(ValueError, match="inconsistent primary contrast names"):
        summarize_benchmark_worlds([first, mixed_primary], resamples=100, seed=3)
    second["world_provenance"]["world_parameters_sha256"] = first["world_provenance"]["world_parameters_sha256"]
    with pytest.raises(ValueError, match="world parameter fingerprint reused"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3, contrast="minimax")
    second = world(2)
    second["world_provenance"]["family"] = "mixture"
    with pytest.raises(ValueError, match="world families differ"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3, contrast="minimax")
    second = world(2)
    second["world_provenance"]["scan_loss"] = {"miss_p": 0.2}
    with pytest.raises(ValueError, match="scan-loss settings differ"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3, contrast="minimax")
    second = world(2)
    del second["arm_outcomes"]["finite_set_d0_minimax_mpc"]
    with pytest.raises(ValueError, match="required minimax arms"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3, contrast="minimax")
    second = world(2)
    second["input_fingerprints"]["evaluation_scenarios_sha256"] = first["input_fingerprints"]["evaluation_scenarios_sha256"]
    with pytest.raises(ValueError, match="fingerprint reused"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3)
    second = world(2)
    second["paired_contrasts"]["loo"]["robust_utility_difference"] = [2, 2]
    with pytest.raises(ValueError, match="recorded loo paired difference"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3, contrast="loo")
    second = world(2)
    second["signal_generation"]["mri_template_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="mri_template_sha256 differs"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3)
    second = world(2)
    second["calibration_signal_generation"]["parameters"] = {"d": 0.9}
    with pytest.raises(ValueError, match="calibration and evaluation"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3)


def test_max_t_world_intervals_are_seeded_and_order_independent():
    single = np.linspace(8.0, 12.0, 40)
    mixture = np.linspace(5.0, 9.0, 40)
    first = simultaneous_max_t_world_intervals(
        {"single-response": single, "mixture-response": mixture},
        resamples=2000,
        seed=20261108,
        precision_target_half_width=14,
    )
    reversed_order = simultaneous_max_t_world_intervals(
        {"mixture-response": mixture, "single-response": single},
        resamples=2000,
        seed=20261108,
        precision_target_half_width=14,
    )

    assert first["method"] == "world-resampled observed-SE max-T bootstrap"
    assert first["decision"] == "pass_if_adopted"
    assert first["all_adjusted_lower_bounds_above_zero"] is True
    assert first["precision_pass"] is True
    assert first["rng"].startswith("numpy.random.Generator(PCG64DXSM)")
    assert "independent within and across" in first["dependence_assumption"]
    assert first["precision_target_half_width_unit"] == "simulated days"
    assert first["max_t_critical_value"] == reversed_order["max_t_critical_value"]
    assert first["intervals_by_family"] == reversed_order["intervals_by_family"]


def test_max_t_requires_both_families_to_pass_positive_rule():
    result = simultaneous_max_t_world_intervals(
        {
            "single-response": np.linspace(8.0, 12.0, 40),
            "mixture-response": np.linspace(-6.0, -2.0, 40),
        },
        resamples=1000,
        seed=7,
    )
    assert result["precision_pass"] is True
    assert result["all_adjusted_lower_bounds_above_zero"] is False
    assert result["decision"] == "does_not_pass_positive_rule"


def test_max_t_marks_failed_precision_target_inconclusive():
    result = simultaneous_max_t_world_intervals(
        {
            "single-response": np.linspace(-100.0, 100.0, 40),
            "mixture-response": np.linspace(-90.0, 110.0, 40),
        },
        resamples=1000,
        seed=8,
        precision_target_half_width=14,
    )
    assert result["precision_pass"] is False
    assert result["decision"] == "inconclusive_precision_target_missed"
    assert all(not row["precision_pass"] for row in result["intervals_by_family"].values())


@pytest.mark.parametrize(
    "families,kwargs",
    [
        ({"single-response": [1.0, 2.0]}, {"resamples": 10, "seed": 1}),
        ({"single-response": [1.0, 2.0], "mixture-response": [1.0, 2.0]}, {"resamples": True, "seed": 1}),
        ({"single-response": [1.0, 2.0], "mixture-response": [1.0, 2.0]}, {"resamples": 10, "seed": -1}),
        ({"single-response": [1.0, 2.0], "mixture-response": [1.0, 2.0]}, {"resamples": 10, "seed": 1, "confidence_level": 1.0}),
        ({"single-response": [1.0], "mixture-response": [1.0, 2.0]}, {"resamples": 10, "seed": 1}),
        ({"single-response": [1.0, 2.0], "mixture-response": [1.0, np.nan]}, {"resamples": 10, "seed": 1}),
        ({"single-response": [1.0, 1.0], "mixture-response": [1.0, 2.0]}, {"resamples": 10, "seed": 1}),
        ({"single-response": [1.0, 2.0], 3: [1.0, 2.0]}, {"resamples": 10, "seed": 1}),
        ({"single-response": [1.0, 2.0], "other": [1.0, 2.0]}, {"resamples": 10, "seed": 1}),
        ({"single-response": [1.0, 2.0], "mixture-response": [1.0, 2.0], "other": [1.0, 2.0]}, {"resamples": 10, "seed": 1}),
    ],
)
def test_max_t_rejects_invalid_inputs(families, kwargs):
    with pytest.raises(ValueError):
        simultaneous_max_t_world_intervals(families, **kwargs)
