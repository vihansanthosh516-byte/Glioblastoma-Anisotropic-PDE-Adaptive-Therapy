import numpy as np
import pytest

from trackc_world_summary import paired_world_summary, summarize_benchmark_worlds


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


@pytest.mark.parametrize("candidate,comparator,ids", [
    ([[1]], [[0]], ["a"]),
    ([[1], [2]], [[0]], ["a", "b"]),
    ([[1], [np.nan]], [[0], [0]], ["a", "b"]),
    ([[1], [2]], [[0], [0]], ["a", "a"]),
])
def test_invalid_pairs_are_rejected(candidate, comparator, ids):
    with pytest.raises(ValueError):
        paired_world_summary(candidate, comparator, ids, resamples=10, seed=1)


def test_benchmark_worlds_reject_duplicate_inputs_and_inconsistent_contrast():
    def world(i):
        return {
            "analysis_type": "synthetic_virtual_world_development", "world_id": str(i),
            "settings": {"stop": "observed"}, "inference_model_parameters": [{"d": 0}],
            "signal_generation": {"generator": "test", "parameters": {"d": 0.2}, "seed": 100 + i,
                                  "f0_source": "evaluation_scenarios['f0']", "mri_template_sha256": None},
            "calibration_signal_generation": {"generator": "test", "parameters": {"d": 0.2},
                                              "f0_source": "calibration_scenarios['f0']"},
            "input_fingerprints": {name: format(i * 3 + j, "064x") for j, name in enumerate(
                ("selection_scenarios_sha256", "calibration_scenarios_sha256", "evaluation_scenarios_sha256"))},
            "seeds": {"selection": [200 + i], "evaluation": [300 + i]},
            "arm_outcomes": {"loo_dependence_aware_mpc": {"robust_utility": [1, 2]},
                             "loo_conditional_independence_mpc": {"robust_utility": [0, 0]}},
            "primary_paired_contrast": {"contrast": "loo_dependence_aware_mpc minus loo_conditional_independence_mpc",
                                        "robust_utility_difference": [1, 2]},
            "n_evaluation_patients": 2,
        }

    first, second = world(1), world(2)
    result = summarize_benchmark_worlds([first, second], resamples=100, seed=3)
    assert result["equal_world_mean_difference"] == pytest.approx(1.5)
    second["input_fingerprints"]["evaluation_scenarios_sha256"] = first["input_fingerprints"]["evaluation_scenarios_sha256"]
    with pytest.raises(ValueError, match="fingerprint reused"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3)
    second = world(2)
    second["primary_paired_contrast"]["robust_utility_difference"] = [2, 2]
    with pytest.raises(ValueError, match="primary contrast"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3)
    second = world(2)
    second["signal_generation"]["mri_template_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="mri_template_sha256 differs"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3)
    second = world(2)
    second["calibration_signal_generation"]["parameters"] = {"d": 0.9}
    with pytest.raises(ValueError, match="calibration and evaluation"):
        summarize_benchmark_worlds([first, second], resamples=100, seed=3)
