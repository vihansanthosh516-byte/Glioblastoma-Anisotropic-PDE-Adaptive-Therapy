import numpy as np
import pytest

from src.trackc_integrated_signal_model import (
    _calibrate_mri_score,
    _uncalibrate_mri_score,
    generate_integrated_signals,
)


def test_signal_draw_is_deterministic_and_policy_view_has_no_truth_fields():
    kwargs = dict(n=100, q=0.72, a=0.5, b=-0.5, d=0.0, seed=20261105)
    observed_a, truth_a = generate_integrated_signals(**kwargs)
    observed_b, truth_b = generate_integrated_signals(**kwargs)

    np.testing.assert_array_equal(observed_a.molecular_report, observed_b.molecular_report)
    np.testing.assert_array_equal(observed_a.mri_score, observed_b.mri_score)
    np.testing.assert_array_equal(truth_a.resistance_z, truth_b.resistance_z)
    assert set(observed_a.__dataclass_fields__) == {"molecular_report", "mri_score"}
    assert not hasattr(observed_a, "resistance_z")


@pytest.mark.parametrize("q", [0.45, 0.69, 0.95])
def test_a5_value_sets_world_level_correct_report_probability(q):
    observed, truth = generate_integrated_signals(
        n=100_000, q=q, a=0.0, b=0.0, d=0.0, seed=20261106
    )
    accuracy = np.mean(observed.molecular_report == truth.molecular_state)
    sensitivity = np.mean(observed.molecular_report[truth.molecular_state])
    specificity = np.mean(~observed.molecular_report[~truth.molecular_state])
    assert abs(accuracy - q) < 0.01
    assert abs(sensitivity - q) < 0.01
    assert abs(specificity - q) < 0.01
    assert abs(np.mean(truth.molecular_state) - 0.5) < 0.01


def test_gaussian_copula_dependence_matches_frozen_factor_formula():
    values = (-0.5, 0.0, 0.5)
    n = 100_000
    case = 0
    for a in values:
        for b in values:
            for d in values:
                observed, truth = generate_integrated_signals(
                    n=n, q=0.75, a=a, b=b, d=d, seed=20261110 + case
                )
                empirical = np.corrcoef(
                    np.vstack(
                        [
                            truth.resistance_z,
                            truth.molecular_state_z,
                            truth.mri_latent_z,
                        ]
                    )
                )
                expected_bio_mri = (
                    a * b + np.sqrt((1.0 - a * a) * (1.0 - b * b)) * d
                )
                expected = np.array(
                    [[1.0, a, b], [a, 1.0, expected_bio_mri], [b, expected_bio_mri, 1.0]]
                )
                assert np.max(np.abs(empirical - expected)) < 0.015
                assert np.linalg.eigvalsh(expected).min() > -1e-10
                assert observed.molecular_report.shape == (n,)
                case += 1


@pytest.mark.parametrize(
    ("a", "b", "d"),
    [(-1.0, -1.0, -1.0), (1.0, 1.0, 1.0), (-1.0, 1.0, 0.0), (1.0, -1.0, 0.5)],
)
def test_correlation_boundary_settings_remain_finite_and_valid(a, b, d):
    observed, truth = generate_integrated_signals(
        n=100, q=0.7, a=a, b=b, d=d, seed=20261111
    )
    assert np.isfinite(truth.resistance_z).all()
    assert np.isfinite(truth.molecular_state_z).all()
    assert np.isfinite(truth.mri_latent_z).all()
    assert np.isfinite(observed.mri_score).all()


def test_history_template_sets_mri_score_location_and_scale():
    template = np.array([-3.0, -1.0, 0.0, 2.0, 5.0])
    observed, _ = generate_integrated_signals(
        n=20_000, q=0.7, a=0.0, b=0.0, d=0.0, seed=20261107,
        mri_template=template,
    )
    assert np.isfinite(observed.mri_score).all()
    assert abs(np.mean(observed.mri_score) - np.mean(template)) < 0.1
    assert abs(np.std(observed.mri_score, ddof=1) - np.std(template, ddof=1)) < 0.1


@pytest.mark.parametrize(
    "template",
    [np.array([-4.0, -1.0, 0.0, 2.0, 7.0]), np.array([0.0, 0.0, 1.0])],
)
def test_history_template_transform_round_trips_even_with_ties(template):
    latent = np.array([-1.5, -0.25, 0.0, 0.4, 1.7])
    scores = _calibrate_mri_score(latent, template)

    np.testing.assert_allclose(
        _uncalibrate_mri_score(scores, template), latent, rtol=0.0, atol=1e-14
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n": 0, "q": 0.7, "a": 0, "b": 0, "d": 0, "seed": 1},
        {"n": True, "q": 0.7, "a": 0, "b": 0, "d": 0, "seed": 1},
        {"n": 10, "q": [0.6, 0.8], "a": 0, "b": 0, "d": 0, "seed": 1},
        {"n": 10, "q": 1.1, "a": 0, "b": 0, "d": 0, "seed": 1},
        {"n": 10, "q": 0.7, "a": 1.1, "b": 0, "d": 0, "seed": 1},
        {"n": 10, "q": 0.7, "a": 0, "b": 0, "d": np.nan, "seed": 1},
        {"n": 10, "q": 0.7, "a": 0, "b": 0, "d": 0, "seed": -1},
    ],
)
def test_invalid_settings_are_rejected(kwargs):
    with pytest.raises(ValueError):
        generate_integrated_signals(**kwargs)


@pytest.mark.parametrize(
    "template",
    [
        np.array([]),
        np.array([1.0]),
        np.array([1.0, 1.0]),
        np.array([1.0, np.inf]),
    ],
)
def test_invalid_mri_template_is_rejected(template):
    with pytest.raises(ValueError, match="mri_template"):
        generate_integrated_signals(
            n=5, q=0.7, a=0, b=0, d=0, seed=1, mri_template=template
        )

