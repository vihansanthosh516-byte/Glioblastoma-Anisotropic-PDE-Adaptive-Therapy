import numpy as np
import pytest

from src.track_b_transition_estimand import (history_band_counts, summarize_window,
                                              summarize_patients, bootstrap_patients)


def test_weighted_bands_and_label_distribution():
    # Two bands: weights 2/3 and 1/3; differences 1/2 and -1.
    state = np.array([1, 1, 2, 2, 1, 2])
    bands = np.array([0, 0, 0, 0, 1, 1])
    target = np.array([1, 0, 0, 2, 0, 3])
    result = summarize_window(state, bands, target, min_voxels=1)
    assert result['difference'] == pytest.approx(0)
    assert result['bands'][0]['weight'] == pytest.approx(2 / 3)
    assert result['bands'][0]['new_labels'] == {0: 1, 1: 1, 2: 0, 3: 0, 4: 0}


def test_history_eligibility_and_no_conversion():
    state = np.array([1, 2, 1, 1])
    bands = np.array([0, 0, 1, 1])
    target = np.array([0, 4, 3, 3])
    result = summarize_window(state, bands, target, min_voxels=1)
    assert result['difference'] == 0
    assert len(result['bands']) == 1
    assert summarize_window(state, bands, target, min_voxels=2) is None


def test_patient_equal_weight_not_window_equal_weight():
    assert summarize_patients({'a': [1.0, 1.0, 1.0], 'b': [-1.0]}) == 0.0


@pytest.mark.parametrize('state,bands,target', [
    ([1], [0, 1], [0]), ([3], [0], [0]), ([1], [4], [0]),
    ([1], [0], [5]), ([1], [0], [float('nan')]),
])
def test_invalid_input_rejected(state, bands, target):
    with pytest.raises(ValueError):
        summarize_window(np.array(state), np.array(bands), np.array(target), min_voxels=1)


def test_invalid_patient_results_rejected():
    for value in ({}, {'a': []}, {'a': [float('nan')]}, {'a': [1.1]}):
        with pytest.raises(ValueError):
            summarize_patients(value)


def test_history_eligibility_requires_no_target():
    state = np.array([1, 1, 2, 2, 1, 2])
    bands = np.array([0, 0, 0, 0, 1, 1])
    rows = history_band_counts(state, bands, min_voxels=2)
    assert rows == [{'band': 0, 'n_new': 2, 'n_persistent': 2}]
    assert history_band_counts(state, bands, min_voxels=3) == []


def test_bootstrap_resamples_patient_means_not_windows():
    data = {'b': [-1.0], 'a': [1.0, 1.0, 1.0]}
    result = bootstrap_patients(data, n_resamples=10000, seed=20261104)
    assert result['mean'] == 0
    assert result['interval'] == (-1.0, 1.0)
    assert result == bootstrap_patients(dict(reversed(list(data.items()))), n_resamples=10000, seed=20261104)


def test_bootstrap_refuses_invalid_settings():
    for options in ({'n_resamples': 0}, {'n_resamples': True}, {'seed': -1}):
        with pytest.raises(ValueError):
            bootstrap_patients({'a': [0.0]}, **options)
