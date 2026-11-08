import numpy as np
import pytest

from src.track_b_transition_estimand import (construct_history_states, history_band_counts, geometric_coverage, output_to_input_affine, target_field_coverage, summarize_window,
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


def test_history_states_use_physical_distance_and_cavity():
    previous = np.zeros((24, 3, 3), dtype=int)
    middle = np.zeros_like(previous)
    middle[1, 1, 1] = 3
    middle[5:24, 1, 1] = 2
    middle[9, 2, 1] = 4
    previous[5, 1, 1] = 2
    state, bands = construct_history_states(previous, middle)
    assert state[5, 1, 1] == 2
    assert state[6, 1, 1] == 2  # within 2 mm of prior FLAIR
    assert state[7, 1, 1] == 1
    assert state[9, 1, 1] == 0  # within 2 mm of cavity
    assert bands[6, 1, 1] == 0  # exactly 10 mm
    assert bands[11, 1, 1] == 1  # exactly 20 mm
    assert bands[16, 1, 1] == 2  # exactly 30 mm
    assert bands[17, 1, 1] == 3


def test_empty_previous_flair_and_empty_core():
    previous = np.zeros((5, 5, 5), dtype=int)
    middle = previous.copy()
    middle[0, 0, 0] = 1
    middle[3, 3, 3] = 2
    state, _ = construct_history_states(previous, middle)
    assert state[3, 3, 3] == 1
    middle[0, 0, 0] = 0
    with pytest.raises(ValueError, match='empty core'):
        construct_history_states(previous, middle)


def test_invalid_history_label_and_shape_rejected():
    with pytest.raises(ValueError):
        construct_history_states(np.zeros((2, 2, 2)), np.full((2, 2, 2), 5))
    with pytest.raises(ValueError):
        construct_history_states(np.zeros((2, 2)), np.zeros((2, 2)))


def test_synthetic_history_to_window_score():
    previous = np.zeros((8, 8, 8), dtype=int)
    middle = previous.copy()
    middle[0, 0, 0] = 1
    middle[3, 3, 3] = middle[3, 3, 4] = 2
    previous[3, 3, 3] = 2
    # Move the new site away from the 2-mm prior dilation.
    middle[3, 3, 4] = 0
    middle[3, 5, 3] = 2
    state, bands = construct_history_states(previous, middle)
    target = np.zeros_like(previous)
    target[3, 5, 3] = 3
    result = summarize_window(state, bands, target, min_voxels=1)
    assert result['difference'] == 1.0
    assert result['bands'][0]['n_new'] == result['bands'][0]['n_persistent'] == 1


def test_descriptive_label_zero_sensitivity_keeps_history_weights():
    state = np.array([1, 1, 2, 2, 1, 2])
    bands = np.array([0, 0, 0, 0, 1, 1])
    target = np.array([1, 0, 0, 2, 0, 3])
    result = summarize_window(state, bands, target, min_voxels=1)
    # New denominator in band 1 becomes zero: sensitivity is incomplete,
    # rather than silently dropping a history-selected band.
    assert result['label_zero_sensitivity']['difference'] is None
    assert result['label_zero_sensitivity']['incomplete_bands'] == [1]


def test_label_zero_sensitivity_valid_difference():
    state = np.array([1, 1, 2, 2])
    bands = np.zeros(4, dtype=int)
    target = np.array([1, 0, 2, 0])
    result = summarize_window(state, bands, target, min_voxels=1)
    assert result['difference'] == pytest.approx(0.5)
    assert result['label_zero_sensitivity']['difference'] == pytest.approx(1.0)
    assert result['label_zero_sensitivity']['incomplete_bands'] == []


def test_target_field_coverage_keeps_distinct_history_groups():
    state = np.array([1, 1, 2, 2, 1, 2])
    bands = np.array([0, 0, 0, 0, 1, 1])
    valid = np.array([True, False, True, True, False, True])
    rows = target_field_coverage(state, bands, valid, min_voxels=1)
    assert rows[0]['new_outside'] == 1
    assert rows[0]['persistent_outside'] == 0
    assert rows[1]['new_outside'] == 1
    assert rows[1]['persistent_outside'] == 0
    assert rows[0]['n_new'] == rows[0]['n_persistent'] == 2


def test_target_field_coverage_refuses_invalid_mask():
    state, bands = np.array([1, 2]), np.array([0, 0])
    for valid in (np.array([1, 1]), np.array([True])):
        with pytest.raises(ValueError):
            target_field_coverage(state, bands, valid, min_voxels=1)


def test_geometric_coverage_identity_and_translation():
    identity = np.eye(4)
    covered = geometric_coverage((4, 3, 2), (4, 3, 2), identity)
    assert covered.dtype == np.bool_
    assert covered.all()
    identity[0, 3] = 1  # fixed output index maps to moving input index + 1
    shifted = geometric_coverage((4, 3, 2), (4, 3, 2), identity)
    assert shifted[:3].all()
    assert not shifted[3].any()


def test_geometric_coverage_refuses_bad_mapping():
    for transform in (np.zeros((4, 4)), np.full((4, 4), np.nan), np.eye(3)):
        with pytest.raises(ValueError):
            geometric_coverage((2, 2, 2), (2, 2, 2), transform)
    with pytest.raises(ValueError):
        geometric_coverage((2, 0, 2), (2, 2, 2), np.eye(4))


def test_physical_affines_map_output_to_moving_indices():
    moving = np.diag([2., 3., 4., 1.])
    moving[:3, 3] = [10., -3., 4.]
    fixed = np.diag([1., 2., 2., 1.])
    fixed[:3, 3] = [12., 1., 8.]
    moving_to_fixed = np.eye(4)
    moving_to_fixed[:3, 3] = [2., 1., -4.]
    mapping = output_to_input_affine(moving, fixed, moving_to_fixed)
    fixed_index = np.array([4., 3., 2., 1.])
    # Fixed world = (16,7,12); inverse motion = (14,6,16);
    # inverse moving affine = (2,3,3).
    np.testing.assert_allclose(mapping @ fixed_index, [2, 3, 3, 1])


def test_physical_affines_reject_singular_and_invalid_homogeneous():
    bad_bottom = np.eye(4)
    bad_bottom[3, 0] = 1
    for bad in (np.zeros((4, 4)), bad_bottom, np.eye(3), np.full((4, 4), np.nan)):
        with pytest.raises(ValueError):
            output_to_input_affine(bad, np.eye(4), np.eye(4))


def test_physical_rotation_uses_inverse_registration():
    moving_to_fixed = np.eye(4)
    moving_to_fixed[:2, :2] = [[0, -1], [1, 0]]
    mapping = output_to_input_affine(np.eye(4), np.eye(4), moving_to_fixed)
    np.testing.assert_allclose(mapping @ [0, 2, 1, 1], [2, 0, 1, 1])


def test_geometric_coverage_matches_nibabel_physical_resampling():
    nib = pytest.importorskip('nibabel')
    from nibabel.processing import resample_from_to

    def rotation(axis, angle):
        c, s = np.cos(angle), np.sin(angle)
        if axis == 'x':
            return np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])
        if axis == 'y':
            return np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]])
        return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])

    moving_shape = (7, 6, 5)
    output_shape = (9, 8, 6)
    moving_affine = np.eye(4)
    moving_affine[:3, :3] = np.diag([1.7, 2.2, 3.0])
    moving_affine[:3, 3] = [-8., 4., 11.]
    fixed_affine = np.eye(4)
    fixed_affine[:3, :3] = rotation('z', .23) @ np.diag([1.3, 1.8, 2.5])
    fixed_affine[:3, 3] = [-4., 2., 8.]
    moving_to_fixed = np.eye(4)
    moving_to_fixed[:3, :3] = rotation('y', .17) @ rotation('x', -.11)
    moving_to_fixed[:3, 3] = [2.1, -1.5, .7]

    cases = [(moving_shape, output_shape, moving_affine, fixed_affine, moving_to_fixed)]

    translated_moving = np.eye(4)
    translated_moving[:3, :3] = np.diag([2., 1.5, 1.])
    translated_moving[:3, 3] = [3., -5., 7.]
    translated_fixed = np.eye(4)
    translated_fixed[:3, 3] = [2., -7., 6.]
    translation = np.eye(4)
    translation[:3, 3] = [.5, 1., -1.]
    cases.append(((6, 7, 8), (8, 7, 6), translated_moving,
                  translated_fixed, translation))

    for src_shape, dst_shape, src_affine, dst_affine, transform in cases:
        output_to_input = output_to_input_affine(src_affine, dst_affine, transform)
        actual = geometric_coverage(src_shape, dst_shape, output_to_input)

        # The registration maps moving-world coordinates to fixed-world coordinates.
        registered = nib.Nifti1Image(
            np.ones(src_shape, dtype=np.uint8), transform @ src_affine)
        resampled = resample_from_to(
            registered, (dst_shape, dst_affine), order=0).get_fdata() > 0.5
        np.testing.assert_array_equal(actual, resampled)
