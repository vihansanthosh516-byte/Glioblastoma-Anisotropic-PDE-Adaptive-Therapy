"""Synthetic-development calculation for the proposed Duke transition estimand.

History labels must already be registered to the same 2-mm grid. This helper
does not perform registration, assess image QC, or authorize scoring.
It implements the older 2-mm draft, not the distinct 10-mm private Duke run.
State codes are 0 outside risk, 1 newly appearing, and 2 persistent. Band codes
0--3 represent the four prespecified distance bands. Target codes are 0--4.
"""

import numpy as np
from scipy.ndimage import affine_transform, distance_transform_edt


GRID_SPACING_MM = (2.0, 2.0, 2.0)


def construct_history_states(previous_labels, middle_labels, *, voxel_spacing_mm):
    """Build draft risk states and distance bands on a registered 2-mm grid.

    Only history labels are accepted. No support, registration, or QC is inferred.
    ``voxel_spacing_mm`` must be the grid's actual spacing and must equal 2 mm on every
    axis; distances below are in millimetres on that grid. An empty previous-visit
    FLAIR mask is refused, because it would label every risk voxel as new.
    """
    spacing = tuple(float(s) for s in voxel_spacing_mm)
    if len(spacing) != 3 or not np.allclose(spacing, GRID_SPACING_MM, rtol=0, atol=1e-6):
        raise ValueError(f'voxel spacing {spacing} is not the required 2-mm grid')
    previous, middle = map(np.asarray, (previous_labels, middle_labels))
    if previous.ndim != 3 or previous.shape != middle.shape or not previous.size:
        raise ValueError('history labels must have the same nonempty 3D shape')
    if not np.isin(previous, [0, 1, 2, 3, 4]).all() or not np.isin(middle, [0, 1, 2, 3, 4]).all():
        raise ValueError('invalid history label')
    core = np.isin(middle, [1, 3])
    if not core.any():
        raise ValueError('empty core')
    core_distance = distance_transform_edt(~core, sampling=2.0)
    excluded = core | (middle == 4)
    excluded_distance = distance_transform_edt(~excluded, sampling=2.0)
    risk = (middle == 2) & (excluded_distance > 2.0)
    previous_flair = previous == 2
    if not previous_flair.any():
        raise ValueError('empty previous FLAIR mask')
    near_previous = distance_transform_edt(~previous_flair, sampling=2.0) <= 2.0
    state = np.zeros(middle.shape, dtype=np.uint8)
    state[risk & ~near_previous] = 1
    state[risk & near_previous] = 2
    bands = np.searchsorted([10.0, 20.0, 30.0], core_distance, side='left').astype(np.uint8)
    return state, bands


def history_band_counts(state, bands, *, min_voxels=10):
    """Identify eligible bands using history-derived arrays only."""
    state, bands = map(np.asarray, (state, bands))
    if state.shape != bands.shape or state.size == 0:
        raise ValueError('inputs must have the same nonempty shape')
    if isinstance(min_voxels, bool) or not isinstance(min_voxels, int) or min_voxels < 1:
        raise ValueError('min_voxels must be a positive integer')
    if not np.isin(state, [0, 1, 2]).all() or not np.isin(bands, [0, 1, 2, 3]).all():
        raise ValueError('invalid state or distance band')
    result = []
    for band in range(4):
        n_new = int(np.count_nonzero((state == 1) & (bands == band)))
        n_persistent = int(np.count_nonzero((state == 2) & (bands == band)))
        if min(n_new, n_persistent) >= min_voxels:
            result.append({'band': band, 'n_new': n_new, 'n_persistent': n_persistent})
    return result


def output_to_input_affine(moving_affine, fixed_affine, moving_to_fixed):
    """Map fixed output indices to moving input indices via physical coordinates.

    Affines map voxel centers to a shared physical coordinate system; the
    registration transform maps moving physical points to fixed physical points.
    """
    matrices = [np.asarray(value, dtype=float) for value in
                (moving_affine, fixed_affine, moving_to_fixed)]
    for matrix in matrices:
        if (matrix.shape != (4, 4) or not np.isfinite(matrix).all()
                or not np.allclose(matrix[3], [0, 0, 0, 1])
                or np.linalg.matrix_rank(matrix[:3, :3]) < 3):
            raise ValueError('invalid homogeneous affine')
    moving, fixed, transform = matrices
    return np.linalg.solve(moving, np.linalg.solve(transform, fixed))


def geometric_coverage(moving_shape, output_shape, output_to_input):
    """Resample a geometric all-ones field using scipy nearest-neighbor rules.

    The matrix maps output indices to input indices, not physical coordinates.
    The caller must supply the exact transform and grids used for label resampling.
    """
    def valid_shape(shape):
        return (len(shape) == 3 and all(isinstance(x, (int, np.integer))
                and not isinstance(x, (bool, np.bool_)) and x > 0 for x in shape))
    if not valid_shape(moving_shape) or not valid_shape(output_shape):
        raise ValueError('shapes must contain three positive integers')
    matrix = np.asarray(output_to_input, dtype=float)
    if (matrix.shape != (4, 4) or not np.isfinite(matrix).all()
            or not np.allclose(matrix[3], [0, 0, 0, 1])
            or np.linalg.matrix_rank(matrix[:3, :3]) < 3):
        raise ValueError('invalid output-to-input affine')
    return affine_transform(np.ones(moving_shape, dtype=np.uint8), matrix[:3, :3],
                            offset=matrix[:3, 3], output_shape=tuple(output_shape),
                            order=0, mode='constant', cval=0, prefilter=False).astype(bool)


def target_field_coverage(state, bands, target_valid, *, min_voxels=10):
    """Count risk voxels outside target geometric coverage; no target labels used.

    Caller supplies an independently computed coverage mask from the target grid
    and registration transform. This diagnostic does not validate that mask.
    """
    state, bands, target_valid = map(np.asarray, (state, bands, target_valid))
    eligible = history_band_counts(state, bands, min_voxels=min_voxels)
    if target_valid.shape != state.shape or target_valid.dtype != np.bool_:
        raise ValueError('target_valid must be a same-shape boolean mask')
    return [{**row,
             'new_outside': int(np.count_nonzero((state == 1) & (bands == row['band']) & ~target_valid)),
             'persistent_outside': int(np.count_nonzero((state == 2) & (bands == row['band']) & ~target_valid))}
            for row in eligible]


def summarize_window(state, bands, target, *, min_voxels=10):
    """Return matched-band contrast and target-label counts, or None if ineligible."""
    state, bands, target = map(np.asarray, (state, bands, target))
    eligible = history_band_counts(state, bands, min_voxels=min_voxels)
    if state.shape != target.shape:
        raise ValueError('inputs must have the same nonempty shape')
    if not np.isin(target, [0, 1, 2, 3, 4]).all():
        raise ValueError('invalid target label')
    rows = []
    for entry in eligible:
        band = entry['band']
        new = (state == 1) & (bands == band)
        persistent = (state == 2) & (bands == band)
        n_new, n_persistent = entry['n_new'], entry['n_persistent']
        def counts(mask):
            return {label: int(np.count_nonzero(target[mask] == label)) for label in range(5)}
        new_labels, persistent_labels = counts(new), counts(persistent)
        delta = ((new_labels[1] + new_labels[3]) / n_new
                 - (persistent_labels[1] + persistent_labels[3]) / n_persistent)
        rows.append({'band': band, 'n_new': n_new, 'n_persistent': n_persistent,
                     'new_labels': new_labels, 'persistent_labels': persistent_labels,
                     'difference': delta, 'matched_count': min(n_new, n_persistent)})
    if not rows:
        return None
    total = sum(row['matched_count'] for row in rows)
    for row in rows:
        row['weight'] = row['matched_count'] / total
    incomplete = []
    sensitivity = 0.0
    for row in rows:
        new_count = row['n_new'] - row['new_labels'][0]
        persistent_count = row['n_persistent'] - row['persistent_labels'][0]
        if not new_count or not persistent_count:
            incomplete.append(row['band'])
            continue
        new_rate = (row['new_labels'][1] + row['new_labels'][3]) / new_count
        persistent_rate = (row['persistent_labels'][1] + row['persistent_labels'][3]) / persistent_count
        sensitivity += row['weight'] * (new_rate - persistent_rate)
    return {'difference': sum(row['weight'] * row['difference'] for row in rows), 'bands': rows,
            'label_zero_sensitivity': {'difference': None if incomplete else sensitivity,
                                       'incomplete_bands': incomplete}}


def summarize_patients(patient_windows):
    """Average window contrasts within patients and then equally across patients."""
    if not patient_windows:
        raise ValueError('no eligible patients')
    patient_means = []
    for values in patient_windows.values():
        values = np.asarray(values, dtype=float)
        if values.ndim != 1 or not values.size or not np.isfinite(values).all() or (np.abs(values) > 1).any():
            raise ValueError('invalid patient window contrasts')
        patient_means.append(float(values.mean()))
    return float(np.mean(patient_means))


def bootstrap_patients(patient_windows, *, n_resamples=10000, seed=20261104):
    """Percentile bootstrap of equally weighted patient means, not windows."""
    if isinstance(n_resamples, bool) or not isinstance(n_resamples, int) or n_resamples < 1:
        raise ValueError('n_resamples must be a positive integer')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    mean = summarize_patients(patient_windows)
    patient_means = np.array([summarize_patients({key: patient_windows[key]})
                              for key in sorted(patient_windows)])
    rng = np.random.default_rng(seed)
    replicates = np.empty(n_resamples)
    for i in range(n_resamples):
        replicates[i] = patient_means[rng.integers(0, len(patient_means), len(patient_means))].mean()
    lo, hi = np.quantile(replicates, [0.025, 0.975], method='linear')
    return {'mean': mean, 'interval': (float(lo), float(hi)),
            'n_patients': len(patient_means), 'n_resamples': n_resamples, 'seed': seed}
