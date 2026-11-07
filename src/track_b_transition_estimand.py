"""Synthetic-development calculation for the proposed Duke transition estimand.

Inputs are already registered, history-defined state and distance-band arrays.
This helper does not establish eligibility, construct masks, or authorize scoring.
State codes are 0 outside risk, 1 newly appearing, and 2 persistent. Band codes
0--3 represent the four prespecified distance bands. Target codes are 0--4.
"""

import numpy as np


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
    return {'difference': sum(row['weight'] * row['difference'] for row in rows), 'bands': rows}


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
