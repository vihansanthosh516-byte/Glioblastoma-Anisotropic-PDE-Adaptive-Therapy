"""Synthetic-development calculation for the proposed Duke transition estimand.

Inputs are already registered, history-defined state and distance-band arrays.
This helper does not establish eligibility, construct masks, or authorize scoring.
State codes are 0 outside risk, 1 newly appearing, and 2 persistent. Band codes
0--3 represent the four prespecified distance bands. Target codes are 0--4.
"""

import numpy as np


def summarize_window(state, bands, target, *, min_voxels=10):
    """Return matched-band contrast and target-label counts, or None if ineligible."""
    state, bands, target = map(np.asarray, (state, bands, target))
    if state.shape != bands.shape or state.shape != target.shape or state.size == 0:
        raise ValueError('inputs must have the same nonempty shape')
    if isinstance(min_voxels, bool) or not isinstance(min_voxels, int) or min_voxels < 1:
        raise ValueError('min_voxels must be a positive integer')
    if not np.isin(state, [0, 1, 2]).all() or not np.isin(bands, [0, 1, 2, 3]).all():
        raise ValueError('invalid state or distance band')
    if not np.isin(target, [0, 1, 2, 3, 4]).all():
        raise ValueError('invalid target label')
    rows = []
    for band in range(4):
        new = (state == 1) & (bands == band)
        persistent = (state == 2) & (bands == band)
        n_new, n_persistent = int(new.sum()), int(persistent.sum())
        if min(n_new, n_persistent) < min_voxels:
            continue
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
