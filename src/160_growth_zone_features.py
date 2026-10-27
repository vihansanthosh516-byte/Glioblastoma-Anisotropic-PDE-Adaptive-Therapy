#!/usr/bin/env python3
"""Features for an exploratory, patient-level GBM growth-zone forecast.

This module builds only from the previous and current visits.  The future
visit belongs in ``make_transition_target`` and is deliberately not accepted
by the feature builder.  Labels 1 and 3 are treated as annotated core and
enhancing tissue, not as verified viable tumor or clinical recurrence.

The spatial prior is a screened graph-diffusion approximation to
``-div(A grad u) + kappa**2 u = source``.  Edge conductances depend on
multisequence MRI contrast and label-2 tissue, so propagation is spatially
directional and cavity voxels are removed from the graph.  It is a ranking
prior, not a cellular-density estimate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, distance_transform_edt, gaussian_filter, shift
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, cg


MODALITIES = ("t1n", "t1c", "t2w", "t2f")
FEATURE_NAMES = (
    "distance_to_core_mm",
    "distance_to_cavity_mm",
    "t1n_z",
    "t1c_z",
    "t2w_z",
    "t2f_z",
    "label2_indicator",
    "prior_addition_density",
)
ORIGIN_FEATURE_NAMES = tuple(name for name in FEATURE_NAMES if name != "prior_addition_density")


@dataclass(frozen=True)
class VisitVolume:
    """One preprocessed visit on a 3-D NIfTI grid."""

    images: Mapping[str, np.ndarray]
    labels: np.ndarray
    affine: np.ndarray


@dataclass(frozen=True)
class DiffusionParameters:
    """Low-dimensional graph-PDE settings fixed for a run or selected in-fold."""

    edge_sensitivity: float = 0.5
    screening_length_mm: float = 10.0
    label2_log_permeability: float = 0.69314718056
    candidate_radius_mm: float = 20.0
    prior_smoothing_mm: float = 2.0
    block_cavity: bool = True
    solver_tolerance: float = 1e-6
    solver_max_iterations: int = 1000

    def __post_init__(self) -> None:
        if self.edge_sensitivity < 0:
            raise ValueError("edge_sensitivity must be nonnegative")
        if self.screening_length_mm <= 0 or self.candidate_radius_mm <= 0:
            raise ValueError("screening length and candidate radius must be positive")
        if self.solver_tolerance <= 0 or self.solver_max_iterations < 1:
            raise ValueError("invalid iterative-solver settings")


@dataclass(frozen=True)
class OriginDiffusionParameters:
    """Parameters for the whole-brain, one-origin-scan spatial prior."""

    edge_sensitivity: float = 0.5
    screening_length_mm: float = 10.0
    label2_log_permeability: float = 0.69314718056
    block_cavity: bool = True
    solver_tolerance: float = 1e-6
    solver_max_iterations: int = 1000

    def __post_init__(self) -> None:
        if self.edge_sensitivity < 0:
            raise ValueError("edge sensitivity must be nonnegative")
        if self.screening_length_mm <= 0 or self.solver_tolerance <= 0 or self.solver_max_iterations < 1:
            raise ValueError("invalid origin diffusion settings")


@dataclass(frozen=True)
class VisitFeatures:
    """Voxel maps for one forecast origin; no future-visit information."""

    candidate: np.ndarray
    source: np.ndarray
    recent_additions: np.ndarray
    prior_addition_count: int
    prior_candidate_voxels: int
    maps: Mapping[str, np.ndarray]
    pde_field: np.ndarray
    log_pde_field: np.ndarray


def _as_3d_float(array: np.ndarray, name: str) -> np.ndarray:
    value = np.asarray(array, dtype=np.float32)
    if value.ndim != 3:
        raise ValueError(f"{name} must be a 3-D array")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} contains non-finite values")
    return value


def validate_visit(visit: VisitVolume, name: str = "visit") -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    """Check modality, label, and affine contracts without resampling data."""
    if set(visit.images) != set(MODALITIES):
        missing = sorted(set(MODALITIES) - set(visit.images))
        extra = sorted(set(visit.images) - set(MODALITIES))
        raise ValueError(f"{name} modalities mismatch; missing={missing}, extra={extra}")
    labels = np.asarray(visit.labels)
    if labels.ndim != 3:
        raise ValueError(f"{name} labels must be a 3-D array")
    if not np.issubdtype(labels.dtype, np.integer):
        rounded = np.rint(labels)
        if not np.all(np.isfinite(labels)) or not np.allclose(labels, rounded):
            raise ValueError(f"{name} labels must contain integer codes")
        labels = rounded.astype(np.int16)
    labels = labels.astype(np.int16, copy=False)
    unknown = np.setdiff1d(np.unique(labels), np.array([0, 1, 2, 3, 4]))
    if unknown.size:
        raise ValueError(f"{name} contains undocumented label codes: {unknown.tolist()}")
    images = {key: _as_3d_float(visit.images[key], f"{name}.{key}") for key in MODALITIES}
    if any(image.shape != labels.shape for image in images.values()):
        raise ValueError(f"{name} image and label shapes differ")
    affine = np.asarray(visit.affine, dtype=np.float64)
    if affine.shape != (4, 4) or not np.all(np.isfinite(affine)):
        raise ValueError(f"{name} affine must be a finite 4x4 matrix")
    spacing = np.linalg.norm(affine[:3, :3], axis=0)
    if np.any(spacing <= 0) or not np.all(np.isfinite(spacing)):
        raise ValueError(f"{name} has invalid voxel spacing")
    directions = affine[:3, :3] / spacing[np.newaxis, :]
    axis_gram = directions.T @ directions
    if not np.allclose(axis_gram, np.eye(3), rtol=0.0, atol=1e-5):
        raise ValueError(
            f"{name} affine has non-orthogonal voxel axes; "
            "distance and graph operators require an orthogonal grid"
        )
    return images, labels, affine, spacing


def validate_aligned_visits(*visits: VisitVolume, atol_mm: float = 1e-3) -> np.ndarray:
    """Require common voxel arrays and affines; never silently register scans."""
    if len(visits) < 2:
        raise ValueError("at least two visits are required")
    checked = [validate_visit(visit, f"visit[{index}]") for index, visit in enumerate(visits)]
    shape = checked[0][1].shape
    affine = checked[0][2]
    for index, (_, labels, other_affine, _) in enumerate(checked[1:], start=1):
        if labels.shape != shape:
            raise ValueError(f"visit[{index}] grid shape differs; explicit registration is required")
        if not np.allclose(other_affine, affine, rtol=0.0, atol=atol_mm):
            raise ValueError(f"visit[{index}] affine differs; explicit registration is required")
    return checked[0][3]


def annotated_core_mask(labels: np.ndarray) -> np.ndarray:
    """Return labels 1+3, the operational endpoint used by this prototype."""
    return (np.asarray(labels) == 1) | (np.asarray(labels) == 3)


def brain_support(images: Mapping[str, np.ndarray], labels: np.ndarray) -> np.ndarray:
    """Approximate the skull-stripped image support, retaining labeled voxels."""
    nonzero = np.zeros(np.asarray(labels).shape, dtype=bool)
    for modality in MODALITIES:
        image = np.asarray(images[modality])
        nonzero |= np.isfinite(image) & (np.abs(image) > 1e-6)
    return nonzero | (np.asarray(labels) != 0)


def _robust_z(image: np.ndarray, support: np.ndarray) -> np.ndarray:
    values = np.asarray(image, dtype=np.float32)[support]
    if values.size < 100:
        raise ValueError("too few supported voxels to normalize MRI intensities")
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median))) * 1.4826
    if not np.isfinite(mad) or mad <= 1e-6:
        lo, hi = np.percentile(values, [10.0, 90.0])
        mad = float((hi - lo) / 2.563)
    if not np.isfinite(mad) or mad <= 1e-6:
        raise ValueError("MRI modality has no usable intensity variation")
    result = (np.asarray(image, dtype=np.float32) - median) / mad
    return np.clip(result, -8.0, 8.0).astype(np.float32, copy=False)


def _candidate_region(
    labels: np.ndarray,
    support: np.ndarray,
    spacing: Sequence[float],
    radius_mm: float,
    block_cavity: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    source = annotated_core_mask(labels)
    if not np.any(source):
        raise ValueError("current visit has no label-1/3 source voxels")
    distance = distance_transform_edt(~source, sampling=spacing).astype(np.float32)
    cavity = np.asarray(labels) == 4
    in_radius = distance <= radius_mm
    candidate = in_radius & support & ~source & ~cavity
    graph_domain = in_radius & support & (~cavity if block_cavity else np.ones_like(cavity))
    if not np.any(candidate):
        raise ValueError("candidate region is empty after support and cavity masking")
    return source, candidate, graph_domain


def _screened_graph_field(
    images_z: Mapping[str, np.ndarray],
    labels: np.ndarray,
    source: np.ndarray,
    candidate: np.ndarray,
    graph_domain: np.ndarray,
    spacing: Sequence[float],
    parameters: DiffusionParameters,
) -> np.ndarray:
    """Solve a screened, image-guided 6-neighbor graph diffusion problem."""
    flat_domain = np.asarray(graph_domain, dtype=bool)
    node_index = np.full(flat_domain.shape, -1, dtype=np.int32)
    node_index[flat_domain] = np.arange(int(flat_domain.sum()), dtype=np.int32)
    node_count = int(flat_domain.sum())
    if node_count == 0:
        raise ValueError("diffusion graph has no nodes")

    edge_i_parts: list[np.ndarray] = []
    edge_j_parts: list[np.ndarray] = []
    edge_w_parts: list[np.ndarray] = []
    label2 = np.asarray(labels) == 2
    z = np.stack([np.asarray(images_z[key], dtype=np.float32) for key in MODALITIES], axis=0)

    for axis, voxel_size in enumerate(spacing):
        left = [slice(None)] * 3
        right = [slice(None)] * 3
        left[axis] = slice(0, -1)
        right[axis] = slice(1, None)
        left_t, right_t = tuple(left), tuple(right)
        i_grid = node_index[left_t]
        j_grid = node_index[right_t]
        valid = (i_grid >= 0) & (j_grid >= 0)
        if not np.any(valid):
            continue
        i_nodes = i_grid[valid].astype(np.int32, copy=False)
        j_nodes = j_grid[valid].astype(np.int32, copy=False)
        delta = z[(slice(None),) + left_t][:, valid] - z[(slice(None),) + right_t][:, valid]
        squared_difference = np.mean(delta * delta, axis=0, dtype=np.float32)
        conductance = np.exp(-parameters.edge_sensitivity * squared_difference)
        label2_pair = 0.5 * (
            label2[left_t][valid].astype(np.float32) + label2[right_t][valid].astype(np.float32)
        )
        conductance *= np.exp(parameters.label2_log_permeability * label2_pair)
        conductance /= float(voxel_size) ** 2
        conductance = np.maximum(conductance, 1e-8).astype(np.float64, copy=False)
        edge_i_parts.append(i_nodes)
        edge_j_parts.append(j_nodes)
        edge_w_parts.append(conductance)

    if not edge_w_parts:
        raise ValueError("diffusion graph has no adjacent voxel pairs")
    edge_i = np.concatenate(edge_i_parts)
    edge_j = np.concatenate(edge_j_parts)
    edge_w = np.concatenate(edge_w_parts)
    degree = np.bincount(edge_i, weights=edge_w, minlength=node_count)
    degree += np.bincount(edge_j, weights=edge_w, minlength=node_count)
    screening = 1.0 / float(parameters.screening_length_mm) ** 2
    diagonal = degree + screening
    all_nodes = np.arange(node_count, dtype=np.int32)
    row = np.concatenate((edge_i, edge_j, all_nodes))
    col = np.concatenate((edge_j, edge_i, all_nodes))
    values = np.concatenate((-edge_w, -edge_w, diagonal))
    matrix = coo_matrix((values, (row, col)), shape=(node_count, node_count)).tocsr()

    rhs = np.asarray(source[flat_domain], dtype=np.float64)
    preconditioner = LinearOperator(matrix.shape, matvec=lambda x: x / np.maximum(diagonal, 1e-12))
    try:
        solution, info = cg(
            matrix,
            rhs,
            rtol=parameters.solver_tolerance,
            atol=0.0,
            maxiter=parameters.solver_max_iterations,
            M=preconditioner,
        )
    except TypeError:  # SciPy before the rtol/atol API.
        solution, info = cg(
            matrix,
            rhs,
            tol=parameters.solver_tolerance,
            maxiter=parameters.solver_max_iterations,
            M=preconditioner,
        )
    if info != 0 or not np.all(np.isfinite(solution)):
        raise RuntimeError(f"screened diffusion solver did not converge (info={info})")

    field = np.zeros(flat_domain.shape, dtype=np.float32)
    field[flat_domain] = np.maximum(solution, 0.0).astype(np.float32)
    scale = float(np.max(field[candidate]))
    if not np.isfinite(scale) or scale <= 0:
        raise RuntimeError("screened diffusion field has no positive candidate values")
    field /= scale
    return field


def build_visit_features(
    previous: VisitVolume,
    current: VisitVolume,
    parameters: DiffusionParameters,
) -> VisitFeatures:
    """Build forecast-origin features from the previous and current scans only."""
    spacing = validate_aligned_visits(previous, current)
    previous_images, previous_labels, _, _ = validate_visit(previous, "previous")
    current_images, current_labels, _, _ = validate_visit(current, "current")
    support = brain_support(current_images, current_labels)
    source, candidate, graph_domain = _candidate_region(
        current_labels,
        support,
        spacing,
        parameters.candidate_radius_mm,
        parameters.block_cavity,
    )

    previous_support = brain_support(previous_images, previous_labels)
    previous_source, previous_candidate, _ = _candidate_region(
        previous_labels,
        previous_support,
        spacing,
        parameters.candidate_radius_mm,
        parameters.block_cavity,
    )
    # Count the completed previous transition inside the ROI that was
    # available at its forecast origin. This matches the target definition
    # applied to the current transition and avoids counting label changes from
    # outside that historical candidate region.
    prior_additions = source & ~previous_source & previous_candidate
    prior_count = int(np.count_nonzero(prior_additions))
    sigma_voxels = tuple(float(parameters.prior_smoothing_mm / value) for value in spacing)
    prior_density = gaussian_filter(prior_additions.astype(np.float32), sigma=sigma_voxels)
    peak = float(np.max(prior_density))
    if peak > 0:
        prior_density /= peak

    normalized = {key: _robust_z(current_images[key], support) for key in MODALITIES}
    field = _screened_graph_field(
        normalized,
        current_labels,
        source,
        candidate,
        graph_domain,
        spacing,
        parameters,
    )
    log_field = np.zeros(field.shape, dtype=np.float32)
    log_field[candidate] = np.log(np.maximum(field[candidate], 1e-8))

    cavity = current_labels == 4
    if np.any(cavity):
        distance_to_cavity = distance_transform_edt(~cavity, sampling=spacing).astype(np.float32)
        distance_to_cavity = np.minimum(distance_to_cavity, parameters.candidate_radius_mm)
    else:
        distance_to_cavity = np.full(current_labels.shape, parameters.candidate_radius_mm, dtype=np.float32)

    maps = {
        "distance_to_core_mm": distance_transform_edt(~source, sampling=spacing).astype(np.float32),
        "distance_to_cavity_mm": distance_to_cavity,
        "t1n_z": normalized["t1n"],
        "t1c_z": normalized["t1c"],
        "t2w_z": normalized["t2w"],
        "t2f_z": normalized["t2f"],
        "label2_indicator": (current_labels == 2).astype(np.float32),
        "prior_addition_density": prior_density.astype(np.float32, copy=False),
    }
    return VisitFeatures(
        candidate=candidate,
        source=source,
        recent_additions=prior_additions,
        prior_addition_count=prior_count,
        prior_candidate_voxels=int(np.count_nonzero(previous_candidate)),
        maps=maps,
        pde_field=field,
        log_pde_field=log_field,
    )


def build_origin_features(
    current: VisitVolume,
    parameters: OriginDiffusionParameters,
) -> VisitFeatures:
    """Build next-visit features from one forecast-origin scan only.

    The risk set is the full skull-stripped support in the current MRI. Current
    core and resection-cavity voxels are excluded from the addition risk set;
    no fixed-radius crop or previous-transition feature is used. The spatial
    prior solves on the same support with cavity voxels removed from the graph.
    """
    images, labels, _, spacing = validate_visit(current, "origin")
    support = brain_support(images, labels)
    source = annotated_core_mask(labels)
    cavity = labels == 4
    candidate = support & ~source & ~cavity
    graph_domain = support & (~cavity if parameters.block_cavity else np.ones_like(cavity))
    if not np.any(source):
        raise ValueError("origin visit has no label-1/3 source voxels")
    if not np.any(candidate):
        raise ValueError("origin visit has no full-brain addition risk voxels")

    normalized = {key: _robust_z(images[key], support) for key in MODALITIES}
    field = _screened_graph_field(
        normalized,
        labels,
        source,
        candidate,
        graph_domain,
        spacing,
        parameters,
    )
    log_field = np.zeros(field.shape, dtype=np.float32)
    log_field[candidate] = np.log(np.maximum(field[candidate], 1e-8))
    if np.any(cavity):
        distance_to_cavity = distance_transform_edt(~cavity, sampling=spacing).astype(np.float32)
    else:
        maximum_distance = float(np.linalg.norm(np.asarray(labels.shape) * np.asarray(spacing)))
        distance_to_cavity = np.full(labels.shape, maximum_distance, dtype=np.float32)

    maps = {
        "distance_to_core_mm": distance_transform_edt(~source, sampling=spacing).astype(np.float32),
        "distance_to_cavity_mm": distance_to_cavity,
        "t1n_z": normalized["t1n"],
        "t1c_z": normalized["t1c"],
        "t2w_z": normalized["t2w"],
        "t2f_z": normalized["t2f"],
        "label2_indicator": (labels == 2).astype(np.float32),
        # Retain the common VisitFeatures contract while making history absent.
        "prior_addition_density": np.zeros(labels.shape, dtype=np.float32),
    }
    return VisitFeatures(
        candidate=candidate,
        source=source,
        recent_additions=np.zeros(labels.shape, dtype=bool),
        prior_addition_count=0,
        prior_candidate_voxels=int(np.count_nonzero(candidate)),
        maps=maps,
        pde_field=field,
        log_pde_field=log_field,
    )


def feature_matrix(features: VisitFeatures) -> tuple[np.ndarray, tuple[str, ...]]:
    """Return a compact design matrix for matched non-PDE and residual arms."""
    columns = [np.asarray(features.maps[name][features.candidate], dtype=np.float32) for name in FEATURE_NAMES]
    matrix = np.column_stack(columns).astype(np.float32, copy=False)
    if not np.all(np.isfinite(matrix)):
        raise ValueError("feature matrix contains non-finite values")
    return matrix, FEATURE_NAMES


def origin_feature_matrix(features: VisitFeatures) -> tuple[np.ndarray, tuple[str, ...]]:
    """Return current-scan predictors only; no prior-transition column."""
    columns = [
        np.asarray(features.maps[name][features.candidate], dtype=np.float32)
        for name in ORIGIN_FEATURE_NAMES
    ]
    matrix = np.column_stack(columns).astype(np.float32, copy=False)
    if not np.all(np.isfinite(matrix)):
        raise ValueError("origin feature matrix contains non-finite values")
    return matrix, ORIGIN_FEATURE_NAMES


def fixed_origin_perturbation_draws(
    current: VisitVolume,
    reference: VisitFeatures,
    parameters: OriginDiffusionParameters,
) -> np.ndarray:
    """Recompute the spatial prior after fixed one-voxel core erosions/dilation.

    These two deterministic boundary stresses define sensitivity, not a
    probability distribution over segmentation or registration uncertainty.
    """
    images, labels, _, _ = validate_visit(current, "origin")
    source = annotated_core_mask(labels)
    structure = np.ones((3, 3, 3), dtype=bool)
    eroded = binary_erosion(source, structure=structure, border_value=0)
    dilated = binary_dilation(source, structure=structure, border_value=0)

    def replace_source(new_source: np.ndarray) -> VisitVolume:
        changed = labels.copy()
        changed[source] = 0
        changed[new_source & (changed != 4)] = 1
        return VisitVolume(images=images, labels=changed, affine=current.affine)

    if np.any(eroded):
        eroded_log = build_origin_features(replace_source(eroded), parameters).log_pde_field
    else:
        eroded_log = np.full(reference.log_pde_field.shape, -18.0, dtype=np.float32)
    dilated_log = build_origin_features(replace_source(dilated), parameters).log_pde_field
    draws = np.stack((eroded_log, dilated_log), axis=0).astype(np.float32, copy=False)
    if draws.shape[1:] != reference.candidate.shape:
        raise RuntimeError("origin perturbation draws do not match reference grid")
    return draws


def origin_perturbation_instability(
    reference_log_field: np.ndarray,
    perturbed_log_fields: np.ndarray,
    candidate: np.ndarray,
) -> np.ndarray:
    """Compute median absolute change in the log prior under the fixed stresses."""
    reference = np.asarray(reference_log_field, dtype=np.float32)
    draws = np.asarray(perturbed_log_fields, dtype=np.float32)
    roi = np.asarray(candidate, dtype=bool)
    if reference.shape != roi.shape or draws.ndim != 4 or draws.shape[1:] != roi.shape:
        raise ValueError("origin prior perturbations and candidate mask must share a grid")
    if draws.shape[0] < 2 or not np.all(np.isfinite(draws)):
        raise ValueError("at least two finite origin perturbations are required")
    instability = np.median(np.abs(draws - reference[None, ...]), axis=0).astype(np.float32)
    instability[~roi] = 0.0
    return instability


def origin_stability_gate(
    instability: np.ndarray,
    candidate: np.ndarray,
    *,
    scale: float,
    strength: float = 1.0,
) -> np.ndarray:
    """Shrink the learned residual where the spatial prior is perturbation-sensitive."""
    values = np.asarray(instability, dtype=np.float32)
    roi = np.asarray(candidate, dtype=bool)
    if values.shape != roi.shape or not np.all(np.isfinite(values)):
        raise ValueError("instability and candidate mask must share a finite grid")
    if not np.isfinite(scale) or scale <= 0 or not np.isfinite(strength) or strength < 0:
        raise ValueError("gate scale must be positive and strength nonnegative")
    gate = np.ones(values.shape, dtype=np.float32)
    gate[roi] = np.exp(-strength * values[roi] / scale).astype(np.float32)
    return gate


def make_transition_target(
    current_labels: np.ndarray,
    future_labels: np.ndarray,
    candidate: np.ndarray,
) -> np.ndarray:
    """Mark newly annotated label-1/3 voxels inside a forecast-origin ROI."""
    current = np.asarray(current_labels)
    future = np.asarray(future_labels)
    roi = np.asarray(candidate, dtype=bool)
    if current.shape != future.shape or current.shape != roi.shape:
        raise ValueError("current labels, future labels, and candidate mask must share a grid")
    return (annotated_core_mask(future) & ~annotated_core_mask(current) & roi)


def fixed_perturbation_draws(
    previous: VisitVolume,
    current: VisitVolume,
    reference: VisitFeatures,
    parameters: DiffusionParameters,
) -> tuple[np.ndarray, np.ndarray]:
    """Create deterministic one-voxel stress perturbations for the gate.

    The PDE draws use one-voxel erosion/dilation of the current source mask.
    The history-map draws use one-voxel shifts along each grid axis.  These are
    sensitivity scenarios, not a probabilistic model of segmentation or
    registration error.  The exact scenario set must be frozen before
    evaluation.
    """
    previous_images, previous_labels, _, _ = validate_visit(previous, "previous")
    current_images, current_labels, _, _ = validate_visit(current, "current")
    validate_aligned_visits(previous, current)
    source = annotated_core_mask(current_labels)
    structure = np.ones((3, 3, 3), dtype=bool)
    eroded = binary_erosion(source, structure=structure, border_value=0)
    dilated = binary_dilation(source, structure=structure, border_value=0)

    def replace_source(new_source: np.ndarray) -> VisitVolume:
        labels = current_labels.copy()
        labels[source] = 0
        labels[new_source & (labels != 4)] = 1
        return VisitVolume(images=current_images, labels=labels, affine=current.affine)

    if np.any(eroded):
        erosion_features = build_visit_features(
            VisitVolume(previous_images, previous_labels, previous.affine),
            replace_source(eroded),
            parameters,
        )
        erosion_log = erosion_features.log_pde_field
    else:
        # An erased seed is maximally different from the reference source and
        # should make the gate cautious, rather than aborting the patient.
        erosion_log = np.full(reference.log_pde_field.shape, -18.0, dtype=np.float32)
    dilation_features = build_visit_features(
        VisitVolume(previous_images, previous_labels, previous.affine),
        replace_source(dilated),
        parameters,
    )
    log_fields = np.stack(
        (reference.log_pde_field, erosion_log, dilation_features.log_pde_field),
        axis=0,
    ).astype(np.float32, copy=False)

    prior_maps = [reference.recent_additions.astype(np.float32, copy=False)]
    for axis in range(3):
        for direction in (-1, 1):
            offset = [0.0, 0.0, 0.0]
            offset[axis] = float(direction)
            prior_maps.append(
                shift(reference.recent_additions, shift=offset, order=1, mode="constant", cval=0.0)
            )
    prior_draws = np.stack(prior_maps, axis=0).astype(np.float32, copy=False)
    if log_fields.shape[1:] != reference.candidate.shape or prior_draws.shape[1:] != reference.candidate.shape:
        raise RuntimeError("perturbation draws do not match reference grid")
    return log_fields, prior_draws


def perturbation_stability_gate(
    reference_log_field: np.ndarray,
    perturbed_log_fields: np.ndarray,
    reference_prior_map: np.ndarray,
    perturbed_prior_maps: np.ndarray,
    candidate: np.ndarray,
    log_scale: float,
    prior_scale: float,
    strength: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a heuristic residual gate and local instability score.

    The supplied perturbations must be fixed before evaluation.  The returned
    instability is not a calibrated probability or biological uncertainty.
    """
    if strength < 0 or log_scale <= 0 or prior_scale <= 0:
        raise ValueError("strength must be nonnegative and training scales positive")
    candidate = np.asarray(candidate, dtype=bool)
    base_log = np.asarray(reference_log_field, dtype=np.float32)
    base_prior = np.asarray(reference_prior_map, dtype=np.float32)
    log_draws = np.asarray(perturbed_log_fields, dtype=np.float32)
    prior_draws = np.asarray(perturbed_prior_maps, dtype=np.float32)
    expected_shape = base_log.shape
    if base_prior.shape != expected_shape or candidate.shape != expected_shape:
        raise ValueError("reference maps and candidate mask must share a grid")
    if log_draws.ndim != 4 or prior_draws.ndim != 4:
        raise ValueError("perturbation maps must have shape (draw, x, y, z)")
    if log_draws.shape[1:] != expected_shape or prior_draws.shape[1:] != expected_shape:
        raise ValueError("perturbation maps must share the reference grid")
    if log_draws.shape[0] < 2 or prior_draws.shape[0] < 2:
        raise ValueError("at least two fixed perturbations are required for stability")

    log_instability, prior_instability = perturbation_instability_maps(
        base_log, log_draws, base_prior, prior_draws, candidate
    )

    return stability_gate_from_maps(
        log_instability,
        prior_instability,
        candidate,
        log_scale=log_scale,
        prior_scale=prior_scale,
        strength=strength,
    )


def perturbation_instability_maps(
    reference_log_field: np.ndarray,
    perturbed_log_fields: np.ndarray,
    reference_prior_map: np.ndarray,
    perturbed_prior_maps: np.ndarray,
    candidate: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return raw local median-deviation maps for training-fold scaling."""
    candidate = np.asarray(candidate, dtype=bool)
    base_log = np.asarray(reference_log_field, dtype=np.float32)
    base_prior = np.asarray(reference_prior_map, dtype=np.float32)
    log_draws = np.asarray(perturbed_log_fields, dtype=np.float32)
    prior_draws = np.asarray(perturbed_prior_maps, dtype=np.float32)
    if base_log.shape != base_prior.shape or base_log.shape != candidate.shape:
        raise ValueError("reference maps and candidate mask must share a grid")
    if log_draws.ndim != 4 or prior_draws.ndim != 4:
        raise ValueError("perturbation maps must have shape (draw, x, y, z)")
    if log_draws.shape[1:] != base_log.shape or prior_draws.shape[1:] != base_log.shape:
        raise ValueError("perturbation maps must share the reference grid")
    if log_draws.shape[0] < 2 or prior_draws.shape[0] < 2:
        raise ValueError("at least two fixed perturbations are required for stability")
    log_instability = np.median(np.abs(log_draws - base_log[None, ...]), axis=0)
    prior_instability = np.median(np.abs(prior_draws - base_prior[None, ...]), axis=0)
    log_instability[~candidate] = 0.0
    prior_instability[~candidate] = 0.0
    return log_instability.astype(np.float32), prior_instability.astype(np.float32)


def stability_gate_from_maps(
    log_instability: np.ndarray,
    prior_instability: np.ndarray,
    candidate: np.ndarray,
    *,
    log_scale: float,
    prior_scale: float,
    strength: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Scale raw instability with training-only scales and gate the residual."""
    log_map = np.asarray(log_instability, dtype=np.float32)
    prior_map = np.asarray(prior_instability, dtype=np.float32)
    roi = np.asarray(candidate, dtype=bool)
    if log_map.shape != prior_map.shape or log_map.shape != roi.shape:
        raise ValueError("instability maps and candidate mask must share a grid")
    if strength < 0 or log_scale <= 0 or prior_scale <= 0:
        raise ValueError("strength must be nonnegative and training scales positive")
    if not np.all(np.isfinite(log_map)) or not np.all(np.isfinite(prior_map)):
        raise ValueError("instability maps must be finite")
    instability = np.zeros(roi.shape, dtype=np.float32)
    scaled_log = np.clip(log_map / log_scale, 0.0, 1.0)
    scaled_prior = np.clip(prior_map / prior_scale, 0.0, 1.0)
    instability[roi] = 0.5 * scaled_log[roi] + 0.5 * scaled_prior[roi]
    gate = np.ones(roi.shape, dtype=np.float32)
    gate[roi] = np.exp(-strength * instability[roi])
    return gate, instability
