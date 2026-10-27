from __future__ import annotations

import csv
import importlib.util
import sys
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "166_trackb_ucsd_external_test.py"
SPEC = importlib.util.spec_from_file_location("trackb_ucsd_locked_external_test", MODULE_PATH)
RUNNER = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


def _sample(patient_key: str) -> RUNNER.PatientSample:
    return RUNNER.PatientSample(
        patient_key=patient_key,
        x=np.zeros((2, len(RUNNER.FEATURES.ORIGIN_FEATURE_NAMES)), dtype=np.float32),
        prior=np.zeros(2, dtype=np.float32),
        instability=np.zeros(2, dtype=np.float32),
        target=np.asarray([0, 1], dtype=np.uint8),
        risk_voxels=2,
        additions_total=1,
        additions_outside_risk=0,
        additions_outside_current_support=0,
        interval_days=30.0,
    )


def test_visit_file_matching_does_not_confuse_timepoint_one_with_ten(tmp_path: Path) -> None:
    visit = RUNNER.ClinicalVisit("case", 1, "1", "BraTS-case-1", None)
    paths = []
    for visit_number in (1, 10):
        for suffix in (
            "T1pre",
            "T1post",
            "T2",
            "FLAIR",
            "BraTS_tumor_seg",
            "enhancing_cellular_tumor_seg",
            "non_enhancing_cellular_tumor_seg",
            "total_cellular_tumor_seg",
        ):
            paths.append(tmp_path / f"case_{visit_number}_BraTS-case-{visit_number}_{suffix}.nii.gz")

    bundle = RUNNER._files_for_visit(paths, visit)

    assert bundle is not None
    assert len(bundle.paths) == 5
    assert all("case_1_" in path.name and "case_10_" not in path.name for path in bundle.paths.values())
    assert bundle.paths["t2w"].name.endswith("_T2.nii.gz")
    assert bundle.paths["seg"].name.endswith("_BraTS_tumor_seg.nii.gz")


def test_ucsd_standard_brats_segmentation_excludes_cellular_masks() -> None:
    assert RUNNER._modality(Path("case_T2.nii.gz")) == "t2w"
    assert RUNNER._modality(Path("case_BraTS_tumor_seg.nii.gz")) == "seg"
    assert RUNNER._modality(Path("case_enhancing_cellular_tumor_seg.nii.gz")) is None
    assert RUNNER._modality(Path("case_non_enhancing_cellular_tumor_seg.nii.gz")) is None
    assert RUNNER._modality(Path("case_total_cellular_tumor_seg.nii.gz")) is None


def test_preflight_catches_flat_origin_modality_without_future_masks(monkeypatch) -> None:
    shape = (12, 12, 12)
    base = np.indices(shape).sum(axis=0).astype(np.float32) + 1.0
    images = {
        name: base * float(index + 1)
        for index, name in enumerate(RUNNER.FEATURES.MODALITIES)
    }
    images["t2w"] = np.zeros(shape, dtype=np.float32)
    pair = SimpleNamespace(first_files=object(), second_files=object())
    loaded = []
    monkeypatch.setattr(RUNNER, "select_ucsd_pairs", lambda *_args: ([pair], {"gbm": 1}))
    monkeypatch.setattr(
        RUNNER,
        "_load_ucsd_origin_images",
        lambda scan: (loaded.append(scan), images)[1],
    )

    pairs, counts, modality_failures = RUNNER._preflight(
        Namespace(ucsd_clinical=Path("unused"), ucsd_root=Path("unused"))
    )

    assert pairs == [pair]
    assert loaded == [pair.first_files]
    assert counts["origin_inputs_quality_checked"] == 1
    assert counts["origin_inputs_quality_failed"] == 1
    assert modality_failures == {"t2w": 1}


def test_segmentation_loader_rejects_fractional_or_unknown_codes(tmp_path: Path) -> None:
    affine = np.eye(4)
    valid = tmp_path / "valid.nii.gz"
    nib.save(nib.Nifti1Image(np.asarray([[[0, 1, 2, 3, 4]]], dtype=np.int16), affine), valid)
    assert RUNNER._load_labels(valid).tolist() == [[[0, 1, 2, 3, 4]]]

    fractional = tmp_path / "fractional.nii.gz"
    nib.save(nib.Nifti1Image(np.asarray([[[0.0, 1.2]]], dtype=np.float32), affine), fractional)
    with pytest.raises(ValueError, match="noninteger"):
        RUNNER._load_labels(fractional)

    unknown = tmp_path / "unknown.nii.gz"
    nib.save(nib.Nifti1Image(np.asarray([[[0, 5]]], dtype=np.int16), affine), unknown)
    with pytest.raises(ValueError, match="outside the documented"):
        RUNNER._load_labels(unknown)


def test_locked_scoring_rejects_noncanonical_source_and_dirty_protocol(monkeypatch) -> None:
    args = Namespace(
        ucsd_root=RUNNER.DEFAULT_UCSD_ROOT,
        ucsd_clinical=RUNNER.DEFAULT_UCSD_CLINICAL,
        mu_data_root=RUNNER.DEFAULT_MU_ROOT,
        patient_manifest=RUNNER.DEFAULT_PATIENT_MANIFEST,
        pair_manifest=RUNNER.DEFAULT_PAIR_MANIFEST,
        protocol=RUNNER.DEFAULT_PROTOCOL,
    )
    monkeypatch.setattr(RUNNER, "_sha256", lambda _path: RUNNER.EXPECTED_UCSD_CLINICAL_SHA256)
    RUNNER._validate_locked_source_paths(args)

    args.pair_manifest = args.pair_manifest.parent / "unfrozen_pairs.csv"
    with pytest.raises(ValueError, match="source roots and protocol"):
        RUNNER._validate_locked_source_paths(args)

    monkeypatch.setattr(
        RUNNER.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(stdout=" M src/166_trackb_ucsd_external_test.py\n"),
    )
    with pytest.raises(RuntimeError, match="uncommitted"):
        RUNNER._assert_locked_inputs_committed(RUNNER.DEFAULT_PROTOCOL)


def test_preflight_affine_gate_rejects_shear() -> None:
    orthogonal = np.eye(4)
    shear = np.eye(4)
    shear[0, 1] = 0.1
    assert RUNNER._valid_orthogonal_affine(orthogonal)
    assert not RUNNER._valid_orthogonal_affine(shear)

    non_one_mm = np.diag([2.0, 1.0, 1.0, 1.0])
    assert not RUNNER._valid_orthogonal_affine(non_one_mm)

    image = nib.Nifti1Image(np.zeros((2, 2, 2), dtype=np.uint8), orthogonal)
    assert not RUNNER._image_grid_valid(image)
    image.header.set_xyzt_units("mm")
    assert RUNNER._image_grid_valid(image)

def test_added_value_requires_lower_confidence_bound_above_margin() -> None:
    assert not RUNNER._added_value_criterion(0.03, 0.01)
    assert not RUNNER._added_value_criterion(0.01, -0.01)
    assert RUNNER._added_value_criterion(0.04, 0.025)


def test_source_content_manifest_detects_same_size_edits(tmp_path: Path) -> None:
    source = tmp_path / "source.nii.gz"
    source.write_bytes(b"AAAA")
    before = RUNNER._hash_files([source], tmp_path)
    size_before = RUNNER._hash_path_size_manifest([source], tmp_path)

    source.write_bytes(b"BBBB")
    after = RUNNER._hash_files([source], tmp_path)
    size_after = RUNNER._hash_path_size_manifest([source], tmp_path)

    assert before != after
    assert size_before == size_after


def test_mu_manifest_alias_loads_only_mu_development_patients(tmp_path: Path, monkeypatch) -> None:
    data_root = tmp_path / "mu"
    patients_path = tmp_path / "patients.csv"
    pairs_path = tmp_path / "pairs.csv"
    patients = [(f"SYN-{index}", "MU" if index < 3 else "MU-Glioma-Post") for index in range(5)]
    patients.append(("SYN-NONMU", "CGGA"))
    with patients_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["patient_id", "dataset", "eligible", "is_gbm"])
        writer.writerows((patient_id, dataset, "True", "True") for patient_id, dataset in patients)

    with pairs_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUNNER.BASE.PAIR_MANIFEST_COLUMNS)
        writer.writeheader()
        for patient_id, _dataset in patients:
            writer.writerow({
                "patient_id": patient_id,
                "tp_in": 1,
                "tp_out": 2,
                "day_in": 0,
                "day_out": 30,
                "dt_days": 30,
                "in_primary": "True",
                "input_core_nonempty": "True",
            })

    for patient_id, _dataset in patients:
        for timepoint in (1, 2):
            directory = data_root / patient_id / f"Timepoint_{timepoint}"
            directory.mkdir(parents=True)
            (directory / f"{patient_id}_tumorMask.nii.gz").touch()
            if timepoint == 1:
                for modality in RUNNER.FEATURES.MODALITIES:
                    (directory / f"{patient_id}_brain_{modality}.nii.gz").touch()

    monkeypatch.setattr(RUNNER.BASE, "_load_visit", lambda *_args: object())
    monkeypatch.setattr(RUNNER.BASE, "_load_mask", lambda *_args: (np.zeros((1, 1, 1), dtype=np.int16), np.eye(4)))
    monkeypatch.setattr(RUNNER, "_patient_sample", lambda patient_key, *_args: _sample(patient_key))
    monkeypatch.setattr(RUNNER, "_validate_one_mm_file_headers", lambda *_args: None)

    samples, source_paths = RUNNER.load_mu_development_samples(
        data_root,
        patients_path,
        pairs_path,
        RUNNER.FEATURES.OriginDiffusionParameters(),
    )

    assert len(samples) == 5
    assert {sample.patient_key for sample in samples} == {patient_id for patient_id, _ in patients[:5]}
    assert len(source_paths) == 30


def test_origin_features_use_only_current_scan_and_whole_support() -> None:
    shape = (18, 18, 18)
    i, j, k = np.indices(shape, dtype=np.float32)
    support = (i >= 2) & (i < 16) & (j >= 2) & (j < 16) & (k >= 2) & (k < 16)
    images = {
        "t1n": np.where(support, i + 0.2 * j + 0.1 * k, 0.0).astype(np.float32),
        "t1c": np.where(support, 0.1 * i + j + 0.3 * k, 0.0).astype(np.float32),
        "t2w": np.where(support, 0.4 * i + 0.1 * j + k, 0.0).astype(np.float32),
        "t2f": np.where(support, i + j + k, 0.0).astype(np.float32),
    }
    labels = np.zeros(shape, dtype=np.int16)
    labels[8:10, 8:10, 8:10] = 1
    labels[10:12, 8:10, 8:10] = 2
    labels[5, 5, 5] = 4
    current = RUNNER.FEATURES.VisitVolume(images=images, labels=labels, affine=np.eye(4))
    parameters = RUNNER.FEATURES.OriginDiffusionParameters(solver_tolerance=1e-5)

    visit = RUNNER.FEATURES.build_origin_features(current, parameters)
    matrix, names = RUNNER.FEATURES.origin_feature_matrix(visit)
    draws = RUNNER.FEATURES.fixed_origin_perturbation_draws(current, visit, parameters)
    instability = RUNNER.FEATURES.origin_perturbation_instability(visit.log_pde_field, draws, visit.candidate)
    gate = RUNNER.FEATURES.origin_stability_gate(instability, visit.candidate, scale=1.0)

    assert names == RUNNER.FEATURES.ORIGIN_FEATURE_NAMES
    assert len(names) == 7
    assert matrix.shape == (int(visit.candidate.sum()), 7)
    assert not np.any(visit.candidate & (visit.maps["prior_addition_density"] != 0))
    assert np.any(visit.candidate & (visit.maps["label2_indicator"] == 1))
    assert np.all(np.isfinite(draws)) and np.all(np.isfinite(instability))
    assert np.all((gate[visit.candidate] >= 0) & (gate[visit.candidate] <= 1))
