from __future__ import annotations

import importlib.util
import csv
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "cfb_origin_preflight.py"
SPEC = importlib.util.spec_from_file_location("cfb_origin_preflight", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
PREFLIGHT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PREFLIGHT
SPEC.loader.exec_module(PREFLIGHT)
PATIENT_ID = "deidentified-case-001"
ORIGIN_NAME = f"{PATIENT_ID}_t1_t1gd.nii.gz"


class FakeImage:
    def __init__(self, *, dataobj=None, affine=None, allow_voxels=False, accessed=None):
        self.shape = (12, 12, 12)
        self.affine = np.eye(4) if affine is None else np.asarray(affine, dtype=np.float64)
        self.header = SimpleNamespace(get_xyzt_units=lambda: ("mm", "sec"))
        self._dataobj = dataobj
        self.allow_voxels = allow_voxels
        self.accessed = accessed

    @property
    def dataobj(self):
        if not self.allow_voxels:
            raise AssertionError("preflight read non-origin MRI or mask voxels")
        if self.accessed is not None:
            self.accessed.append("t1_t1gd")
        return self._dataobj


def _row(tmp_path: Path) -> dict[str, str]:
    row = {
        "patient_id": PATIENT_ID,
        "identity_mapping_source": "v4_crosswalk_row_7",
        "provenance_t0": "humanSegmentation",
        "provenance_source_t0": "source-paper-and-mask-ledger-row-7",
        "provenance_t1": "humanSegmentation",
        "provenance_source_t1": "source-paper-and-mask-ledger-row-8",
        "provenance_t2": "humanSegmentation",
        "provenance_source_t2": "source-paper-and-mask-ledger-row-9",
        "review_status_t0": "not_applicable_manual",
        "review_source_t0": "mask-review-ledger-row-7",
        "review_status_t1": "not_applicable_manual",
        "review_source_t1": "mask-review-ledger-row-8",
        "review_status_t2": "not_applicable_manual",
        "review_source_t2": "mask-review-ledger-row-9",
        "t0_to_t1_transform_documented": "true",
        "t0_to_t1_transform_source": "transform-ledger-row-7",
        "t0_to_t1_qc_pass": "true",
        "t0_to_t1_qc_source": "registration-qc-row-7",
        "t2_to_t1_transform_documented": "true",
        "t2_to_t1_transform_source": "transform-ledger-row-9",
        "t2_to_t1_qc_pass": "true",
        "t2_to_t1_qc_source": "registration-qc-row-9",
    }
    for key, role in PREFLIGHT.IMAGE_COLUMNS.items():
        visit, modality = role.split("_")
        visit_index = int(visit[1:])
        path = tmp_path / PATIENT_ID / f"visit_{visit_index:02d}" / f"{PATIENT_ID}_{visit}_{modality}.nii.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        row[key] = str(path.relative_to(tmp_path))
    return row


def test_preflight_reads_only_t1_origin_mri_voxels(tmp_path: Path) -> None:
    row = _row(tmp_path)
    accessed: list[str] = []
    intensity = np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1

    def loader(path: Path) -> FakeImage:
        is_origin = path.name == ORIGIN_NAME
        return FakeImage(
            dataobj=intensity,
            allow_voxels=is_origin,
            accessed=accessed if is_origin else None,
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert result["preflight_pass"] is True
    assert result["scoring_authorized"] is False
    assert result["target_mask_voxels_read"] is False
    assert result["t0_t2_mri_voxels_read"] is False
    assert accessed == ["t1_t1gd"]


def test_preflight_flags_flat_origin_without_reading_masks(tmp_path: Path) -> None:
    row = _row(tmp_path)

    def loader(path: Path) -> FakeImage:
        origin = path.name == ORIGIN_NAME
        return FakeImage(
            dataobj=np.full((12, 12, 12), 2.0, dtype=np.float32),
            allow_voxels=origin,
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert result["preflight_pass"] is False
    assert result["origin_qc_pass"] is False
    assert "origin_no_usable_intensity_variation" in result["errors"]


def test_preflight_rejects_unaligned_header_and_mixed_provenance(tmp_path: Path) -> None:
    row = _row(tmp_path)
    row["provenance_t2"] = "genSegmentation"
    shifted = np.eye(4)
    shifted[0, 3] = 1.0

    def loader(path: Path) -> FakeImage:
        affine = shifted if path.name == f"{PATIENT_ID}_t2_gtv.nii.gz" else np.eye(4)
        return FakeImage(
            dataobj=np.indices((12, 12, 12)).sum(axis=0).astype(np.float32) + 1,
            affine=affine,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert result["headers_pass"] is False
    assert result["provenance_pass"] is False
    assert "t2_gtv_not_aligned_to_t1_reference" in result["errors"]
    assert "provenance_classes_not_matched" in result["errors"]


def test_preflight_rejects_unknown_generated_mask_review_status(tmp_path: Path) -> None:
    row = _row(tmp_path)
    row["provenance_t0"] = row["provenance_t1"] = row["provenance_t2"] = "genSegmentation"
    row["review_status_t0"] = row["review_status_t1"] = row["review_status_t2"] = "unverified"

    def loader(path: Path) -> FakeImage:
        return FakeImage(
            dataobj=np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert result["provenance_pass"] is True
    assert result["mask_review_pass"] is False
    assert "mask_review_status_t0_unknown" in result["errors"]
    assert "mask_review_status_t1_unknown" in result["errors"]
    assert "mask_review_status_t2_unknown" in result["errors"]


def test_preflight_rejects_review_status_incompatible_with_human_mask(tmp_path: Path) -> None:
    row = _row(tmp_path)
    row["review_status_t0"] = row["review_status_t1"] = row["review_status_t2"] = "reviewed"

    def loader(path: Path) -> FakeImage:
        return FakeImage(
            dataobj=np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert result["mask_review_pass"] is False
    assert "mask_review_status_t0_invalid_for_provenance" in result["errors"]


def test_qc_documented_value_does_not_count_as_qc_pass(tmp_path: Path) -> None:
    row = _row(tmp_path)
    row["t0_to_t1_qc_pass"] = "documented"

    def loader(path: Path) -> FakeImage:
        return FakeImage(
            dataobj=np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert result["t0_to_t1_documented_qc_pass"] is False
    assert "t0_to_t1_qc_not_passed" in result["errors"]


def test_swapped_role_path_stops_before_any_image_loader_call(tmp_path: Path) -> None:
    row = _row(tmp_path)
    row["t0_gtv_path"] = row["t0_t1gd_path"]
    calls: list[Path] = []

    result = PREFLIGHT.preflight_rows(
        [row],
        tmp_path,
        image_loader=lambda path: (calls.append(path), FakeImage())[1],
    )[0]

    assert calls == []
    assert "t0_gtv_filename_identity_mismatch" in result["errors"]
    assert result["target_mask_voxels_read"] is False


def test_package_numeric_id_layout_is_accepted(tmp_path: Path) -> None:
    row = _row(tmp_path)
    row["patient_id"] = "5"
    for column, role in PREFLIGHT.IMAGE_COLUMNS.items():
        visit, modality = role.split("_")
        path = tmp_path / "005" / visit / f"5_{visit}_{modality}.nii.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        row[column] = str(path.relative_to(tmp_path))

    values = np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1

    def loader(path: Path) -> FakeImage:
        return FakeImage(
            dataobj=values if path.name == "5_t1_t1gd.nii.gz" else None,
            allow_voxels=path.name == "5_t1_t1gd.nii.gz",
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert not any("_path_" in error or "_directory_mismatch" in error or "_filename_identity" in error for error in result["errors"])
    assert result["headers_pass"] is True
    assert result["origin_qc_pass"] is True
    assert result["metadata_gates_pass"] is True
    assert result["scoring_authorized"] is False


def test_preflight_rejects_missing_crosswalk_files(tmp_path: Path) -> None:
    row = _row(tmp_path)
    (tmp_path / row["t2_gtv_path"]).unlink()
    result = PREFLIGHT.preflight_rows(
        [row],
        tmp_path,
        image_loader=lambda _path: FakeImage(
            dataobj=np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1,
            allow_voxels=True,
        ),
    )[0]

    assert result["preflight_pass"] is False
    assert "t2_gtv_file_not_found" in result["errors"]
    assert result["target_mask_voxels_read"] is False


def test_preflight_rejects_paths_outside_trusted_package_root_before_loading(tmp_path: Path) -> None:
    row = _row(tmp_path)
    package_root = tmp_path / "trusted_package"
    package_root.mkdir()
    calls: list[Path] = []

    result = PREFLIGHT.preflight_rows(
        [row],
        tmp_path,
        package_root=package_root,
        image_loader=lambda path: (calls.append(path), FakeImage())[1],
    )[0]

    assert calls == []
    assert "t1_t1gd_outside_package_root" in result["errors"]
    assert result["target_mask_voxels_read"] is False


def test_preflight_rejects_nonorthogonal_geometry(tmp_path: Path) -> None:
    row = _row(tmp_path)
    sheared = np.eye(4)
    sheared[0, 1] = 0.2

    def loader(path: Path) -> FakeImage:
        affine = sheared if path.name == f"{PATIENT_ID}_t0_gtv.nii.gz" else np.eye(4)
        return FakeImage(
            dataobj=np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1,
            affine=affine,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert result["headers_pass"] is False
    assert "t0_gtv_non_orthogonal_affine" in result["errors"]


def test_native_t0_t2_mri_affines_are_allowed_but_gtv_masks_must_match_t1(tmp_path: Path) -> None:
    row = _row(tmp_path)
    values = np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1
    native_t0 = np.eye(4)
    native_t0[0, 3] = 7.0
    native_t2 = np.eye(4)
    native_t2[1, 3] = 9.0

    def loader(path: Path) -> FakeImage:
        affine = {
            f"{PATIENT_ID}_t0_t1gd.nii.gz": native_t0,
            f"{PATIENT_ID}_t2_t1gd.nii.gz": native_t2,
        }.get(path.name, np.eye(4))
        return FakeImage(
            dataobj=values,
            affine=affine,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    native_result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=loader)[0]

    assert native_result["headers_pass"] is True
    assert native_result["preflight_pass"] is True
    assert native_result["file_headers"]["t0_t1gd"]["aligned_to_t1_reference"] is None
    assert native_result["file_headers"]["t0_gtv"]["aligned_to_t1_reference"] is True

    def misaligned_mask_loader(path: Path) -> FakeImage:
        if path.name == f"{PATIENT_ID}_t2_gtv.nii.gz":
            return FakeImage(dataobj=values, affine=native_t2)
        return loader(path)

    mask_result = PREFLIGHT.preflight_rows([row], tmp_path, image_loader=misaligned_mask_loader)[0]

    assert mask_result["headers_pass"] is False
    assert "t2_gtv_not_aligned_to_t1_reference" in mask_result["errors"]


def test_run_report_hashes_crosswalk_script_and_protocol(tmp_path: Path, monkeypatch) -> None:
    row = _row(tmp_path)
    crosswalk = tmp_path / "candidate_crosswalk.csv"
    with crosswalk.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREFLIGHT.REQUIRED_COLUMNS)
        writer.writeheader()
        writer.writerow(row)

    protocol = tmp_path / "protocol.md"
    protocol.write_text("Frozen synthetic protocol\n", encoding="utf-8")
    output_dir = tmp_path / "preflight_report"
    values = np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1

    def image_loader(path: Path) -> FakeImage:
        return FakeImage(
            dataobj=values if path.name == ORIGIN_NAME else None,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    monkeypatch.setattr(PREFLIGHT, "_load_nifti", image_loader)
    summary = PREFLIGHT.run_preflight(
        crosswalk,
        output_dir,
        package_root=tmp_path,
        protocol_path=protocol,
    )
    saved = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))

    assert summary["preflight_script_sha256"] == hashlib.sha256(MODULE_PATH.read_bytes()).hexdigest()
    assert summary["analysis_protocol_sha256"] == hashlib.sha256(protocol.read_bytes()).hexdigest()
    assert saved["source_crosswalk_sha256"] == hashlib.sha256(crosswalk.read_bytes()).hexdigest()
    assert saved["candidate_rows"] == 1
    assert saved["scoring_authorized"] is False


def test_run_preflight_uses_default_loader_when_not_injected(tmp_path: Path, monkeypatch) -> None:
    row = _row(tmp_path)
    crosswalk = tmp_path / "candidate_crosswalk.csv"
    with crosswalk.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREFLIGHT.REQUIRED_COLUMNS)
        writer.writeheader()
        writer.writerow(row)

    values = np.arange(12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12) + 1
    calls: list[Path] = []

    def loader(path: Path) -> FakeImage:
        calls.append(path)
        return FakeImage(
            dataobj=values if path.name == ORIGIN_NAME else None,
            allow_voxels=path.name == ORIGIN_NAME,
        )

    monkeypatch.setattr(PREFLIGHT, "_load_nifti", loader)
    protocol = tmp_path / "protocol.md"
    protocol.write_text("Frozen synthetic protocol\n", encoding="utf-8")
    report = PREFLIGHT.run_preflight(
        crosswalk,
        tmp_path / "default_loader_report",
        package_root=tmp_path,
        protocol_path=protocol,
    )

    assert len(calls) == len(PREFLIGHT.IMAGE_COLUMNS)
    assert report["candidate_rows"] == 1
    assert report["preflight_pass_rows"] == 1
