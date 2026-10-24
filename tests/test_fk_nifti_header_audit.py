from __future__ import annotations

import numpy as np
import nibabel as nib

from src.fk_nifti_header_audit import audit_nifti_headers


def _write_timepoint(root, *, segmentation_affine=None):
    folder = root / "PatientID_test" / "Timepoint_1"
    folder.mkdir(parents=True)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    for name in ("t1n", "t1c", "t2w", "t2f"):
        image = nib.Nifti1Image(np.zeros((4, 5, 6), dtype=np.uint8), affine)
        nib.save(image, folder / f"{name}.nii.gz")
    mask_image = nib.Nifti1Image(
        np.zeros((4, 5, 6), dtype=np.uint8),
        affine if segmentation_affine is None else segmentation_affine,
    )
    nib.save(mask_image, folder / "segmentation.nii.gz")


def test_header_audit_checks_alignment_without_reading_voxel_arrays(tmp_path, monkeypatch):
    _write_timepoint(tmp_path)

    def no_voxel_access(*args, **kwargs):
        raise AssertionError("header audit accessed voxel values")

    monkeypatch.setattr(nib.Nifti1Image, "get_fdata", no_voxel_access)
    report = audit_nifti_headers(tmp_path)

    assert report["status"] == "consistent_by_header_metadata"
    assert report["summary"]["nifti_files_seen"] == 5
    assert report["summary"]["nifti_headers_read"] == 5
    assert report["summary"]["timepoints_without_shape_spacing_or_affine_mismatch"] == 1
    assert "PatientID_test" not in str(report)


def test_header_audit_blocks_mask_affine_mismatch(tmp_path):
    mask_affine = np.diag([1.0, 1.0, 1.0, 1.0])
    mask_affine[0, 3] = 2.0
    _write_timepoint(tmp_path, segmentation_affine=mask_affine)

    report = audit_nifti_headers(tmp_path)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert report["status"] == "blocked"
    assert codes["within_timepoint_affine_mismatches"] == 1
