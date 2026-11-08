from __future__ import annotations

import csv
import io
import json
import zipfile

import pytest

from src.fk_duke_landmark_inventory import build_landmark_inventory


def _csv(points):
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(["x", "y", "z"])
    writer.writerows(points)
    return stream.getvalue()


def _archive(path, *, preoperative_near=((1, 2, 3), (7, 8, 9)), recurrence_near=((4, 5, 6), (8, 9, 10))):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        visits = ((1, "preoperative"), (4, "recurrence"))
        for scan, visit in visits:
            for zone in ("near", "far"):
                points = ((1, 2, 3), (7, 8, 9)) if zone == "near" else ((10, 11, 12),)
                if visit == "preoperative" and zone == "near":
                    points = preoperative_near
                if visit == "recurrence" and zone == "near":
                    points = recurrence_near
                archive.writestr(
                    f"DATASET_FINAL/Patient1/Patient1_scan{scan}_landmarks_{zone}.csv",
                    _csv(points),
                )
            archive.writestr(
                f"DATASET_FINAL/Patient1/scan{scan}_{visit}/Patient1_scan{scan}_t1c.nii.gz",
                b"not an image",
            )
        archive.writestr(
            "DATASET_FINAL/Patient1/scan4_recurrence/Patient1_scan4_tumor-seg.nii.gz", b"not a mask"
        )


def test_row_order_diagnostic_is_reproducible_and_lower_than_shuffled_pairs(tmp_path):
    path = tmp_path / "paired_landmarks.zip"
    _archive(
        path,
        preoperative_near=((1, 2, 3), (7, 8, 9), (1, 8, 3), (7, 2, 9)),
        recurrence_near=((2, 2, 3), (8, 8, 9), (2, 8, 3), (8, 2, 9)),
    )

    first = build_landmark_inventory(path, volume_shape=(20, 20, 20), permutation_replicates=200, seed=17)
    second = build_landmark_inventory(path, volume_shape=(20, 20, 20), permutation_replicates=200, seed=17)
    diagnostic = first["row_order_pairing_diagnostic"]

    assert diagnostic == second["row_order_pairing_diagnostic"]
    assert diagnostic["patient_zone_pair_count"] == 1
    assert diagnostic["median_row_order_distance_voxels"] < diagnostic["median_shuffled_distance_voxels"]
    assert diagnostic["fraction_row_order_median_below_shuffled_median"] == 1.0
    assert diagnostic["permutation_replicates_per_pair"] == 200
    assert diagnostic["seed"] == 17
    assert diagnostic["identity_permutations_excluded"] is True


def test_row_order_diagnostic_excludes_identity_shuffle_for_two_points(tmp_path):
    path = tmp_path / "two_point_landmarks.zip"
    _archive(path, recurrence_near=((2, 2, 3), (8, 8, 9)))

    diagnostic = build_landmark_inventory(
        path, volume_shape=(20, 20, 20), permutation_replicates=100, seed=5
    )["row_order_pairing_diagnostic"]

    assert diagnostic["patient_zone_pair_count"] == 1
    assert diagnostic["median_row_order_distance_voxels"] < diagnostic["median_shuffled_distance_voxels"]


def test_row_order_diagnostic_skips_pairs_with_invalid_coordinates(tmp_path):
    path = tmp_path / "invalid_pair.zip"
    _archive(path, recurrence_near=((4, 5, 6), ("bad", 9, 10)))

    result = build_landmark_inventory(path, volume_shape=(20, 20, 20), permutation_replicates=20)

    diagnostic = result["row_order_pairing_diagnostic"]
    assert diagnostic["patient_zone_pair_count"] == 0
    assert diagnostic["median_row_order_distance_voxels"] is None
    assert result["near_pair_point_count_mismatches"] == 1


def test_row_order_diagnostic_skips_equal_count_pairs_with_invalid_rows(tmp_path):
    path = tmp_path / "invalid_equal_count_pair.zip"
    _archive(
        path,
        preoperative_near=((1, 2, 3), (5, 6, 7), ("bad", 0, 0)),
        recurrence_near=((2, 2, 3), (6, 6, 7), ("nan", 0, 0)),
    )

    result = build_landmark_inventory(path, volume_shape=(20, 20, 20), permutation_replicates=20)

    diagnostic = result["row_order_pairing_diagnostic"]
    assert result["near_pair_point_count_mismatches"] == 0
    assert diagnostic["patient_zone_pair_count"] == 0
    assert diagnostic["pairs_skipped_coordinate_issues"] == 1


@pytest.mark.parametrize("replicates", [0, -1, True, 1.5])
def test_row_order_diagnostic_rejects_invalid_permutation_replicate_counts(tmp_path, replicates):
    path = tmp_path / "invalid_replicates.zip"
    _archive(path)

    with pytest.raises(ValueError, match="permutation_replicates"):
        build_landmark_inventory(path, permutation_replicates=replicates)


@pytest.mark.parametrize("seed", [-1, True, 1.5])
def test_row_order_diagnostic_rejects_invalid_seed(tmp_path, seed):
    path = tmp_path / "invalid_seed.zip"
    _archive(path)

    with pytest.raises(ValueError, match="seed"):
        build_landmark_inventory(path, seed=seed)


def test_landmark_inventory_checks_pairs_and_never_opens_image_or_segmentation(tmp_path, monkeypatch):
    path = tmp_path / "landmarks.zip"
    _archive(path)
    opened = []
    original_open = zipfile.ZipFile.open

    def tracked_open(self, name, *args, **kwargs):
        member = name.filename if isinstance(name, zipfile.ZipInfo) else str(name)
        opened.append(member)
        if member.casefold().endswith((".nii.gz", "tumor-seg.nii.gz")):
            raise AssertionError("landmark inventory must not open image or segmentation members")
        return original_open(self, name, *args, **kwargs)

    monkeypatch.setattr(zipfile.ZipFile, "open", tracked_open)
    result = build_landmark_inventory(path, volume_shape=(20, 20, 20))

    assert result["archive_landmark_csv_count"] == 4
    assert result["patient_count_with_any_parsed_landmarks"] == 1
    assert result["patients_with_all_four_visit_zone_files"] == 1
    assert result["missing_patient_visit_zone_groups"] == 0
    assert result["near_pair_point_count_mismatches"] == 0
    assert result["far_pair_point_count_mismatches"] == 0
    assert result["coordinate_quality_totals"] == {
        "invalid_coordinate_rows": 0,
        "fractional_coordinate_rows": 0,
        "out_of_bounds_rows": 0,
        "duplicate_points_within_file": 0,
    }
    assert result["read_scope"]["nifti_member_bodies_opened"] is False
    assert "Patient1" not in json.dumps(result)
    assert not any(name.casefold().endswith(".nii.gz") for name in opened)


def test_landmark_inventory_reports_pair_count_mismatch(tmp_path):
    path = tmp_path / "mismatch.zip"
    _archive(path, recurrence_near=((4, 5, 6),))
    result = build_landmark_inventory(path, volume_shape=(20, 20, 20))
    assert result["near_pair_point_count_mismatches"] == 1
    assert result["far_pair_point_count_mismatches"] == 0


def test_landmark_inventory_reports_invalid_noninteger_out_of_bounds_and_duplicate_points(tmp_path):
    path = tmp_path / "quality.zip"
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        archive.writestr(
            "DATASET_FINAL/Patient1/Patient1_scan1_landmarks_near.csv",
            _csv(((1, 2, 3), (1, 2, 3), (20, 0, 0), (1.5, 2, 3), ("nan", 1, 1), ("bad", 1, 1))),
        )
        archive.writestr(
            "DATASET_FINAL/Patient1/scan1_preoperative/Patient1_scan1_t1c.nii.gz", b"not an image"
        )
    result = build_landmark_inventory(path, volume_shape=(20, 20, 20))
    assert result["coordinate_quality_totals"] == {
        "invalid_coordinate_rows": 2,
        "fractional_coordinate_rows": 1,
        "out_of_bounds_rows": 1,
        "duplicate_points_within_file": 1,
    }


def test_landmark_inventory_requires_xyz_columns(tmp_path):
    path = tmp_path / "bad_header.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(
            "DATASET_FINAL/Patient1/Patient1_scan1_landmarks_near.csv",
            "i,j,k\n1,2,3\n",
        )
        archive.writestr(
            "DATASET_FINAL/Patient1/scan1_preoperative/Patient1_scan1_t1c.nii.gz", b"not an image"
        )
    with pytest.raises(ValueError, match="x, y, and z"):
        build_landmark_inventory(path)
