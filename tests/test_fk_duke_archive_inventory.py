from __future__ import annotations

import hashlib
import zipfile

import pytest

from src import fk_duke_archive_inventory as inventory


def make_zip(path, members):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, content in members:
            zf.writestr(name, content)
    return path


def test_inventory_lists_only_central_directory_metadata(tmp_path, monkeypatch):
    archive = make_zip(
        tmp_path / "synthetic.zip",
        [
            ("Patient01/Scan1/T1c.nii.gz", b"synthetic image payload"),
            ("Patient01/Scan1/segmentation.nii.gz", b"synthetic mask payload"),
            ("patient_metadata.csv", b"Patient,IDH\nPatient01,WT\n"),
        ],
    )

    def no_member_open(*_args, **_kwargs):
        raise AssertionError("inventory must not open or decompress ZIP members")

    monkeypatch.setattr(zipfile.ZipFile, "open", no_member_open)
    record = inventory.inventory_zip_members(archive)

    assert record["member_count"] == 3
    assert record["read_scope"] == "ZIP central directory only; no members opened or decompressed"
    assert record["csv_member_names"] == ["patient_metadata.csv"]
    assert {item["name"] for item in record["members"]} == {
        "Patient01/Scan1/T1c.nii.gz",
        "Patient01/Scan1/segmentation.nii.gz",
        "patient_metadata.csv",
    }


def test_archive_integrity_checks_expected_size_and_md5(tmp_path):
    archive = make_zip(tmp_path / "DATASET_FINAL.zip", [("metadata.csv", b"id,scan\nP1,1\n")])
    payload = archive.read_bytes()
    expected_md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()

    record = inventory.verify_archive_integrity(
        archive, expected_size_bytes=len(payload), expected_md5=expected_md5
    )

    assert record["size_bytes"] == len(payload)
    assert record["md5"] == expected_md5
    assert len(record["sha256"]) == 64


def test_inspect_csv_header_reads_only_named_metadata_header(tmp_path, monkeypatch):
    archive = make_zip(
        tmp_path / "synthetic.zip",
        [
            ("patient_metadata.csv", b"Patient ID,Scan Date,IDH status\nP1,0,WT\n"),
            ("Patient01/Scan1/landmarks.csv", b"x,y,z\n1,2,3\n"),
        ],
    )
    original_open = zipfile.ZipFile.open
    opened = []

    def track_open(self, name, mode="r", pwd=None, force_zip64=False):
        opened.append(name.filename if isinstance(name, zipfile.ZipInfo) else name)
        return original_open(self, name, mode=mode, pwd=pwd, force_zip64=force_zip64)

    monkeypatch.setattr(zipfile.ZipFile, "open", track_open)
    record = inventory.inspect_csv_header(archive, "patient_metadata.csv")
    assert record["column_names"] == ["Patient ID", "Scan Date", "IDH status"]
    assert opened == ["patient_metadata.csv"]
    with pytest.raises(ValueError, match="landmark/coordinate CSVs"):
        inventory.inspect_csv_header(archive, "Patient01/Scan1/landmarks.csv")
    with pytest.raises(ValueError, match="not found"):
        inventory.inspect_csv_header(archive, "missing.csv")


def test_archive_integrity_rejects_partial_size_before_hashing(tmp_path, monkeypatch):
    archive = make_zip(tmp_path / "partial.zip", [("x.csv", b"small")])

    def no_hash(*_args, **_kwargs):
        raise AssertionError("wrong-size partial file must be rejected before hashing")

    monkeypatch.setattr(inventory, "_hash_file", no_hash)
    with pytest.raises(ValueError, match="archive size mismatch"):
        inventory.verify_archive_integrity(archive, expected_size_bytes=archive.stat().st_size + 1)


def test_archive_integrity_rejects_checksum_mismatch(tmp_path):
    archive = make_zip(tmp_path / "wrong-checksum.zip", [("metadata.csv", b"id\nP1\n")])
    with pytest.raises(ValueError, match="archive MD5 mismatch"):
        inventory.verify_archive_integrity(
            archive, expected_size_bytes=archive.stat().st_size, expected_md5="0" * 32
        )


def test_inventory_rejects_path_traversal_and_duplicate_member_names(tmp_path):
    traversal = make_zip(tmp_path / "traversal.zip", [("../outside.csv", b"x")])
    with pytest.raises(ValueError, match="unsafe ZIP member path"):
        inventory.inventory_zip_members(traversal)

    duplicate = tmp_path / "duplicate.zip"
    with zipfile.ZipFile(duplicate, "w") as zf:
        with pytest.warns(UserWarning, match="Duplicate name"):
            zf.writestr("metadata.csv", b"first")
            zf.writestr("metadata.csv", b"second")
    with pytest.raises(ValueError, match="duplicate ZIP member name"):
        inventory.inventory_zip_members(duplicate)


def test_inventory_rejects_non_zip_and_bad_expected_md5(tmp_path):
    not_zip = tmp_path / "not.zip"
    not_zip.write_bytes(b"not a zip")
    with pytest.raises(ValueError, match="invalid ZIP archive"):
        inventory.inventory_zip_members(not_zip)

    with pytest.raises(ValueError, match="expected_md5 must be 32 hexadecimal characters"):
        inventory.verify_archive_integrity(not_zip, expected_md5="not-a-checksum")
