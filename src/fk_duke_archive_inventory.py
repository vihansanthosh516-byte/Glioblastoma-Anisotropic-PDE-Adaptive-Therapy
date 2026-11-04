"""Metadata-only inventory for the restricted Duke GBM archive.

The integrity check reads compressed archive bytes. The ZIP inventory reads only the central directory and
never opens or decompresses a member. The separate CSV-header helper opens only one explicitly named,
non-landmark CSV and returns its first row. Neither function reads MRI or mask members.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path as _Path, PurePosixPath, PureWindowsPath
import zipfile

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
EXPECTED_ARCHIVE_NAME = "DATASET_FINAL.zip"
EXPECTED_ARCHIVE_SIZE_BYTES = 17_043_569_844
EXPECTED_ARCHIVE_MD5 = "99580fb8ccb3daed186dc6afe5b5bf2f"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "track_b_duke_inventory_20261103" / "archive_inventory.json"


def _hash_file(path, *, chunk_bytes: int = 8 * 1024 * 1024) -> dict[str, str]:
    path = _Path(path)
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_bytes):
            md5.update(chunk)
            sha256.update(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def verify_archive_integrity(path, *, expected_size_bytes: int = EXPECTED_ARCHIVE_SIZE_BYTES,
                             expected_md5: str = EXPECTED_ARCHIVE_MD5) -> dict:
    """Check completed archive size and Zenodo MD5, then record local SHA-256."""
    path = _Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"archive file does not exist: {path}")
    if (isinstance(expected_size_bytes, bool) or not isinstance(expected_size_bytes, int)
            or expected_size_bytes < 1):
        raise ValueError("expected_size_bytes must be a positive integer")
    if (not isinstance(expected_md5, str) or len(expected_md5) != 32
            or any(char not in "0123456789abcdefABCDEF" for char in expected_md5)):
        raise ValueError("expected_md5 must be 32 hexadecimal characters")

    actual_size = path.stat().st_size
    if actual_size != expected_size_bytes:
        raise ValueError(
            f"archive size mismatch: expected {expected_size_bytes} bytes, found {actual_size}"
        )
    hashes = _hash_file(path)
    if hashes["md5"].lower() != expected_md5.lower():
        raise ValueError(f"archive MD5 mismatch: expected {expected_md5.lower()}, found {hashes['md5']}")
    return {
        "archive_name": path.name,
        "size_bytes": actual_size,
        "md5": hashes["md5"],
        "sha256": hashes["sha256"],
        "integrity_status": "passed",
        "hash_scope": "compressed archive bytes; no ZIP members opened",
    }


def _validate_member_name(name: str) -> None:
    """Reject ZIP names that would be unsafe or ambiguous on Windows."""
    posix = PurePosixPath(name)
    windows = PureWindowsPath(name)
    if (not name or "\x00" in name or "\\" in name or posix.is_absolute()
            or windows.is_absolute() or windows.drive
            or any(part in {".", ".."} for part in name.split("/"))):
        raise ValueError(f"unsafe ZIP member path: {name!r}")


def inventory_zip_members(path) -> dict:
    """List ZIP central-directory entries and CSV filenames without opening members."""
    path = _Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"archive file does not exist: {path}")
    try:
        with zipfile.ZipFile(path, "r") as archive:
            entries = archive.infolist()
    except (OSError, zipfile.BadZipFile, zipfile.LargeZipFile) as exc:
        raise ValueError(f"invalid ZIP archive: {path.name}") from exc

    seen_names = set()
    members = []
    csv_names = []
    total_uncompressed = 0
    for entry in entries:
        name = entry.filename
        _validate_member_name(name)
        folded_name = name.casefold()
        if folded_name in seen_names:
            raise ValueError(f"duplicate ZIP member name: {name!r}")
        seen_names.add(folded_name)
        is_directory = entry.is_dir()
        lower_name = name.lower()
        is_csv = not is_directory and lower_name.endswith(".csv")
        row = {
            "name": name,
            "is_directory": is_directory,
            "compressed_bytes": int(entry.compress_size),
            "uncompressed_bytes": int(entry.file_size),
            "crc32": f"{entry.CRC:08x}",
            "is_csv_candidate": is_csv,
            "is_nifti_candidate": (not is_directory and
                                   (lower_name.endswith(".nii") or lower_name.endswith(".nii.gz"))),
        }
        members.append(row)
        if is_csv:
            csv_names.append(name)
        if not is_directory:
            total_uncompressed += int(entry.file_size)

    return {
        "archive_name": path.name,
        "archive_size_bytes": path.stat().st_size,
        "member_count": len(members),
        "file_count": sum(not item["is_directory"] for item in members),
        "csv_member_names": sorted(csv_names, key=str.casefold),
        "nifti_member_count": sum(item["is_nifti_candidate"] for item in members),
        "declared_uncompressed_bytes": total_uncompressed,
        "read_scope": "ZIP central directory only; no members opened or decompressed",
        "members": members,
    }


def inspect_csv_header(path, member_name: str) -> dict:
    """Read column names from one explicitly named CSV member; refuse landmark/coordinate tables."""
    if not isinstance(member_name, str) or not member_name.lower().endswith(".csv"):
        raise ValueError("member_name must name one CSV file")
    lower_name = member_name.casefold()
    if any(term in lower_name for term in ("landmark", "coordinate", "coordinates", "point")):
        raise ValueError("landmark/coordinate CSVs are outside the metadata-header scope")
    _validate_member_name(member_name)
    try:
        with zipfile.ZipFile(path, "r") as archive:
            info = archive.getinfo(member_name)
            if info.is_dir():
                raise ValueError("member_name must name a CSV file, not a directory")
            with archive.open(info, "r") as binary:
                with io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as text:
                    header = next(csv.reader(text), None)
    except KeyError as exc:
        raise ValueError(f"CSV member not found in archive: {member_name}") from exc
    except (OSError, csv.Error, zipfile.BadZipFile, zipfile.LargeZipFile) as exc:
        raise ValueError(f"could not read CSV header: {member_name}") from exc
    if not header or any(not str(name).strip() for name in header):
        raise ValueError(f"CSV has missing or blank column names: {member_name}")
    normalized = [str(name).strip() for name in header]
    if len({name.casefold() for name in normalized}) != len(normalized):
        raise ValueError(f"CSV has duplicate column names: {member_name}")
    return {
        "member_name": member_name,
        "column_names": normalized,
        "read_scope": "one explicitly named non-landmark CSV header row only",
    }


def build_archive_record(path, *, expected_size_bytes: int = EXPECTED_ARCHIVE_SIZE_BYTES,
                         expected_md5: str = EXPECTED_ARCHIVE_MD5) -> dict:
    """Verify archive bytes, then list central-directory metadata only."""
    integrity = verify_archive_integrity(
        path, expected_size_bytes=expected_size_bytes, expected_md5=expected_md5
    )
    inventory = inventory_zip_members(path)
    return {"record_id": 20548528, "integrity": integrity, **inventory}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=_Path, help="completed local DATASET_FINAL.zip path")
    parser.add_argument("--output", type=_Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expected-size-bytes", type=int, default=EXPECTED_ARCHIVE_SIZE_BYTES)
    parser.add_argument("--expected-md5", default=EXPECTED_ARCHIVE_MD5)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit(f"{args.output} exists; refusing to overwrite")
    record = build_archive_record(
        args.archive, expected_size_bytes=args.expected_size_bytes, expected_md5=args.expected_md5
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print({key: record[key] for key in ("archive_name", "member_count", "file_count", "csv_member_names",
                                        "nifti_member_count", "read_scope")})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
