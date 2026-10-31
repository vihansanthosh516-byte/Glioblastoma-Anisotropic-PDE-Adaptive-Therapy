"""Inventory IvyGAP histopathology feature names in GBM-Space Visium H5AD files.

This reads HDF5 metadata only: /obs/_index (spot count), /var/_index,
/var/feature_types, and the feature-type category/code vectors. It never reads
/X, /layers, image data, or annotation values. HTTP reads use byte ranges.

Dependencies: h5py, fsspec[http]. Example:
  python src/tracka_gbmspace_feature_inventory.py --output output/inventory.json
Resume an existing output explicitly with --resume.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import csv
import json
import os
import tempfile
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import fsspec
import h5py
import numpy as np


ACCESSION = "E-MTAB-17164"
API_URL = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{ACCESSION}"
FILES_ROOT = f"https://www.ebi.ac.uk/biostudies/files/{ACCESSION}/"
SDRF_URL = f"{FILES_ROOT}{ACCESSION}.sdrf.txt"
SCHEMA_VERSION = 2
MAX_WORKERS = 4
BLOCK_SIZE = 1_048_576
REQUEST_TIMEOUT_SECONDS = 90
PREFLIGHT_REQUEST_TIMEOUT_SECONDS = 60
HISTOPATH_TYPE = "Histopath annotation overlap"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def get_json(url: str) -> dict[str, Any]:
    import aiohttp

    with fsspec.open(
        url,
        "rt",
        encoding="utf-8",
        client_kwargs={"timeout": aiohttp.ClientTimeout(total=PREFLIGHT_REQUEST_TIMEOUT_SECONDS)},
    ) as stream:
        return json.load(stream)


def walk_h5ad_files(value: Any):
    if isinstance(value, dict):
        if isinstance(value.get("path"), str) and value["path"].endswith(".h5ad"):
            yield {"name": value["path"], "size_bytes": value.get("size")}
        for key in ("files", "subsections", "section"):
            if key in value:
                yield from walk_h5ad_files(value[key])
    elif isinstance(value, list):
        for child in value:
            yield from walk_h5ad_files(child)


def get_sdrf_map() -> dict[str, dict[str, str]]:
    import aiohttp

    with fsspec.open(
        SDRF_URL,
        "rt",
        encoding="utf-8",
        client_kwargs={"timeout": aiohttp.ClientTimeout(total=PREFLIGHT_REQUEST_TIMEOUT_SECONDS)},
    ) as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("SDRF is empty")
    header = rows[0]
    try:
        donor_i = header.index("Characteristics[individual]")
    except ValueError as exc:
        raise ValueError("SDRF lacks Characteristics[individual]") from exc
    try:
        sample_i = header.index("Characteristics[sample_id]")
    except ValueError:
        sample_i = -1

    mapping: dict[str, dict[str, str]] = {}
    for row in rows[1:]:
        if len(row) != len(header):
            continue  # Blank trailing line is not a sample.
        h5ad_name = next((cell for cell in row if cell.endswith(".h5ad")), None)
        if not h5ad_name:
            continue
        donor = row[donor_i].strip()
        sample = row[sample_i].strip() if sample_i >= 0 else ""
        if not donor:
            raise ValueError(f"Missing donor for {h5ad_name}")
        if h5ad_name in mapping:
            raise ValueError(f"Duplicate SDRF H5AD mapping: {h5ad_name}")
        mapping[h5ad_name] = {"donor_id": donor, "sample_id": sample}
    return mapping


def decode_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "item"):
        value = value.item()
        if isinstance(value, bytes):
            return value.decode("utf-8")
    return str(value)


def inspect_one(entry: dict[str, Any]) -> dict[str, Any]:
    url = entry["source_url"]
    started = time.monotonic()
    try:
        # aiohttp timeout applies to each HTTP request, not the full file task.
        # Block caching bounds transfer to the HDF5 metadata needed below.
        import aiohttp

        options = {
            "block_size": BLOCK_SIZE,
            "cache_type": "readahead",
            "client_kwargs": {
                "timeout": aiohttp.ClientTimeout(total=REQUEST_TIMEOUT_SECONDS)
            },
        }
        with fsspec.open(url, "rb", **options) as remote:
            declared_size = remote.size
            with h5py.File(remote, "r") as h5:
                obs_index = h5["obs"]["_index"]
                feature_types = h5["var"]["feature_types"]
                categories = [decode_text(x) for x in feature_types["categories"][:]]
                codes = feature_types["codes"][:]
                var_count = len(codes)
                if "Histopath annotation overlap" in categories:
                    histopath_code = categories.index(HISTOPATH_TYPE)
                    histopath_indexes = [int(i) for i in np.flatnonzero(codes == histopath_code)]
                    if histopath_indexes:
                        # HDF5 fancy indexing requires sorted, unique indexes.
                        histopath_names = sorted(
                            decode_text(name)
                            for name in h5["var"]["_index"][histopath_indexes]
                        )
                    else:
                        histopath_names = []
                else:
                    histopath_names = []
                return {
                    **entry,
                    "status": "complete",
                    "http_content_length_bytes": declared_size,
                    "spot_count": len(obs_index),
                    "var_feature_count": var_count,
                    "feature_type_names": sorted(set(categories)),
                    "histopath_feature_names": histopath_names,
                    "histopath_feature_count": len(histopath_names),
                    "elapsed_seconds": round(time.monotonic() - started, 3),
                    "read_scope": [
                        "/obs/_index length only",
                        "/var/feature_types/categories",
                        "/var/feature_types/codes",
                        "/var/_index at histopath-feature indexes only (if present)",
                    ],
                    "feature_names_read": True,
                    "feature_types_read": True,
                    "annotation_matrix_values_read": False,
                    "expression_matrix_values_read": False,
                    "images_read": False,
                    "error": None,
                }
    except Exception as exc:  # Capture per-file failures so run can resume.
        return {
            **entry,
            "status": "error",
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "feature_names_read": False,
            "feature_types_read": False,
            "annotation_matrix_values_read": False,
            "expression_matrix_values_read": False,
            "images_read": False,
            "error": f"{type(exc).__name__}: {exc}",
        }


def atomic_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_name, path)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)


@contextlib.contextmanager
def exclusive_writer(output: Path):
    """Acquire an exclusive sidecar lock; fail closed if any lock already exists."""
    lock_path = output.with_name(output.name + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    lock = {"pid": os.getpid(), "token": token, "created_at_utc": utc_now()}
    encoded = (json.dumps(lock) + "\n").encode("utf-8")
    try:
        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise RuntimeError(
            f"Lock exists; refusing concurrent write. Inspect owner and remove lock explicitly only after it stops: {lock_path}"
        ) from exc
    with os.fdopen(fd, "wb") as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())

    try:
        yield
    finally:
        try:
            current = json.loads(lock_path.read_text(encoding="utf-8"))
            if current.get("token") == token:
                lock_path.unlink()
        except FileNotFoundError:
            pass


def summary(records: list[dict[str, Any]], expected_donors: list[str]) -> dict[str, Any]:
    by_donor: dict[str, list[dict[str, Any]]] = {donor: [] for donor in expected_donors}
    for record in records:
        by_donor.setdefault(record["donor_id"], []).append(record)
    result = {}
    for donor, donor_records in sorted(by_donor.items()):
        complete = [r for r in donor_records if r.get("status") == "complete"]
        features = [set(r.get("histopath_feature_names", [])) for r in complete]
        result[donor] = {
            "section_count": len(donor_records),
            "completed_section_count": len(complete),
            "sections_with_any_histopath_feature": sum(bool(s) for s in features),
            # Presence counts only; values and nonzero overlap were not read.
            "sections_with_cellular_tumor_and_leading_edge_features": sum(
                "Cellular tumor" in s
                and any(name.startswith("Leading edge") for name in s)
                for s in features
            ),
            "sections_with_cellular_tumor_and_infiltrating_tumor_features": sum(
                "Cellular tumor" in s
                and any(name.startswith("Infiltrating tumor") for name in s)
                for s in features
            ),
        }
    return result


def run_inventory(args: argparse.Namespace, output: Path, parser: argparse.ArgumentParser) -> int:
    listing = get_json(API_URL)
    file_entries = sorted(walk_h5ad_files(listing), key=lambda item: item["name"])
    sdrf = get_sdrf_map()
    api_names = [item["name"] for item in file_entries]
    api_name_set = set(api_names)
    sdrf_name_set = set(sdrf)
    if (
        len(file_entries) != 97
        or len(api_name_set) != len(api_names)
        or len(sdrf) != 97
        or api_name_set != sdrf_name_set
    ):
        raise SystemExit(
            "Refusing nonmatching inventory: "
            f"API H5AD rows={len(file_entries)}, unique={len(api_name_set)}, "
            f"SDRF rows={len(sdrf)}; API-only={sorted(api_name_set-sdrf_name_set)[:5]}, "
            f"SDRF-only={sorted(sdrf_name_set-api_name_set)[:5]}"
        )
    records = []
    for item in file_entries:
        if item["name"] not in sdrf:
            raise SystemExit(f"H5AD missing from SDRF: {item['name']}")
        records.append(
            {
                **item,
                **sdrf[item["name"]],
                "source_url": FILES_ROOT + item["name"],
                "status": "pending",
            }
        )

    if args.resume:
        payload = json.loads(output.read_text(encoding="utf-8"))
        if (
            payload.get("schema_version") != SCHEMA_VERSION
            or payload.get("accession") != ACCESSION
            or payload.get("source_api_url") != API_URL
            or payload.get("sdrf_url") != SDRF_URL
            or payload.get("declared_h5ad_file_count") != 97
        ):
            parser.error("Existing JSON schema/accession mismatch; refusing resume")
        old = {r["name"]: r for r in payload.get("files", [])}
        if set(old) != {r["name"] for r in records}:
            parser.error("Existing JSON file inventory differs; refusing resume")
        for current in records:
            saved = old[current["name"]]
            for field in ("donor_id", "sample_id", "source_url", "size_bytes"):
                if saved.get(field) != current.get(field):
                    parser.error(
                        f"Existing JSON metadata mismatch for {current['name']} ({field}); refusing resume"
                    )
            if saved.get("status") == "complete":
                features = saved.get("histopath_feature_names")
                if (
                    not isinstance(features, list)
                    or saved.get("histopath_feature_count") != len(features)
                    or not isinstance(saved.get("spot_count"), int)
                    or not isinstance(saved.get("var_feature_count"), int)
                ):
                    parser.error(
                        f"Existing completed record is incomplete: {current['name']}"
                    )
        records = [
            old[r["name"]] if old[r["name"]].get("status") == "complete" else r
            for r in records
        ]
        created_at = payload.get("created_at_utc", utc_now())
    else:
        created_at = utc_now()
        payload = {}

    payload.update(
        {
            "schema_version": SCHEMA_VERSION,
            "accession": ACCESSION,
            "created_at_utc": created_at,
            "updated_at_utc": utc_now(),
            "status": "running",
            "source_api_url": API_URL,
            "sdrf_url": SDRF_URL,
            "declared_h5ad_file_count": len(records),
            "donor_count_from_sdrf": len({r["donor_id"] for r in records}),
            "max_concurrent_files": MAX_WORKERS,
            "per_request_timeout_seconds": REQUEST_TIMEOUT_SECONDS,
            "preflight_request_timeout_seconds": PREFLIGHT_REQUEST_TIMEOUT_SECONDS,
            "method": "HDF5 metadata-only; fsspec HTTP byte-range reads; reads feature-type codes before selecting only histopath feature names; no expression or annotation matrix values, no images",
            "metadata_scope": "Reads spot-index length, var feature types, and histopath feature names only; no expression/annotation matrix values",
            "feature_semantics": {
                "strict_edge": "Leading edge feature names",
                "intermediate_zone": "Infiltrating tumor feature names; not equivalent to leading edge",
                "core": "Cellular tumor feature name",
                "counts_are_presence_only": True,
                "nonzero_values_verified": False,
                "donor_eligibility_assessed": False,
            },
            "files": records,
            "per_donor_summary": summary(records, sorted({r["donor_id"] for r in records})),
        }
    )
    payload.pop("per_file_timeout_seconds", None)
    atomic_write(output, payload)

    pending = [r for r in records if r.get("status") != "complete"]
    print(f"Inventory: {len(records)} files; complete={len(records)-len(pending)}; pending={len(pending)}; concurrency={MAX_WORKERS}", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(inspect_one, record): record["name"] for record in pending}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            records = [result if r["name"] == result["name"] else r for r in records]
            payload["files"] = records
            payload["per_donor_summary"] = summary(records, sorted({r["donor_id"] for r in records}))
            payload["updated_at_utc"] = utc_now()
            atomic_write(output, payload)
            print(
                f"{result['status']}: {result['name']} donor={result['donor_id']} "
                f"features={result.get('histopath_feature_count', 'NA')} "
                f"error={result.get('error') or '-'}",
                flush=True,
            )

    errors = sum(r.get("status") == "error" for r in records)
    payload["status"] = "complete" if errors == 0 else "partial_with_errors"
    payload["completed_at_utc"] = utc_now()
    payload["error_file_count"] = errors
    payload["per_donor_summary"] = summary(records, sorted({r["donor_id"] for r in records}))
    atomic_write(output, payload)
    print(f"Finished: status={payload['status']} errors={errors} output={output}", flush=True)
    return 0 if errors == 0 else 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--resume", action="store_true", help="Resume this exact output JSON")
    args = parser.parse_args()
    output = args.output.resolve()

    if output.exists() and not args.resume:
        parser.error(f"Output already exists; choose a new path or pass --resume: {output}")
    if not output.exists() and args.resume:
        parser.error(f"Cannot resume missing output: {output}")
    with exclusive_writer(output):
        return run_inventory(args, output, parser)


if __name__ == "__main__":
    raise SystemExit(main())
