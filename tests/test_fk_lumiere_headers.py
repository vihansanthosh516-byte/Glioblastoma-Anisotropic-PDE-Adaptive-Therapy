"""Synthetic tests for the header-only LUMIERE reader. No patient data."""
import gzip
import io
import struct
import zipfile

import numpy as np

from src import fk_lumiere_headers as lh


def _nifti(shape=(4, 5, 6), spacing=(1.0, 1.0, 1.0), datatype=2, payload=b""):
    hdr = bytearray(348)
    struct.pack_into("<i", hdr, 0, 348)
    struct.pack_into("<8h", hdr, 40, 3, *shape, 1, 1, 1, 1)
    struct.pack_into("<h", hdr, 70, datatype)
    struct.pack_into("<8f", hdr, 76, 1.0, *spacing, 0, 0, 0, 0)
    struct.pack_into("<2h", hdr, 252, 1, 1)
    hdr[344:348] = b"n+1\x00"
    return gzip.compress(bytes(hdr) + b"\x00" * 4 + payload)


def _zip(tmp_path, weeks, *, poison=False, missing=None):
    path = tmp_path / "Imaging.zip"
    with zipfile.ZipFile(path, "w") as zf:
        for w in weeks:
            base = f"Imaging/Patient-001/week-{w:03d}/"
            for s in lh.SEQUENCES:
                zf.writestr(base + s, _nifti())
            for f in lh.HEADER_FILES:
                if f == missing:
                    continue
                # A poisoned payload would break any reader that decodes voxels; headers must still parse.
                zf.writestr(base + f, _nifti(payload=b"\xff" * 64 if poison else b""))
    return path


def test_header_parse_reads_shape_and_spacing():
    h = lh.parse_nifti1_header(gzip.decompress(_nifti((182, 218, 182), (1, 1, 1)))[:348])
    assert h["shape"] == [182, 218, 182] and h["spacing_mm"] == [1.0, 1.0, 1.0] and h["magic"] == "n+1"


def test_reader_reads_only_header_bytes(tmp_path, monkeypatch):
    path = _zip(tmp_path, [4, 8, 21], poison=True)
    sizes = []
    real = lh.read_header

    def spy(zf, name):
        with zf.open(name) as fh, gzip.GzipFile(fileobj=fh) as gz:
            sizes.append(len(gz.read(lh.NIFTI1_HEADER_BYTES)))
        return real(zf, name)
    monkeypatch.setattr(lh, "read_header", spy)
    rec = lh.header_check(path)
    assert rec["n_eligible"] == 1 and not rec["problems"]
    assert set(sizes) == {lh.NIFTI1_HEADER_BYTES}
    assert rec["patients"]["Patient-001"]["gap_weeks"] == 13


def test_missing_support_mask_makes_visit_unusable(tmp_path):
    path = _zip(tmp_path, [4, 8, 21], missing=lh.BRAIN_MASK)
    rec = lh.header_check(path)
    assert rec["n_eligible"] == 0 and rec["n_excluded"] == 1


def test_gap_outside_window_excluded(tmp_path):
    rec = lh.header_check(_zip(tmp_path, [4, 8, 30]))
    assert rec["n_eligible"] == 0
