"""Synthetic checks for the minimal RDS reader: a hand-serialised R data.frame (XDR, version 3)."""
import gzip
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recurrence_reliability.rds_reader import column, column_names, read_rds_bytes  # noqa: E402

ATTR, TAG, OBJ = 1 << 9, 1 << 10, 1 << 8


def i32(x):
    return struct.pack(">i", x)


def charsxp(s):
    if s is None:
        return i32(9) + i32(-1)
    b = s.encode()
    return i32(9 | (1 << 15)) + i32(len(b)) + b   # UTF8 flag bit is irrelevant to the reader


def strsxp(vals, attrs=b""):
    return i32(16 | (ATTR if attrs else 0)) + i32(len(vals)) + b"".join(charsxp(v) for v in vals) + attrs


def intsxp(vals, attrs=b"", obj=False):
    return i32(13 | (ATTR if attrs else 0) | (OBJ if obj else 0)) + i32(len(vals)) + b"".join(i32(v) for v in vals) + attrs


def realsxp(vals):
    return i32(14) + i32(len(vals)) + b"".join(struct.pack(">d", v) for v in vals)


def pairlist(items, refs):
    out = b""
    for tag, val in items:
        if tag in refs:
            tagbytes = i32(255 | (refs[tag] << 8))
        else:
            refs[tag] = len(refs) + 1
            tagbytes = i32(1) + charsxp(tag)
        out += i32(2 | TAG) + tagbytes + val
    return out + i32(254)


def make_df():
    refs = {}
    levels = pairlist([("levels", strsxp(["Malignant", "TAM"])), ("class", strsxp(["factor"]))], refs)
    cols = [strsxp(["NL022-AAACCCAAGAAGGTAG-1", "NL022-AAACCCAAGAAGGTAA-1", None]),
            intsxp([1, 2, -2147483648], attrs=levels, obj=True),
            realsxp([1.5, 2.0, 3.25])]
    attrs = pairlist([("names", strsxp(["CellID", "CellType", "OSmonth"])),
                      ("class", strsxp(["data.frame"])),
                      ("row.names", intsxp([-2147483648, -3]))], refs)
    body = i32(19 | ATTR | OBJ) + i32(3) + b"".join(cols) + attrs
    header = b"X\n" + i32(3) + i32(0x40300) + i32(0x30500) + i32(5) + b"UTF-8"
    return gzip.compress(header + body)


def test_reads_data_frame_columns_and_factor_levels():
    df = read_rds_bytes(make_df())
    assert column_names(df) == ["CellID", "CellType", "OSmonth"]
    assert column(df, "CellID") == ["NL022-AAACCCAAGAAGGTAG-1", "NL022-AAACCCAAGAAGGTAA-1", None]
    assert column(df, "CellType") == ["Malignant", "TAM", None]
    assert list(df.attrs["class"].value) == ["data.frame"]
    assert np.allclose(column(df, "OSmonth"), [1.5, 2.0, 3.25])


def test_rejects_non_xdr_and_truncated_input():
    with pytest.raises(ValueError):
        read_rds_bytes(gzip.compress(b"A\n" + i32(3)))
    with pytest.raises(EOFError):
        read_rds_bytes(gzip.decompress(make_df())[:-10])
