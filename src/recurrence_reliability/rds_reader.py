"""Minimal reader for R serialization (RDS, XDR format versions 2 and 3), standard library + numpy only.

Covers the object types found in data-frame metadata files: NULL, symbols, pairlists, character,
logical, integer, double, complex, raw and list vectors, attributes, references, environments, and
ALTREP compact sequences / wrappers. Values are returned as Python objects and never printed here.
Format reference: R Internals, "Serialization formats".
"""
from __future__ import annotations

import gzip
import struct
from dataclasses import dataclass, field

import numpy as np

NA_INT = -2147483648


@dataclass
class RObj:
    rtype: int
    value: object
    attrs: dict = field(default_factory=dict)


class _Reader:
    def __init__(self, data: bytes):
        self.b, self.p, self.refs = data, 0, []

    def take(self, n: int) -> bytes:
        out = self.b[self.p:self.p + n]
        if len(out) != n:
            raise EOFError("truncated RDS stream")
        self.p += n
        return out

    def i32(self) -> int:
        return struct.unpack(">i", self.take(4))[0]

    def length(self) -> int:
        n = self.i32()
        if n == -1:
            hi, lo = self.i32(), self.i32()
            n = (hi << 32) + lo
        return n

    def attrs_from_pairlist(self, pl) -> dict:
        if pl is None:
            return {}
        return {tag: val for tag, val in pl.value}

    def item(self):
        flags = self.i32()
        t = flags & 0xFF
        has_attr, has_tag = bool(flags & (1 << 9)), bool(flags & (1 << 10))
        if t == 254:
            return None
        if t in (242, 241, 253, 252, 251, 247):
            return RObj(t, None)
        if t == 255:
            idx = flags >> 8
            if idx == 0:
                idx = self.i32()
            return self.refs[idx - 1]
        if t == 1:
            name = self.item()
            sym = RObj(1, name)
            self.refs.append(sym)
            return sym
        if t == 4:
            env = RObj(4, None)
            self.refs.append(env)
            self.i32()
            enclos, frame, hashtab, attrib = self.item(), self.item(), self.item(), self.item()
            env.value = {"enclos": enclos, "frame": frame, "hashtab": hashtab}
            env.attrs = self.attrs_from_pairlist(attrib) if isinstance(attrib, RObj) and attrib.rtype == 2 else {}
            return env
        if t in (2, 3, 5, 6, 17):
            items = []
            attrs = {}
            while True:
                a = self.item() if has_attr else None
                if a is not None and not attrs:
                    attrs = self.attrs_from_pairlist(a)
                tag = self.item() if has_tag else None
                car = self.item()
                tagname = tag.value if isinstance(tag, RObj) and tag.rtype == 1 else None
                items.append((tagname, car))
                nxt_flags = struct.unpack(">i", self.b[self.p:self.p + 4])[0]
                if (nxt_flags & 0xFF) == t:
                    self.p += 4
                    has_attr, has_tag = bool(nxt_flags & (1 << 9)), bool(nxt_flags & (1 << 10))
                    continue
                cdr = self.item()
                if cdr is not None:
                    items.append((None, cdr))
                break
            return RObj(t, items, attrs)
        if t == 9:
            n = self.i32()
            if n == -1:
                return None
            enc = "latin-1" if flags & (1 << 2) else "utf-8"
            return self.take(n).decode(enc, errors="replace")
        if t == 238:
            info, state, attr = self.item(), self.item(), self.item()
            cls = info.value[0][1].value if isinstance(info, RObj) else None
            if cls in ("compact_intseq", "compact_realseq"):
                n, start, inc = [float(x) for x in state.value]
                vals = start + inc * np.arange(int(n))
                obj = RObj(13 if cls == "compact_intseq" else 14, vals.astype(np.int64 if cls == "compact_intseq" else float))
            elif cls and cls.startswith("wrap_"):
                obj = state.value[0] if isinstance(state.value, list) else state
                obj = obj[1] if isinstance(obj, tuple) else obj
            elif cls == "deferred_string":
                src = state.value[0][1] if isinstance(state.value, list) else state
                obj = RObj(16, [None if v == NA_INT else str(v) for v in np.asarray(src.value).tolist()])
            else:
                raise ValueError(f"unsupported ALTREP class {cls}")
            if isinstance(attr, RObj):
                obj.attrs = self.attrs_from_pairlist(attr)
            return obj
        if t in (10, 13):
            n = self.length()
            vals = np.frombuffer(self.take(4 * n), dtype=">i4").astype(np.int64)
            obj = RObj(t, vals)
        elif t == 14:
            n = self.length()
            obj = RObj(t, np.frombuffer(self.take(8 * n), dtype=">f8").astype(float))
        elif t == 15:
            n = self.length()
            obj = RObj(t, np.frombuffer(self.take(16 * n), dtype=">f8").astype(float))
        elif t == 24:
            n = self.length()
            obj = RObj(t, self.take(n))
        elif t == 16:
            n = self.length()
            obj = RObj(t, [self.item() for _ in range(n)])
        elif t in (19, 20):
            n = self.length()
            obj = RObj(t, [self.item() for _ in range(n)])
        elif t == 25:
            obj = RObj(t, None)
        else:
            raise ValueError(f"unsupported SEXP type {t}")
        if has_attr:
            obj.attrs = self.attrs_from_pairlist(self.item())
        return obj


def read_rds_bytes(data: bytes) -> RObj:
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    r = _Reader(data)
    if r.take(2) != b"X\n":
        raise ValueError("only XDR (binary big-endian) RDS is supported")
    version = r.i32()
    r.i32(), r.i32()
    if version == 3:
        r.take(r.i32())
    elif version != 2:
        raise ValueError(f"unsupported serialization version {version}")
    return r.item()


def read_rds(path) -> RObj:
    with open(path, "rb") as f:
        return read_rds_bytes(f.read())


def column_names(df: RObj) -> list[str]:
    return list(df.attrs["names"].value)


def column(df: RObj, name: str) -> list:
    """Decode one data-frame column to Python values (factor levels applied)."""
    col = df.value[column_names(df).index(name)]
    levels = col.attrs.get("levels")
    if levels is not None:
        lv = levels.value
        return [None if v == NA_INT else lv[v - 1] for v in col.value.tolist()]
    if col.rtype == 16:
        return list(col.value)
    return np.asarray(col.value).tolist()
