"""Track A TMZ: frozen repair-axis score, human test on GSE7696 (TMZ+RT vs RT).

Protocol: paper/trackA_tmz_gse7696_prereg_LOCKED_2026-11-24.md (locked before run).
Primary: Cox interaction arm x score on overall survival, adjusted for MGMT status and age.
One run. Refuses to overwrite its result. `--selftest` uses generated data only.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import norm

warnings.filterwarnings("ignore")
PROJECT_ROOT = Path(__file__).resolve().parent.parent
D = PROJECT_ROOT / "data" / "external" / "trackA_gse7696"
MATRIX = D / "GSE7696_series_matrix.txt.gz"
ANNOT = D / "GPL570.annot.gz"
OUT = PROJECT_ROOT / "output" / "trackA_tmz_gse7696" / "result.json"
GENES = ["MGMT", "MSH2", "MSH6", "MLH1", "PMS2"]
MMR = GENES[1:]
N_RAND = 1000
SEED = 20261125


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def zs(v):
    s = v.std()
    return (v - v.mean()) / (s if s > 0 else 1.0)


def repair(expr, mmr_genes):
    return zs(-expr["MGMT"]) + zs(np.mean([zs(expr[g]) for g in mmr_genes], axis=0))


def interaction_z(score, df):
    d = df[["time", "event", "arm", "mgmt", "age"]].copy()
    d["score"] = zs(score)
    d["arm_x_score"] = d["arm"] * d["score"]
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(d, "time", "event")
    return float(cph.params_["arm_x_score"]), float(cph.summary.loc["arm_x_score", "z"])


def analyse(expr, df, rand_pool_expr, rng):
    """expr: dict gene -> array; df: time,event,arm(1=TMZ),mgmt(1=M),age."""
    c2 = repair(expr, MMR)
    coef, z = interaction_z(c2, df)
    p = float(norm.cdf(z))  # one-sided: negative coef = more TMZ benefit at higher score
    _, z_abl = interaction_z(-expr["MGMT"], df)
    pool = list(rand_pool_expr)
    rand = []
    for _ in range(N_RAND):
        pick = rng.choice(len(pool), 4, replace=False)
        e = dict(expr)
        for k, j in enumerate(pick):
            e[f"R{k}"] = rand_pool_expr[pool[j]]
        rand.append(interaction_z(repair(e, [f"R{k}" for k in range(4)]), df)[1])
    rand = np.array(rand)
    pct = float((rand > z).mean())  # share of random swaps with weaker (less negative) z
    ok = bool(coef < 0 and p < 0.05 and z < z_abl and pct >= 0.95)
    return {
        "n": int(len(df)),
        "arm_counts": {"TMZ_RT": int(df["arm"].sum()), "RT": int((1 - df["arm"]).sum())},
        "events": int(df["event"].sum()),
        "c2_interaction_coef": coef,
        "c2_interaction_z": z,
        "c2_one_sided_p": p,
        "ablation_mgmt_expr_interaction_z": z_abl,
        "random_swap_z_median": float(np.median(rand)),
        "random_swap_z_5pct": float(np.percentile(rand, 5)),
        "c2_percentile_vs_random": pct,
        "verdict": "PASS" if ok else "FAIL",
    }


def load():
    probes = {g: [] for g in GENES}
    with gzip.open(ANNOT, "rt", encoding="latin-1") as fh:
        hdr = None
        for line in fh:
            if line[0] in "#!^":
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                si = hdr.index("Gene symbol")
                continue
            syms = f[si].split("///") if len(f) > si else []
            for g in GENES:
                if g in syms:
                    probes[g].append(f[0])
    meta, gsm, data = {}, None, {}
    keep_keys = ("disease status", "patient", "treatment", "mgmt status", "age",
                 "survival status", "survival time in months")
    with gzip.open(MATRIX, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsm = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                k = v[0].split(":")[0]
                if k in keep_keys:
                    meta[k] = [x.split(": ", 1)[-1] for x in v]
            elif line.startswith('"ID_REF"') or line.startswith("!series_matrix_table_end"):
                continue
            elif not line.startswith("!") and line.strip():
                f = line.rstrip("\n").split("\t")
                data[f[0].strip('"')] = np.array([float(x) if x not in ("", "null") else np.nan for x in f[1:]])
    seen, idx = set(), []
    for i in range(len(gsm)):
        if meta["disease status"][i] != "GBM":
            continue
        if meta["treatment"][i] not in ("TMZ/radiotherapy", "radiotherapy"):
            continue
        if meta["mgmt status"][i] not in ("M", "U") or meta["patient"][i] in seen:
            continue
        seen.add(meta["patient"][i])
        idx.append(i)
    idx = np.array(idx)
    expr = {g: np.nanmean(np.vstack([data[p][idx] for p in probes[g]]), axis=0) for g in GENES}
    gene_probes = {}
    with gzip.open(ANNOT, "rt", encoding="latin-1") as fh:
        hdr = None
        for line in fh:
            if line[0] in "#!^":
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                si = hdr.index("Gene symbol")
                continue
            s = f[si] if len(f) > si else ""
            if s and "///" not in s and s not in GENES and f[0] in data:
                gene_probes.setdefault(s, []).append(f[0])
    pool = {s: np.mean([data[p][idx] for p in ps], axis=0) for s, ps in gene_probes.items()}
    pool = {s: v for s, v in pool.items() if np.all(np.isfinite(v))}
    df = pd.DataFrame({
        "time": [float(meta["survival time in months"][i]) for i in idx],
        "event": [int(meta["survival status"][i]) for i in idx],
        "arm": [1 if meta["treatment"][i] == "TMZ/radiotherapy" else 0 for i in idx],
        "mgmt": [1 if meta["mgmt status"][i] == "M" else 0 for i in idx],
        "age": [float(meta["age"][i]) for i in idx],
    })
    return expr, df, pool, {g: probes[g] for g in GENES}, [gsm[i] for i in idx]


def selftest():
    global N_RAND
    N_RAND = 100
    res = []
    for kind in ("signal", "null"):
        rng = np.random.default_rng(7 if kind == "signal" else 8)
        n = 300
        expr = {g: rng.normal(size=n) for g in GENES}
        pool = {f"P{i}": rng.normal(size=n) for i in range(200)}
        arm = rng.integers(0, 2, n)
        mg = rng.integers(0, 2, n)
        age = rng.normal(57, 9, n)
        sc = zs(repair(expr, MMR))
        lp = -0.3 * arm - 0.5 * mg + (-0.8 * arm * sc if kind == "signal" else 0)
        t = rng.exponential(np.exp(-lp) * 15)
        c = rng.uniform(5, 60, n)
        df = pd.DataFrame({"time": np.minimum(t, c), "event": (t <= c).astype(int), "arm": arm, "mgmt": mg, "age": age})
        r = analyse(expr, df, pool, np.random.Generator(np.random.PCG64DXSM(1)))
        print(kind, json.dumps({k: r[k] for k in ("c2_interaction_z", "c2_one_sided_p", "ablation_mgmt_expr_interaction_z", "c2_percentile_vs_random", "verdict")}))
        res.append(r["verdict"])
    assert res == ["PASS", "FAIL"], res
    print("selftest ok")
    return 0


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    if OUT.exists():
        print("refusing to overwrite", OUT)
        return 1
    expr, df, pool, probes, gsms = load()
    res = analyse(expr, df, pool, np.random.Generator(np.random.PCG64DXSM(SEED)))
    res["probes"] = probes
    res["samples"] = gsms
    res["input_sha256"] = {MATRIX.name: sha(MATRIX), ANNOT.name: sha(ANNOT)}
    res["code_sha256"] = sha(Path(__file__))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k not in ("samples",)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
