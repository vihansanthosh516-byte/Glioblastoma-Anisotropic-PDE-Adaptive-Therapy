"""Script 173: Track A sealed test -- SenMayo program, labeled edge vs core.

Protocol: paper/trackA_senmayo/PROTOCOL.md (sealed in commit c02217b).
Preprocessing choices not fixed by the protocol, fixed here before any score is read:
  Ivy GAP: log2(FPKM + 1); genes by symbol (rows-genes.csv), duplicate symbols averaged.
  Angers: Agilent gProcessedSignal, log2(x + 1), quantile-normalized across all
  arrays in the archive, probes collapsed to GeneName by mean.
  Genes kept per cohort: present and non-zero variance across that cohort's used samples.
Writes output/trackA_senmayo/result.json.
"""
import hashlib
import itertools
import json
import zipfile
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
SEAL = PROJECT_ROOT / "paper" / "trackA_senmayo"
EXT = PROJECT_ROOT / "data" / "external"
OUTPUT_DIR = PROJECT_ROOT / "output" / "trackA_senmayo"
SEED = 20261029
N_NULL = 1000

IVY_DONORS = [10865, 10926, 12108, 12109, 12110, 12112, 14734, 292023109]
ANGERS_PTS = ["FT02", "FT04", "FT05", "FT06", "FT08", "FT09", "FT13"]


def check_seal():
    for line in (SEAL / "SHA256SUMS").read_text().splitlines():
        h, name = line.split()
        data = (SEAL / name.lstrip("*")).read_bytes().replace(b"\r\n", b"\n")
        if hashlib.sha256(data).hexdigest() != h:
            raise SystemExit(f"seal mismatch: {name}")


def read_list(name):
    return [g.strip() for g in (SEAL / name).read_text().splitlines() if g.strip()]


def load_ivy():
    m = EXT / "ivygap" / "matrix"
    genes = pd.read_csv(m / "rows-genes.csv")
    fpkm = pd.read_csv(m / "fpkm_table.csv", index_col=0)
    fpkm.columns = fpkm.columns.astype(int)
    det = pd.read_csv(EXT / "ivygap" / "rna_seq_samples_details.csv")
    det = det[det.donor_id.isin(IVY_DONORS)
              & det.structure_acronym.isin(["LE-reference-histology", "CT-reference-histology"])]
    det = det[det.sample_id.isin(fpkm.columns)]
    x = np.log2(fpkm[det.sample_id.tolist()] + 1)
    x.index = genes.set_index("gene_id").loc[x.index, "gene_symbol"].values
    x = x.groupby(level=0).mean()
    zone = {"LE-reference-histology": "edge", "CT-reference-histology": "core"}
    meta = pd.DataFrame({"patient": det.donor_id.astype(str).values,
                         "zone": det.structure_acronym.map(zone).values},
                        index=det.sample_id.values)
    return x, meta


def quantile_normalize(df):
    ranks = df.rank(method="first").astype(int) - 1
    means = np.sort(df.values, axis=0).mean(axis=1)
    return pd.DataFrame(means[ranks.values], index=df.index, columns=df.columns)


def load_angers():
    d = EXT / "angers"
    sdrf = pd.read_csv(d / "E-MTAB-5092.sdrf.txt", sep="\t")
    ind, site = sdrf.columns[4], sdrf.columns[5]
    cols = {}
    zips = [zipfile.ZipFile(d / f"E-MTAB-5092.raw.{i}.zip") for i in (1, 2)]
    if True:
        for f in sdrf["Array Data File"]:
            z = next(z for z in zips if f in z.namelist())
            with z.open(f) as fh:
                lines = fh.read().decode("latin-1").splitlines()
            start = next(i for i, l in enumerate(lines) if l.startswith("FEATURES"))
            rows = [l.split("\t") for l in lines[start:]]
            t = pd.DataFrame(rows[1:], columns=rows[0])
            t = t[t.ControlType == "0"]
            s = pd.Series(t.gProcessedSignal.astype(float).values, index=t.ProbeName.values)
            cols[f] = s.groupby(level=0).mean()  # replicate spots of one probe
            if "genes" not in cols:
                g = pd.Series(t.GeneName.values, index=t.ProbeName.values)
                cols["genes"] = g[~g.index.duplicated()]
    genes = cols.pop("genes")
    raw = pd.DataFrame(cols).dropna()
    qn = quantile_normalize(np.log2(raw.clip(lower=0) + 1))
    qn.index = genes.loc[qn.index].values
    x = qn.groupby(level=0).mean()
    zone = {"interface": "edge", "tumor zone": "core", "peripheral brain zone": "pbz"}
    s = sdrf[sdrf[ind].isin(ANGERS_PTS) & sdrf[site].isin(zone)]
    meta = pd.DataFrame({"patient": s[ind].values, "zone": s[site].map(zone).values},
                        index=s["Array Data File"].values)
    return x[meta.index], meta


def zone_diffs(x, meta, sig, hyp, rng, edge="edge"):
    """Per-patient E for signature and hypoxia, with matched null."""
    x = x[x.var(axis=1) > 0]
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
    pts = sorted(set(meta.patient[meta.zone == edge]) & set(meta.patient[meta.zone == "core"]))
    zone_mean = {}
    for p in pts:
        for zn in (edge, "core"):
            zone_mean[(p, zn)] = z[meta.index[(meta.patient == p) & (meta.zone == zn)]].mean(axis=1)
    delta = pd.DataFrame({p: zone_mean[(p, edge)] - zone_mean[(p, "core")] for p in pts})

    mean_expr = x.mean(axis=1)
    dec = pd.qcut(mean_expr.rank(method="first"), 10, labels=False)
    out = {}
    for name, gs in (("sig", sig), ("hyp", hyp)):
        present = [g for g in gs if g in z.index]
        pool = {b: [g for g in dec.index[dec == b] if g not in set(present)] for b in range(10)}
        need = dec.loc[present].value_counts()
        nulls = np.empty((N_NULL, len(pts)))
        for i in range(N_NULL):
            draw = np.concatenate([rng.choice(pool[b], n, replace=False) for b, n in need.items()])
            nulls[i] = delta.loc[draw].mean().values
        obs = delta.loc[present].mean().values
        out[name] = {"E": obs - np.median(nulls, axis=0), "n_genes": len(present),
                     "D_sig": obs, "D_null_median": np.median(nulls, axis=0)}
    return pts, out


def sign_flip(e):
    e = np.asarray(e, float)
    obs = abs(e.mean())
    stats = [abs((np.array(s) * e).mean()) for s in itertools.product([1, -1], repeat=len(e))]
    return float(np.mean(np.array(stats) >= obs - 1e-12)), float(e.mean())


def scale(e):
    return e / np.sqrt(np.mean(e ** 2))


def main():
    check_seal()
    sig, hyp = read_list("senmayo_genes.txt"), read_list("hallmark_hypoxia_genes.txt")
    rng = np.random.default_rng(SEED)
    cohorts = {"ivygap": load_ivy(), "angers": load_angers()}
    res, E, H = {}, [], []
    for c, (x, meta) in cohorts.items():
        pts, o = zone_diffs(x, meta, sig, hyp, rng)
        res[c] = {"patients": pts, "n_sig_genes": o["sig"]["n_genes"],
                  "n_hyp_genes": o["hyp"]["n_genes"],
                  "E_sig": o["sig"]["E"].tolist(), "E_hyp": o["hyp"]["E"].tolist(),
                  "D_sig": o["sig"]["D_sig"].tolist(),
                  "D_null_median": o["sig"]["D_null_median"].tolist(),
                  "mean_E_sig": float(o["sig"]["E"].mean())}
        E.append(scale(o["sig"]["E"]))
        H.append(o["hyp"]["E"])
    e, h = np.concatenate(E), np.concatenate(H)
    p, m = sign_flip(e)
    b = float((e @ h) / (h @ h))
    e_adj = np.concatenate([scale(a - b * hh) for a, hh in zip(E, H)])
    p_adj, m_adj = sign_flip(e_adj)
    same_sign = np.sign(res["ivygap"]["mean_E_sig"]) == np.sign(res["angers"]["mean_E_sig"])
    primary = {"n_patients": len(e), "p": p, "mean_scaled_E": m, "p_floor": 2 / 2 ** len(e),
               "cohorts_same_sign": bool(same_sign), "hypoxia_b": b,
               "p_hypoxia_adjusted": p_adj, "mean_scaled_E_adjusted": m_adj,
               "direction": "edge > core" if m > 0 else "core > edge"}
    primary["success"] = bool(p < 0.05 and same_sign and p_adj < 0.05
                              and np.sign(m_adj) == np.sign(m))

    sec = {}
    for c, (x, meta) in cohorts.items():
        sec[c + "_only"] = dict(zip(("p", "mean"), sign_flip(scale(np.array(res[c]["E_sig"])))))
    x, meta = cohorts["angers"]
    pts, o = zone_diffs(x, meta, sig, hyp, np.random.default_rng(SEED), edge="pbz")
    sec["angers_pbz_minus_tz"] = {"patients": pts, **dict(zip(("p", "mean"), sign_flip(scale(o["sig"]["E"]))))}

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = {"protocol_commit": "c02217b", "seed": SEED, "primary": primary,
           "cohorts": res, "secondaries": sec}
    (OUTPUT_DIR / "result.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
