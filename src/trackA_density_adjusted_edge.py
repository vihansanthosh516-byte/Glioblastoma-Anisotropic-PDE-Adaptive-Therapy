"""Track A region: density-adjusted edge signature.

Published edge/core signatures built from malignant cells track tumor-cell density in mixed tissue
(see paper/trackA_region_followups_2026-11-24.md). Method: remove a tumor-content component from the
signature score, then read edge vs core from the residual.

--dev  : development on spent GBM-Space sections (learns and freezes the density proxy gene set).
--ivy  : one locked test on Ivy GAP anatomic-structure RNA-seq (LE vs CT; IT vs CT secondary).
Protocol: paper/trackA_density_adjusted_edge_prereg_LOCKED_2026-11-24.md.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import trackA_region_transfer as T  # noqa: E402

PROXY_FILE = ROOT / "output" / "trackA_density_adjusted_edge" / "density_proxy_genes.json"
DEV_OUT = ROOT / "output" / "trackA_density_adjusted_edge" / "dev_result.json"
IVY_OUT = ROOT / "output" / "trackA_density_adjusted_edge" / "ivy_result.json"
IVY = ROOT / "data" / "external" / "ivygap"
N_PROXY = 50
N_PERM = 10000
SEED = 20261125
MAL = re.compile(r"(-like)|(^Hypoxic)")


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def is_malignant(name: str) -> bool:
    return bool(MAL.search(name)) and "TAMs" not in name


def zcols(X):
    sd = X.std(0)
    return (X - X.mean(0)) / np.where(sd == 0, 1, sd)


def sig_core_score(X, genes, sig):
    """Higher = more core. Direction +1 = edge gene, -1 = core gene."""
    gi = {g: i for i, g in enumerate(genes)}
    use = [(gi[g], d) for g, d in sig if g in gi]
    Z = zcols(X[:, [i for i, _ in use]])
    d = np.array([d for _, d in use], float)
    return -(Z * d).mean(1)


def residual(score, dens):
    A = np.column_stack([np.ones(len(dens)), dens])
    return score - A @ np.linalg.lstsq(A, score, rcond=None)[0]


def proxy_score(X, genes, proxy):
    gi = {g: i for i, g in enumerate(genes)}
    idx = [gi[g] for g in proxy if g in gi]
    return zcols(X[:, idx]).mean(1)


# ------------------------------------------------------------------ development
def dev():
    import h5py

    sigs = T.load_sigs()
    sig_genes = {g for s in sigs.values() for g, _ in s}
    rows = [l.rstrip("\n").split("\t") for l in open(T.COVER)]
    names = [r[1] for r in rows if r[2] == "True"]
    # 1) learn density proxy from density only (no labels): pooled Spearman with malignant fraction
    acc = {}
    cnt = {}
    secs = []
    for n in names:
        f = h5py.File(T.DATA / n, "r")
        ft = T._strs(f["var"]["feature_types"])
        nm = T._strs(f["var"]["_index"])
        f.close()
        s = T.load_section(T.DATA / n)
        cs = nm[ft == "Cell state abundances"]
        mal = np.array([is_malignant(c) for c in cs])
        tot = s["cell"].sum(1)
        dens = s["cell"][:, mal].sum(1) / np.where(tot == 0, 1, tot)
        y = T.label(s.pop("hist"))
        rX = np.apply_along_axis(rankdata, 0, s["X"])
        rd = rankdata(dens)
        r = ((rX - rX.mean(0)) * (rd - rd.mean())[:, None]).sum(0) / (
            np.sqrt(((rX - rX.mean(0)) ** 2).sum(0)) * np.sqrt(((rd - rd.mean()) ** 2).sum()) + 1e-12)
        for g, v in zip(s["genes"], r):
            acc[g] = acc.get(g, 0.0) + v
            cnt[g] = cnt.get(g, 0) + 1
        m = y >= 0
        if (y == 1).sum() >= 20 and (y == 0).sum() >= 20:
            secs.append((n, s["X"][m], s["genes"], dens[m], y[m]))
    mean_r = {g: acc[g] / cnt[g] for g in acc if cnt[g] >= len(names) // 2 and g not in sig_genes}
    proxy = [g for g, _ in sorted(mean_r.items(), key=lambda kv: -kv[1])[:N_PROXY]]
    PROXY_FILE.parent.mkdir(parents=True, exist_ok=True)
    PROXY_FILE.write_text(json.dumps({"genes": proxy, "mean_spearman": [mean_r[g] for g in proxy]}, indent=1))
    # 2) development AUROCs (spent labels)
    out = []
    for n, X, genes, dens, y in secs:
        r = {"section": n}
        pr = proxy_score(X, genes, proxy)
        r["density_true_auc"] = float(roc_auc_score(y, dens))
        r["density_proxy_auc"] = float(roc_auc_score(y, pr))
        r["proxy_vs_true_density_rho"] = float(spearmanr(pr, dens).statistic)
        for k, sg in sigs.items():
            raw = sig_core_score(X, genes, sg)
            r[f"{k}_raw"] = float(roc_auc_score(y, raw))
            r[f"{k}_adj_true"] = float(roc_auc_score(y, residual(raw, dens)))
            r[f"{k}_adj_proxy"] = float(roc_auc_score(y, residual(raw, pr)))
        out.append(r)
    keys = [k for k in out[0] if k != "section"]
    summ = {k: float(np.mean([r[k] for r in out])) for k in keys}
    summ["n_sections"] = len(out)
    DEV_OUT.write_text(json.dumps({"status": "development, spent data", "summary": summ, "sections": out}, indent=2))
    print(json.dumps(summ, indent=2))


# ------------------------------------------------------------------ Ivy GAP test
def load_ivy():
    m = IVY / "matrix"
    genes = pd.read_csv(m / "rows-genes.csv")
    fpkm = pd.read_csv(m / "fpkm_table.csv", index_col=0)
    fpkm.columns = fpkm.columns.astype(int)
    det = pd.read_csv(IVY / "rna_seq_samples_details.csv")
    det = det[det.structure_acronym.isin(["LE-reference-histology", "IT-reference-histology", "CT-reference-histology"])]
    det = det[det.sample_id.isin(fpkm.columns)]
    x = np.log2(fpkm[det.sample_id.tolist()] + 1)
    x.index = genes.set_index("gene_id").loc[x.index, "gene_symbol"].values
    x = x.groupby(level=0).mean()
    zone = det.structure_acronym.str[:2].values
    return x.T.values, list(x.index), zone, det.donor_id.values


def perm_p(score, y, donors, rng):
    """One-sided p for AUROC > observed under label shuffling within donor."""
    obs = roc_auc_score(y, score)
    groups = [np.where(donors == d)[0] for d in np.unique(donors)]
    cnt = 0
    for _ in range(N_PERM):
        yy = y.copy()
        for g in groups:
            yy[g] = yy[rng.permutation(g)]
        if len(np.unique(yy)) == 2 and roc_auc_score(yy, score) >= obs - 1e-12:
            cnt += 1
    return float(obs), (1 + cnt) / (1 + N_PERM)


def ivy():
    if IVY_OUT.exists():
        print("refusing to overwrite", IVY_OUT)
        return 1
    proxy = json.loads(PROXY_FILE.read_text())["genes"]
    sigs = T.load_sigs()
    X, genes, zone, donors = load_ivy()
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    res = {"n_samples": {z: int((zone == z).sum()) for z in ("LE", "IT", "CT")},
           "n_donors": int(len(np.unique(donors)))}
    for contrast in ("LE", "IT"):
        keep = (zone == "CT") | (zone == contrast)
        Xk, yk, dk = X[keep], (zone[keep] == "CT").astype(int), donors[keep]
        pr = proxy_score(Xk, genes, proxy)
        block = {"density_proxy_auc": float(roc_auc_score(yk, pr)),
                 "n_ct": int(yk.sum()), "n_other": int((1 - yk).sum()),
                 "n_donors": int(len(np.unique(dk)))}
        for k, sg in sigs.items():
            raw = sig_core_score(Xk, genes, sg)
            adj = residual(raw, pr)
            block[f"{k}_raw_auc"] = float(roc_auc_score(yk, raw))
            if k == "saadawy" and contrast == "LE":
                a, p = perm_p(adj, yk, dk, rng)
                block["saadawy_adj_auc"], block["saadawy_adj_perm_p"] = a, p
            else:
                block[f"{k}_adj_auc"] = float(roc_auc_score(yk, adj))
        res[f"{contrast}_vs_CT"] = block
    b = res["LE_vs_CT"]
    res["verdict"] = "PASS" if (b["saadawy_adj_auc"] > 0.5 and b["saadawy_adj_perm_p"] < 0.05
                                and b["saadawy_adj_auc"] > b["saadawy_raw_auc"]) else "FAIL"
    res["proxy_sha256"] = sha(PROXY_FILE)
    res["code_sha256"] = sha(Path(__file__))
    res["input_sha256"] = {p.name: sha(p) for p in [IVY / "rna_seq_samples_details.csv",
                                                    IVY / "matrix" / "fpkm_table.csv"]}
    IVY_OUT.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    if "--dev" in sys.argv:
        dev()
    elif "--ivy" in sys.argv:
        sys.exit(ivy())
    else:
        print("use --dev or --ivy")
