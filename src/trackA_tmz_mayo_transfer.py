"""Track A TMZ: in-vitro (HGCC) to in-vivo (Mayo PDX) transfer beyond MGMT.

Protocol: paper/trackA_tmz_mayo_transfer_prereg_2026-11-24.md (locked before run).
One run. Refuses to overwrite its result. `--selftest` runs on generated data only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent
HGCC = PROJECT_ROOT / "data" / "hgcc"
MAYO_EXPR = PROJECT_ROOT / "data" / "external" / "trackA_mayo_pdx_cbio" / "data_mrna_seq_rpkm.txt"
SUPP = PROJECT_ROOT / "data" / "external" / "trackA_mayo_pdx_supp"
PANEL = PROJECT_ROOT / "paper" / "trackA_tmz_response_121_gene_panel_2026-11-07.tsv"
EXCLUDED_SOURCE = {"AF131215.2", "PPIAP22"}
OUT = PROJECT_ROOT / "output" / "trackA_tmz_mayo_transfer" / "result.json"

N_TOP = 2000
ALPHAS = [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0]
N_PERM = 10000
N_RAND = 1000
SEED = 20261124


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def col_ranks(mat: np.ndarray) -> np.ndarray:
    """genes x samples -> within-sample percentile ranks."""
    n = mat.shape[0]
    return np.apply_along_axis(lambda c: rankdata(c, method="average") / n, 0, mat)


def ridge_fit(X, y, alpha):
    xm, ym = X.mean(0), y.mean()
    Xc = X - xm
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    w = Vt.T @ ((s / (s**2 + alpha)) * (U.T @ (y - ym)))
    return w, ym - xm @ w


def ridge_cv(X, y, rng):
    folds = np.array_split(rng.permutation(len(y)), 5)
    best, best_err = None, np.inf
    for a in ALPHAS:
        err = 0.0
        for k in range(5):
            te = folds[k]
            tr = np.setdiff1d(np.arange(len(y)), te)
            w, b = ridge_fit(X[tr], y[tr], a)
            err += float(((X[te] @ w + b - y[te]) ** 2).sum())
        if err < best_err:
            best, best_err = a, err
    w, b = ridge_fit(X, y, best)
    return w, b, best


def resid(v, z):
    A = np.column_stack([np.ones(len(z)), z])
    return v - A @ np.linalg.lstsq(A, v, rcond=None)[0]


def partial_rho(score, y, z):
    rs, ry = rankdata(score), rankdata(y)
    a, b = resid(rs, z), resid(ry, z)
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def perm_p(score, y, z, rng):
    obs = partial_rho(score, y, z)
    cnt = 0
    groups = [np.where(z == g)[0] for g in np.unique(z)]
    for _ in range(N_PERM):
        s = score.copy()
        for g in groups:
            s[g] = s[rng.permutation(g)]
        if partial_rho(s, y, z) >= obs - 1e-12:
            cnt += 1
    return obs, (1 + cnt) / (1 + N_PERM)


def pipeline(h_r, h_target, m_r, gene_idx, rng):
    """Fit ridge on HGCC ranks (origins x genes) for gene subset; score Mayo. Returns benefit score."""
    Xh, Xm = h_r[:, gene_idx], m_r[:, gene_idx]
    mu, sd = Xh.mean(0), Xh.std(0)
    sd = np.where(sd == 0, 1.0, sd)
    Xh, Xm = (Xh - mu) / sd, (Xm - mu) / sd
    w, b, a = ridge_cv(Xh, h_target, rng)
    return -(Xm @ w + b), a  # higher AUC = more viable = resistant; negate for benefit


def analyse(h_genes, h_expr, h_auc, m_genes, m_expr, y, mgmt, subtype, panel_pos, panel_neg):
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    shared = sorted(set(h_genes) & set(m_genes))
    hi = {g: i for i, g in enumerate(h_genes)}
    mi = {g: i for i, g in enumerate(m_genes)}
    H = h_expr[[hi[g] for g in shared]]
    M = m_expr[[mi[g] for g in shared]]
    h_r, m_r = col_ranks(H).T, col_ranks(M).T  # samples x genes
    top = np.argsort(-h_r.var(0), kind="stable")[:N_TOP]
    si = {g: i for i, g in enumerate(shared)}
    h_mgmt = H[si["MGMT"]]
    target_res = resid(h_auc, h_mgmt)

    c1, a1 = pipeline(h_r, target_res, m_r, top, rng)
    c5, a5 = pipeline(h_r, h_auc, m_r, top, rng)
    z = mgmt.astype(float)
    p_c1_rho, p_c1 = perm_p(c1, y, z, rng)

    comp = {}
    comp["1_mgmt_status"] = partial_rho(z, y, z)
    comp["2_mgmt_expression"] = partial_rho(-M[si["MGMT"]], y, z)  # MGMT expression = resistance
    comp["3_subtype_best_indicator"] = max(
        partial_rho((subtype == s).astype(float), y, z) for s in np.unique(subtype)
    )
    pp = [si[g] for g in panel_pos if g in si]
    pn = [si[g] for g in panel_neg if g in si]
    comp["4_ntafoulis_102"] = partial_rho(-(m_r[:, pp].mean(1) - m_r[:, pn].mean(1)), y, z)
    comp["5_no_residual_ablation"] = partial_rho(c5, y, z)

    rand = []
    for _ in range(N_RAND):
        g = rng.choice(len(shared), N_TOP, replace=False)
        s, _a = pipeline(h_r, target_res, m_r, g, rng)
        rand.append(partial_rho(s, y, z))
    rand = np.array(rand)
    rand_pct = float((rand < p_c1_rho).mean())

    pass_rule = bool(
        p_c1_rho > 0 and p_c1 < 0.05 and all(p_c1_rho > v for v in comp.values()) and rand_pct >= 0.95
    )
    return {
        "n_shared_genes": len(shared),
        "n_hgcc_origins": int(len(h_auc)),
        "n_mayo_models": int(len(y)),
        "mgmt_counts": {"M": int(z.sum()), "U": int((1 - z).sum())},
        "alpha_c1": a1,
        "alpha_c5": a5,
        "c1_partial_rho": p_c1_rho,
        "c1_perm_p_one_sided": p_c1,
        "c1_raw_spearman": float(spearmanr(c1, y).statistic),
        "comparators_partial_rho": comp,
        "comparator_7_zhao2019": "not run (gene list and weights not obtained before lock)",
        "random_sets_partial_rho_median": float(np.median(rand)),
        "random_sets_partial_rho_95pct": float(np.percentile(rand, 95)),
        "c1_percentile_vs_random": rand_pct,
        "verdict": "PASS" if pass_rule else "FAIL",
    }


def load_hgcc():
    with open(HGCC / "HGCC_RNA_expression.txt", encoding="latin-1") as fh:
        hdr = [x.strip('"') for x in fh.readline().rstrip("\n").split("\t")]
        genes, rows = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genes.append(f[0].strip('"'))
            rows.append([float(v) if v not in ("NA", "") else np.nan for v in f[1:]])
    mat = np.array(rows)
    cnt = {}
    for g in genes:
        cnt[g] = cnt.get(g, 0) + 1
    keep = [i for i, g in enumerate(genes) if g and cnt[g] == 1 and np.all(np.isfinite(mat[i]))]
    genes, mat = [genes[i] for i in keep], mat[keep]
    with open(HGCC / "HGCC_drug_response_AUC.txt", encoding="latin-1") as fh:
        rd = list(csv.reader(fh, delimiter="\t"))
    dh = rd[0]
    tmz = [r for r in rd[1:] if r[0] == "Temozolomide"]
    assert len(tmz) == 1
    cols = {dh[j]: j for j in range(1, len(dh)) if re.fullmatch(r"U\d{4}", dh[j])}
    origins = sorted(set(cols) & set(hdr))
    assert len(origins) == 83, len(origins)
    expr = np.column_stack([mat[:, [k for k, s in enumerate(hdr) if s == o]].mean(1) for o in origins])
    auc = np.array([float(tmz[0][cols[o]]) for o in origins])
    return genes, expr, auc


def load_mayo_labels():
    rows = [r for r in csv.reader(open(SUPP / "table_s7_parsed.csv")) if r and r[0].isdigit()]
    tmz = {r[0]: (float(r[4]), r[1]) for r in rows if re.match(r"^[\d.]+$", r[4])}
    import openpyxl

    s6 = list(
        openpyxl.load_workbook(SUPP / "NIHMS1546970-supplement-7.xlsx", read_only=True, data_only=True)[
            "Table S6A"
        ].iter_rows(values_only=True)
    )[2:]
    sub = {str(r[0]): str(r[1]) for r in s6 if r[0] is not None}
    return tmz, sub


def load_mayo_expr(ids):
    with open(MAYO_EXPR) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        cols = {h.split("-")[-1]: j for j, h in enumerate(hdr) if h.startswith("Mayo-PDX")}
        use = [i for i in ids if i in cols]
        genes, rows = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genes.append(f[0])
            rows.append([float(f[cols[i]]) if f[cols[i]] not in ("NA", "") else np.nan for i in use])
    mat = np.log2(np.array(rows) + 1.0)
    cnt = {}
    for g in genes:
        cnt[g] = cnt.get(g, 0) + 1
    keep = [i for i, g in enumerate(genes) if g and cnt[g] == 1 and np.all(np.isfinite(mat[i]))]
    return use, [genes[i] for i in keep], mat[keep]


def panel_genes():
    rows = [r for r in csv.DictReader(open(PANEL, encoding="utf-8"), delimiter="\t")
            if r["gene_symbol"] not in EXCLUDED_SOURCE]
    pos = [r["gene_symbol"] for r in rows if r["direction"] == "resistance_associated"]
    neg = [r["gene_symbol"] for r in rows if r["direction"] == "sensitivity_associated"]
    return pos, neg


def selftest():
    global N_PERM, N_RAND, N_TOP
    N_PERM, N_RAND, N_TOP = 500, 30, 150
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(300)] + ["MGMT"]
    sig = rng.normal(size=(len(genes), 1))
    h = rng.normal(size=(len(genes), 83))
    m = rng.normal(size=(len(genes), 27))
    w = np.zeros(len(genes)); w[:40] = 1
    h[:40] *= 3; m[:40] *= 3  # signal genes are high-variance, as in real data
    w = w / 3 / np.sqrt(8)
    h_auc = -(w @ h) / 7 + rng.normal(scale=0.5, size=83)
    mg = (rng.random(27) < 0.5).astype(int)
    y = (w @ m) / 7 + mg + rng.normal(scale=0.5, size=27)
    sub = np.array(["Classical", "Proneural", "Mesenchymal"] * 9)
    r = analyse(genes, h, h_auc, genes, m, y, mg, sub, genes[50:60], genes[60:70])
    print(json.dumps({k: r[k] for k in ("c1_partial_rho", "c1_perm_p_one_sided", "c1_percentile_vs_random", "comparators_partial_rho", "verdict")}))
    y0 = mg + rng.normal(scale=0.5, size=27)
    r0 = analyse(genes, h, h_auc, genes, m, y0, mg, sub, genes[50:60], genes[60:70])
    print(json.dumps({k: r0[k] for k in ("c1_partial_rho", "c1_perm_p_one_sided", "verdict")}))
    assert r["verdict"] == "PASS" and r0["verdict"] == "FAIL"
    print("selftest ok")
    return 0


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    if OUT.exists():
        print("refusing to overwrite", OUT)
        return 1
    tmz, sub = load_mayo_labels()
    h_genes, h_expr, h_auc = load_hgcc()
    ids = sorted(tmz, key=int)
    use, m_genes, m_expr = load_mayo_expr(ids)
    y = np.log(np.array([tmz[i][0] for i in use]))
    mgmt = np.array([1 if tmz[i][1] == "M" else 0 for i in use])
    subtype = np.array([sub.get(i, "NA") for i in use])
    pos, neg = panel_genes()
    res = analyse(h_genes, h_expr, h_auc, m_genes, m_expr, y, mgmt, subtype, pos, neg)
    res["mayo_model_ids"] = use
    res["input_sha256"] = {p.name: sha(p) for p in [MAYO_EXPR, HGCC / "HGCC_RNA_expression.txt",
                                                    HGCC / "HGCC_drug_response_AUC.txt",
                                                    SUPP / "table_s7_parsed.csv"]}
    res["code_sha256"] = sha(Path(__file__))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
