"""Track A HGCC TMZ bridge, amendment 1 (102-gene score).

Protocol: paper/trackA_hgcc_tmz_response_bridge_2026-11-12.md
Amendment: paper/trackA_hgcc_amendment1_102genes_2026-11-12.md
Writes aggregate statistics only. Refuses to overwrite an existing result.
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
DATA = PROJECT_ROOT / "data" / "hgcc"
PANEL = PROJECT_ROOT / "paper" / "trackA_tmz_response_121_gene_panel_2026-11-07.tsv"
AMENDMENT = PROJECT_ROOT / "paper" / "trackA_hgcc_amendment1_102genes_2026-11-12.md"
OUT = PROJECT_ROOT / "output" / "trackA_hgcc_tmz_102" / "result.json"
EXCLUDED_SOURCE = {"AF131215.2", "PPIAP22"}
N_PERM, N_BOOT, N_CTRL, K_NEAR = 100_000, 10_000, 1_000, 20


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def percentile_ranks(mat: np.ndarray) -> np.ndarray:
    """Column-wise percentile rank (average ties, rank / n). mat: genes x profiles."""
    n = mat.shape[0]
    return np.apply_along_axis(lambda c: rankdata(c, method="average") / n, 0, mat)


def panel_score(pr: np.ndarray, pos_idx, neg_idx) -> np.ndarray:
    return pr[pos_idx].mean(axis=0) - pr[neg_idx].mean(axis=0)


def perm_p_one_sided(x, y, n_perm, seed) -> tuple[float, float]:
    rx = rankdata(x)
    ry = rankdata(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    cnt = 0
    for _ in range(n_perm):
        if np.corrcoef(rx, rng.permutation(ry))[0, 1] >= rho - 1e-12:
            cnt += 1
    return rho, (1 + cnt) / (1 + n_perm)


def boot_ci(x, y, n_boot, seed) -> tuple[float, float]:
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    n = len(x)
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        r = spearmanr(x[i], y[i]).statistic
        if np.isfinite(r):
            vals.append(r)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def control_panels(feat: np.ndarray, panel_idx, n_pos, n_ctrl, seed):
    """feat: genes x 2 standardized (mean, sd). Returns list of (pos_idx, neg_idx)."""
    panel_set = set(panel_idx)
    nonpanel = np.array([g for g in range(feat.shape[0]) if g not in panel_set])
    near = {}
    for g in panel_idx:
        d = np.sqrt(((feat[nonpanel] - feat[g]) ** 2).sum(axis=1))
        order = np.lexsort((nonpanel, d))  # distance, then index (lexicographic tie-break)
        near[g] = nonpanel[order[:500]]
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    out = []
    for _ in range(n_ctrl):
        used, chosen = set(), {}
        for j in rng.permutation(len(panel_idx)):
            g = panel_idx[j]
            avail = [c for c in near[g] if c not in used][:K_NEAR]
            c = int(avail[rng.integers(0, len(avail))])
            used.add(c)
            chosen[j] = c
        ctrl = [chosen[j] for j in range(len(panel_idx))]
        out.append((ctrl[:n_pos], ctrl[n_pos:]))
    return out


def load_expression():
    with open(DATA / "HGCC_RNA_expression.txt", encoding="latin-1") as fh:
        hdr = [x.strip('"') for x in fh.readline().rstrip("\n").split("\t")]
        genes, rows = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(hdr) + 1:
                raise ValueError("unexpected expression row width")
            g, vals = f[0].strip('"'), f[1:]
            genes.append(g)
            rows.append([float(v) if v not in ("NA", "") else np.nan for v in vals])
    return hdr, genes, np.array(rows, dtype=float)


def main() -> int:
    if OUT.exists():
        print("refusing to overwrite", OUT)
        return 1
    samples, genes, mat = load_expression()
    counts = {}
    for g in genes:
        counts[g] = counts.get(g, 0) + 1
    keep = [i for i, g in enumerate(genes) if g and counts[g] == 1 and np.all(np.isfinite(mat[i]))]
    genes = [genes[i] for i in keep]
    mat = mat[keep]
    gidx = {g: i for i, g in enumerate(genes)}

    panel = list(csv.DictReader(open(PANEL, encoding="utf-8"), delimiter="\t"))
    panel = [r for r in panel if r["gene_symbol"] not in EXCLUDED_SOURCE and r["gene_symbol"] in gidx]
    pos = [gidx[r["gene_symbol"]] for r in panel if r["direction"] == "resistance_associated"]
    neg = [gidx[r["gene_symbol"]] for r in panel if r["direction"] == "sensitivity_associated"]
    assert len(pos) + len(neg) == 102, len(pos) + len(neg)

    with open(DATA / "HGCC_drug_response_AUC.txt", encoding="latin-1") as fh:
        rd = list(csv.reader(fh, delimiter="\t"))
    dh = rd[0]
    tmz_row = [r for r in rd[1:] if r[0] == "Temozolomide"]
    assert len(tmz_row) == 1
    exact_cols = {dh[j]: j for j in range(1, len(dh)) if re.fullmatch(r"U\d{4}", dh[j])}
    origins = sorted(set(exact_cols) & set(samples))
    assert len(origins) == 83, len(origins)

    # Expression-only steps first (scores, controls), outcome loaded after.
    pr = percentile_ranks(mat)
    prof_score = panel_score(pr, pos, neg)
    oidx = {o: [k for k, s in enumerate(samples) if s == o] for o in origins}
    score = np.array([prof_score[oidx[o]].mean() for o in origins])
    origin_expr = np.column_stack([mat[:, oidx[o]].mean(axis=1) for o in origins])
    origin_pr = percentile_ranks(origin_expr)
    gm, gs = origin_expr.mean(axis=1), origin_expr.std(axis=1, ddof=1)
    gs = np.where(gs == 0, 1.0, gs)
    feat = np.column_stack([(gm - gm.mean()) / gm.std(), (gs - gs.mean()) / gs.std()])
    ctrls = control_panels(feat, pos + neg, len(pos), N_CTRL, 20261114)
    ctrl_scores = [
        np.array([panel_score(pr, p, n)[oidx[o]].mean() for o in origins]) for p, n in ctrls
    ]
    mgmt = origin_expr[gidx["MGMT"]] if "MGMT" in gidx else None

    # Outcome.
    auc = np.array([float(tmz_row[0][exact_cols[o]]) for o in origins])
    rho, p = perm_p_one_sided(score, auc, N_PERM, 20261112)
    lo, hi = boot_ci(score, auc, N_BOOT, 20261113)
    ctrl_rho = np.array([spearmanr(c, auc).statistic for c in ctrl_scores])
    if rho > 0 and p < 0.05:
        verdict = "supports_in_vitro_association"
    elif rho < 0 and hi < 0:
        verdict = "opposite_direction"
    else:
        verdict = "no_association_shown"
    res = {
        "status": "completed_amendment1_102_gene_score",
        "n_origins": len(origins),
        "n_expression_profiles_used": int(sum(len(v) for v in oidx.values())),
        "n_genes_ranked": len(genes),
        "panel_genes_used": len(pos) + len(neg),
        "n_resistance_associated": len(pos),
        "n_sensitivity_associated": len(neg),
        "spearman_rho": rho,
        "perm_p_one_sided_positive": p,
        "n_permutations": N_PERM,
        "bootstrap_95ci": [lo, hi],
        "n_bootstrap": N_BOOT,
        "verdict": verdict,
        "secondary_mgmt_expression_rho": None if mgmt is None else float(spearmanr(mgmt, auc).statistic),
        "secondary_control_rho_percentile_of_score": float((ctrl_rho < rho).mean() * 100),
        "secondary_control_rho_median": float(np.median(ctrl_rho)),
        "secondary_control_rho_95range": [float(np.percentile(ctrl_rho, 2.5)), float(np.percentile(ctrl_rho, 97.5))],
        "seeds": {"perm": 20261112, "boot": 20261113, "controls": 20261114},
        "amendment_sha256": sha(AMENDMENT),
        "input_sha256": {f.name: sha(f) for f in [DATA / "HGCC_RNA_expression.txt", DATA / "HGCC_drug_response_AUC.txt"]},
        "code_sha256": sha(Path(__file__)),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
