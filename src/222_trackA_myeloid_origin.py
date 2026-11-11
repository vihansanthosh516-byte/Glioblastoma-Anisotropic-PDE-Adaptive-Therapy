"""Exploratory Track A analysis: myeloid-origin program across Core and Periphery.

Protocol (frozen before scoring): paper/trackA_myeloid_origin_protocol_2026-11-10.md.
Reuses the sealed A5 loader and QC. Writes summary metrics only.
"""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "output" / "trackA_myeloid_origin_20261110"
PROTOCOL = "paper/trackA_myeloid_origin_protocol_2026-11-10.md"

MG = ("P2RY12", "TMEM119", "CX3CR1", "SLC2A5", "P2RY13", "GPR34", "SALL1", "OLFML3")
BMDM = ("TGFBI", "ITGA4", "CD163", "ANXA1", "LYZ", "S100A8", "S100A9", "FCN1")
NEO = ("SOX2", "OLIG1", "PTPRZ1", "BCAN", "NES", "SOX4", "CKB", "VIM")
MIN_PER_CLASS = 10
MIN_DONORS = 4
MIN_SET_GENES = 4
N_NULL = 1000
SEED = 20261110


def load_runner():
    spec = importlib.util.spec_from_file_location("trackA_a5_runner", PROJECT_ROOT / "src" / "189_trackA_yu2020_lodo.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["trackA_a5_runner"] = mod
    spec.loader.exec_module(mod)
    return mod


def symbol_index(gene_ids):
    index = {}
    for i, gid in enumerate(map(str, gene_ids)):
        for token in {gid, gid.split("|")[-1], gid.split("|")[0], gid.split("_")[-1]}:
            index.setdefault(token.upper(), i)
    return index


def set_columns(genes, index):
    found = [g for g in genes if g in index]
    return found, [index[g] for g in found]


def zscore(x):
    mu, sd = x.mean(0), x.std(0)
    z = np.zeros_like(x)
    ok = sd > 0
    z[:, ok] = (x[:, ok] - mu[ok]) / sd[ok]
    return z


def auc(y, s):
    return float(roc_auc_score(y, s))


def main():
    runner = load_runner()
    dataset, strict, _, flow, hashes = runner.load_dataset()
    myeloid = np.array([c == "Myeloid" for c in dataset.cell_types])
    primary = np.isin(dataset.donors, strict)
    keep = myeloid & primary
    if not keep.any():
        raise SystemExit("no Myeloid cells found; check cell_type labels")
    counts = dataset.counts[keep]
    y_all = dataset.labels[keep]
    donors_all = dataset.donors[keep]
    eligibility = {}
    for d in strict:
        m = donors_all == d
        eligibility[d] = {"core": int(np.sum(y_all[m] == 1)), "periphery": int(np.sum(y_all[m] == 0))}
    eligible = [d for d in strict if min(eligibility[d].values()) >= MIN_PER_CLASS]
    index = symbol_index(dataset.gene_ids)
    sets = {name: set_columns(g, index) for name, g in (("MG", MG), ("BMDM", BMDM), ("NEO", NEO))}
    doc = {"protocol": PROTOCOL, "status": None, "eligibility": eligibility, "eligible_donors": eligible,
           "genes_found": {k: v[0] for k, v in sets.items()},
           "genes_missing": {k: [g for g in G if g not in sets[k][0]] for k, G in (("MG", MG), ("BMDM", BMDM), ("NEO", NEO))},
           "input_hashes": hashes}
    if len(eligible) < MIN_DONORS or len(sets["MG"][1]) < MIN_SET_GENES or len(sets["BMDM"][1]) < MIN_SET_GENES:
        doc["status"] = "not_run_gate_failed"
        return write(doc)

    totals = counts.sum(1, dtype=np.float64)
    logx = np.log1p(counts.astype(np.float64) * (10_000.0 / totals[:, None]))
    rng = np.random.default_rng(SEED)
    mg_cols, bm_cols, neo_cols = sets["MG"][1], sets["BMDM"][1], sets["NEO"][1]
    excluded = set(mg_cols) | set(bm_cols)
    per_donor = {}
    for d in eligible:
        m = donors_all == d
        x, y = logx[m], y_all[m]
        detected = np.flatnonzero((counts[m] > 0).any(0))
        z = np.zeros_like(x)
        z[:, detected] = zscore(x[:, detected])
        score = z[:, mg_cols].mean(1) - z[:, bm_cols].mean(1)
        prog = auc(y, -score)
        neo = auc(y, z[:, neo_cols].mean(1)) if neo_cols else None
        means = x[:, detected].mean(0)
        deciles = np.searchsorted(np.quantile(means, np.linspace(0.1, 0.9, 9)), means)
        decile_of = dict(zip(detected.tolist(), deciles.tolist()))
        pool = {k: [g for g, dk in decile_of.items() if dk == k and g not in excluded] for k in range(10)}

        def draw(cols):
            out = []
            for c in cols:
                k = decile_of.get(c, 0)
                cand = pool[k] if pool[k] else [g for g in detected.tolist() if g not in excluded]
                out.append(int(rng.choice(cand)))
            return out

        null = np.empty(N_NULL)
        for b in range(N_NULL):
            s = z[:, draw(mg_cols)].mean(1) - z[:, draw(bm_cols)].mean(1)
            null[b] = auc(y, -s)
        per_donor[d] = {"n_core": int((y == 1).sum()), "n_periphery": int((y == 0).sum()),
                        "program_auc": prog, "null_q95": float(np.quantile(null, 0.95)),
                        "null_mean": float(null.mean()), "neo_auc": neo}

    for d in eligible:
        train = np.isin(donors_all, [e for e in eligible if e != d])
        test = donors_all == d
        var = logx[train].var(0)
        top = np.argsort(-var, kind="stable")[:2000]
        xtr, xte = logx[train][:, top], logx[test][:, top]
        mu, sd = xtr.mean(0), xtr.std(0)
        sd[sd == 0] = 1.0
        lr = LogisticRegression(C=1.0, solver="lbfgs", max_iter=2000)
        lr.fit((xtr - mu) / sd, y_all[train])
        per_donor[d]["logistic_auc"] = auc(y_all[test], lr.predict_proba((xte - mu) / sd)[:, 1])

    n = len(eligible)
    above_half = sum(per_donor[d]["program_auc"] > 0.5 for d in eligible)
    above_null = sum(per_donor[d]["program_auc"] > per_donor[d]["null_q95"] for d in eligible)
    above_neo = sum(per_donor[d]["neo_auc"] is not None and per_donor[d]["program_auc"] > per_donor[d]["neo_auc"]
                    for d in eligible)
    at_or_below = n - above_half
    if above_half >= 4 and above_null >= 4 and above_neo >= 4:
        verdict = "positive"
    elif at_or_below >= 3:
        verdict = "negative"
    else:
        verdict = "inconclusive"
    doc.update({
        "status": "completed", "verdict": verdict, "per_donor": per_donor,
        "counts": {"program_auc_gt_0.5": above_half, "program_gt_null_q95": above_null,
                   "program_gt_neo": above_neo, "n_eligible": n,
                   "neo_auc_gt_0.5": sum(per_donor[d]["neo_auc"] > 0.5 for d in eligible),
                   "program_gt_logistic": sum(per_donor[d]["program_auc"] > per_donor[d]["logistic_auc"] for d in eligible)},
        "equal_donor_means": {k: float(np.mean([per_donor[d][k] for d in eligible]))
                              for k in ("program_auc", "null_q95", "neo_auc", "logistic_auc")},
        "git_head": _git(),
    })
    write(doc)


def _git():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return None


def write(doc):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "result.json"
    path.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    print(f"wrote {path} status={doc['status']} verdict={doc.get('verdict')}")


if __name__ == "__main__":
    main()
