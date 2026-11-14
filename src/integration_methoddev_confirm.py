"""One locked confirmation: does the TCGA-trained MRI model predict the Neftel MES score on new patients?

Lock: paper/integration_methoddev_lock_2026-11-12.md (target neftel_MES, features 'new', model pls2).
Protocol: paper/integration_protocol_methoddev_2026-11-12.md. Run once. Written before any confirm RNA is read.

Choices fixed before running (not in the protocol text):
- CPTAC: GDC STAR 'tpm_unstranded', log2(TPM+1), duplicate gene names averaged.
- IvyGAP: cellular tumor = structure_abbreviation starting with 'CT-' (CT-reference-histology, CT-<gene>,
  CT-control-*). CTmvp / CTpan / CTpnz / CThbv are other structures and are excluded. Patient = tumor_name
  before the first '-'. log2(FPKM+1), mean over that patient's CT samples.
- Score: same Neftel module scoring as development (g1.state_scores), computed within each cohort.

  py -3.14 -I src/integration_methoddev_confirm.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import integration_bamd_gate1 as g1  # noqa: E402
import integration_methoddev_search as ms  # noqa: E402

OUT = PROJECT_ROOT / "output" / "integration_methoddev"
EXT = PROJECT_ROOT / "data" / "external"
LOCK = PROJECT_ROOT / "paper" / "integration_methoddev_lock_2026-11-12.md"
TARGET, MODEL = "neftel_MES", "pls2"
N_PERM, SEED = 10_000, 20261114
PASS_RHO, PASS_P = 0.25, 0.05


def cptac_expr() -> pd.DataFrame:
    d = EXT / "cptac_gbm_star"
    man = [ln.split() for ln in (d / "ids.txt").read_text().splitlines() if ln.strip()]
    cols = {}
    for _fid, fn, pid in man:
        t = pd.read_csv(d / fn, sep="\t", comment="#")
        t = t[~t["gene_id"].str.startswith("N_")]
        cols[pid] = np.log2(t.groupby("gene_name")["tpm_unstranded"].mean() + 1.0)
    if len(cols) != 23:
        raise SystemExit(f"CPTAC: expected 23 RNA files, found {len(cols)}")
    return pd.DataFrame(cols)


def ivygap_expr() -> pd.DataFrame:
    d = EXT / "ivygap" / "matrix"
    c = pd.read_csv(d / "columns-samples.csv")
    c["pid"] = c["tumor_name"].str.split("-").str[0]
    c = c[c["structure_abbreviation"].str.startswith("CT-")]
    genes = pd.read_csv(d / "rows-genes.csv").set_index("gene_id")["gene_symbol"]
    f = pd.read_csv(d / "fpkm_table.csv", index_col=0)
    f.columns = f.columns.astype(int)
    f = np.log2(f + 1.0)
    f.index = genes.reindex(f.index).to_numpy()
    f = f.groupby(level=0).mean()
    return pd.DataFrame({pid: f[g["rna_well_id"].tolist()].mean(axis=1) for pid, g in c.groupby("pid")})


def fisher_combined(rhos, ns):
    w = np.asarray(ns, float) - 3
    return float(np.tanh(np.sum(w * np.arctanh(rhos)) / w.sum()))


def main() -> int:
    if not LOCK.exists():
        raise SystemExit("lock note missing")
    d, T, _old, newcols = ms.load()
    Xtr = d[newcols].to_numpy(float)
    med = np.nanmedian(Xtr, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    model = ms.make_model(MODEL).fit(np.where(np.isnan(Xtr), med, Xtr), T[TARGET].to_numpy(float))

    cohorts = {}
    for name, expr in (("CPTAC-GBM", cptac_expr()), ("IvyGAP", ivygap_expr())):
        feat = pd.read_csv(OUT / f"features_{name}.csv")
        feat["patient_id"] = feat["patient_id"].astype(str)
        feat = feat[feat["patient_id"].isin(expr.columns)].drop_duplicates("patient_id").sort_values("patient_id")
        y = g1.state_scores(expr[feat["patient_id"].tolist()], g1.load_modules())["MES"].to_numpy(float)
        X = feat[newcols].to_numpy(float)
        pred = np.ravel(model.predict(np.where(np.isnan(X), med, X)))
        cohorts[name] = {"pred": pred, "y": y, "ids": feat["patient_id"].tolist()}

    rhos = [spearmanr(c["pred"], c["y"])[0] for c in cohorts.values()]
    ns = [len(c["y"]) for c in cohorts.values()]
    obs = fisher_combined(rhos, ns)
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    null = np.array([fisher_combined([spearmanr(c["pred"], rng.permutation(c["y"]))[0] for c in cohorts.values()], ns)
                     for _ in range(N_PERM)])
    p = float((1 + np.sum(null >= obs)) / (1 + N_PERM))
    passed = bool(obs >= PASS_RHO and p < PASS_P)
    res = {"target": TARGET, "model": MODEL, "features": "new", "n_train_tcga": len(d),
           "cohorts": {k: {"n": n, "spearman_rho": float(r)} for (k, _), n, r in zip(cohorts.items(), ns, rhos)},
           "combined_rho_fisher": obs, "p_one_sided_perm": p, "n_perm": N_PERM,
           "pass_rule": f"combined rho >= {PASS_RHO} and p < {PASS_P}", "passed": passed,
           "label": "one locked confirmation; real patients; MRI vs RNA association, not treatment benefit",
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "lock_sha256": hashlib.sha256(LOCK.read_bytes()).hexdigest()}
    (OUT / "confirmation.json").write_text(json.dumps(res, indent=2))
    pd.concat([pd.DataFrame({"cohort": k, "patient_id": c["ids"], "pred": c["pred"], "mes_score": c["y"]})
               for k, c in cohorts.items()]).to_csv(OUT / "confirmation_predictions.csv", index=False)
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
