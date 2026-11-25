"""Track A region: one labelled Visium section -> other sites of the same GBM (GBM-Space).

Protocol: paper/trackA_region_one_section_transfer_prereg_LOCKED_2026-11-24.md.
One run. Refuses to overwrite its result. `--selftest` uses generated data only.
Histopathology overlap values are read only inside `run()` after the lock.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA = PROJECT_ROOT / "data" / "external" / "gbmspace_visium"
COVER = PROJECT_ROOT / "data" / "external" / "trackA_gbmspace_metadata" / "histopath_coverage.tsv"
SIG = PROJECT_ROOT / "paper" / "trackA_artzi_signatures"
OUT = PROJECT_ROOT / "output" / "trackA_region_transfer" / "result.json"

N_HVG = 2000
N_PC = 30
K_NN = 15
MIN_CLASS = 20
ALPHA = 0.10
SEED = 20261124
NAME_RE = re.compile(r"^(AT\d+)-.*-(\d+)_(\d+)\.h5ad$")
CT, IT, LE = "Cellular tumor", "Infiltrating tumor (grey matter)", "Leading edge (grey matter)"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ---------------- features ----------------
def logcpm(counts: np.ndarray) -> np.ndarray:
    lib = counts.sum(1, keepdims=True)
    lib[lib == 0] = 1
    return np.log1p(counts / lib * 1e4)


def within_rank(x: np.ndarray) -> np.ndarray:
    """Per-gene quantile rank across spots of one section (absorbs section-level shifts)."""
    return np.apply_along_axis(lambda c: rankdata(c) / len(c), 0, x)


def fit_logit(Xs, ys, Xt, rank: bool):
    if rank:
        Xs, Xt = within_rank(Xs), within_rank(Xt)
    mu, sd = Xs.mean(0), Xs.std(0)
    sd[sd == 0] = 1
    Xs, Xt = (Xs - mu) / sd, (Xt - mu) / sd
    pca = PCA(n_components=min(N_PC, Xs.shape[1], Xs.shape[0] - 1), random_state=SEED).fit(Xs)
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000).fit(pca.transform(Xs), ys)
    return clf, pca, mu, sd, rank


def predict(model, X):
    clf, pca, mu, sd, rank = model
    if rank:
        X = within_rank(X)
    return clf.predict_proba(pca.transform((X - mu) / sd))[:, 1]


def knn_transfer(Xs, ys, Xt):
    mu, sd = Xs.mean(0), Xs.std(0)
    sd[sd == 0] = 1
    pca = PCA(n_components=min(N_PC, Xs.shape[1]), random_state=SEED).fit((Xs - mu) / sd)
    nn = NearestNeighbors(n_neighbors=K_NN).fit(pca.transform((Xs - mu) / sd))
    _, idx = nn.kneighbors(pca.transform((Xt - mu) / sd))
    return ys[idx].mean(1)


def marker_core_score(X, genes, sig):
    gi = {g: i for i, g in enumerate(genes)}
    use = [(gi[g], d) for g, d in sig if g in gi]
    if not use:
        return np.zeros(X.shape[0])
    Z = X[:, [i for i, _ in use]]
    Z = (Z - Z.mean(0)) / np.where(Z.std(0) == 0, 1, Z.std(0))
    d = np.array([d for _, d in use], dtype=float)
    return -(Z * d).mean(1)  # direction +1 = edge, -1 = core; higher = core


def best_single(As, ys, At):
    best, sign, j = -1, 1, 0
    for k in range(As.shape[1]):
        if np.ptp(As[:, k]) == 0:
            continue
        a = roc_auc_score(ys, As[:, k])
        if abs(a - 0.5) > best:
            best, sign, j = abs(a - 0.5), (1 if a >= 0.5 else -1), k
    return sign * At[:, j]


def logit_simple(As, ys, At):
    mu, sd = As.mean(0), As.std(0)
    sd[sd == 0] = 1
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000).fit((As - mu) / sd, ys)
    return clf.predict_proba((At - mu) / sd)[:, 1]


def conformal_cov(p_cal, y_cal, p_t, y_t):
    """Split-conformal (LAC) sets at level ALPHA; coverage and singleton rate on target."""
    s = 1 - np.where(y_cal == 1, p_cal, 1 - p_cal)
    n = len(s)
    q = np.quantile(s, min(1.0, np.ceil((n + 1) * (1 - ALPHA)) / n), method="higher")
    in1 = (1 - p_t) <= q
    in0 = p_t <= q
    cov = np.where(y_t == 1, in1, in0).mean()
    single = (in1 ^ in0).mean()
    return float(cov), float(single)


# ---------------- evaluation over sections ----------------
def evaluate(sections, sigs):
    """sections: dict name -> dict(donor, site, rep, genes, X (logcpm), cell, niche, y (1 core, 0 edge, -1 other))."""
    by_donor = {}
    for n, s in sections.items():
        by_donor.setdefault(s["donor"], []).append(n)
    arms = ["candidate_rank_logit", "ablation_plain_logit", "knn_transfer", "marker_saadawy",
            "marker_bastola", "niche_best_single", "cellstate_logit"]
    pairs, donors_out = [], {}
    for d, names in sorted(by_donor.items()):
        ok = [n for n in names if (sections[n]["y"] == 1).sum() >= MIN_CLASS and (sections[n]["y"] == 0).sum() >= MIN_CLASS]
        if len({sections[n]["site"] for n in ok}) < 2:
            continue
        dp = []
        for src in ok:
            S = sections[src]
            ms = S["y"] >= 0
            genes = S["genes"]
            hv = np.argsort(-S["X"][ms].var(0), kind="stable")[:N_HVG]
            hvg = [genes[i] for i in hv]
            Xs = S["X"][ms][:, hv]
            ys = S["y"][ms]
            cand = fit_logit(Xs, ys, Xs, rank=True)
            abl = fit_logit(Xs, ys, Xs, rank=False)
            cal = [n for n in ok if n != src and sections[n]["site"] == S["site"]]
            for tgt in ok:
                T = sections[tgt]
                if T["site"] == S["site"]:
                    continue
                mt = T["y"] >= 0
                ti = {g: i for i, g in enumerate(T["genes"])}
                cols = [ti.get(g, -1) for g in hvg]
                Tm = T["X"][mt]
                Xt = np.zeros((Tm.shape[0], len(cols)), dtype=Tm.dtype)
                ci_ok = [k for k, c in enumerate(cols) if c >= 0]
                Xt[:, ci_ok] = Tm[:, [cols[k] for k in ci_ok]]
                del Tm
                yt = T["y"][mt]
                sc = {
                    "candidate_rank_logit": predict(cand, Xt),
                    "ablation_plain_logit": predict(abl, Xt),
                    "knn_transfer": knn_transfer(Xs, ys, Xt),
                    "marker_saadawy": marker_core_score(T["X"][mt], T["genes"], sigs["saadawy"]),
                    "marker_bastola": marker_core_score(T["X"][mt], T["genes"], sigs["bastola"]),
                    "niche_best_single": best_single(S["niche"][ms], ys, T["niche"][mt]),
                    "cellstate_logit": logit_simple(S["cell"][ms], ys, T["cell"][mt]),
                }
                row = {"donor": d, "source": src, "target": tgt, "n_target": int(len(yt)),
                       "auroc": {k: float(roc_auc_score(yt, v)) for k, v in sc.items()}}
                if cal:
                    C = sections[cal[0]]
                    mc = C["y"] >= 0
                    ci = {g: i for i, g in enumerate(C["genes"])}
                    Cm = C["X"][mc]
                    Xc = np.zeros((Cm.shape[0], len(hvg)), dtype=Cm.dtype)
                    ok_c = [k for k, g in enumerate(hvg) if g in ci]
                    Xc[:, ok_c] = Cm[:, [ci[hvg[k]] for k in ok_c]]
                    del Cm
                    p_cal = predict(cand, Xc)
                    row["conformal_replicate_cal"] = conformal_cov(p_cal, C["y"][mc], sc["candidate_rank_logit"], yt)
                    row["conformal_insample_cal"] = conformal_cov(predict(cand, Xs), ys, sc["candidate_rank_logit"], yt)
                dp.append(row)
        if dp:
            donors_out[d] = {k: float(np.mean([r["auroc"][k] for r in dp])) for k in arms}
            donors_out[d]["n_pairs"] = len(dp)
            pairs.extend(dp)
    comps = [a for a in arms if a != "candidate_rank_logit"]
    wins = {d: all(v["candidate_rank_logit"] > v[c] for c in comps) for d, v in donors_out.items()}
    verdict = "PASS" if donors_out and len(donors_out) >= 4 and all(wins.values()) else "FAIL"
    cov = [r["conformal_replicate_cal"] for r in pairs if "conformal_replicate_cal" in r]
    cov_in = [r["conformal_insample_cal"] for r in pairs if "conformal_insample_cal" in r]
    return {
        "eligible_donors": sorted(donors_out),
        "donor_mean_auroc": donors_out,
        "candidate_beats_all_by_donor": wins,
        "verdict": verdict,
        "comparator_paste2": "not run (package install not permitted in this environment)",
        "secondary_conformal_target_coverage_replicate_cal_mean": float(np.mean([c[0] for c in cov])) if cov else None,
        "secondary_conformal_target_coverage_insample_cal_mean": float(np.mean([c[0] for c in cov_in])) if cov_in else None,
        "secondary_conformal_singleton_rate_replicate_cal_mean": float(np.mean([c[1] for c in cov])) if cov else None,
        "pairs": pairs,
    }


# ---------------- data ----------------
def _strs(g):
    import h5py

    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
        return cats[g["codes"][:]]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]])


def load_section(path: Path):
    import h5py

    f = h5py.File(path, "r")
    ft = _strs(f["var"]["feature_types"])
    names = _strs(f["var"]["_index"])
    X = f["X"]
    if isinstance(X, h5py.Group):
        M = sparse.csr_matrix((X["data"][:].astype(np.float32), X["indices"][:], X["indptr"][:]),
                              shape=tuple(X.attrs["shape"])).toarray()
    else:
        M = X[:].astype(np.float32)
    g = ft == "Gene Expression"
    G = M[:, g]
    keep = (G > 0).mean(0) >= 0.01  # genes detected in >= 1% of spots
    return {
        "genes": list(names[g][keep]),
        "X": logcpm(G)[:, keep].astype(np.float32),
        "cell": M[:, ft == "Cell state abundances"].astype(float),
        "niche": M[:, ft == "Spatial niche abundances"].astype(float),
        "hist": {n: M[:, i].astype(float) for i, n in enumerate(names) if ft[i] == "Histopath annotation overlap"},
    }


def label(hist):
    n = len(next(iter(hist.values()))) if hist else 0
    if hist and max(v.max() for v in hist.values()) <= 1.0:  # fractions -> percent (scale rule fixed before lock)
        hist = {k: v * 100 for k, v in hist.items()}
    ct = hist.get(CT, np.zeros(n))
    ed = hist.get(IT, np.zeros(n)) + hist.get(LE, np.zeros(n))
    y = -np.ones(len(ct), dtype=int)
    y[ct >= 50] = 1
    y[(ed >= 50) & (ct < 50)] = 0
    return y


def load_sigs():
    sa = [(r["Symbol"], int(r["Direction"])) for r in csv.DictReader(open(SIG / "c26_saadawy_le_vs_cellular_core.tsv"), delimiter="\t")]
    ba = [(r["gene_symbol"], int(r["direction"])) for r in csv.DictReader(open(SIG / "c17_bastola_edge_core_markers.tsv"), delimiter="\t")]
    return {"saadawy": sa, "bastola": ba}


def merge(parts):
    donors, pairs = {}, []
    for p in parts:
        donors.update(p["donor_mean_auroc"])
        pairs.extend(p["pairs"])
    arms = ["candidate_rank_logit", "ablation_plain_logit", "knn_transfer", "marker_saadawy",
            "marker_bastola", "niche_best_single", "cellstate_logit"]
    wins = {d: all(v[arms[0]] > v[c] for c in arms[1:]) for d, v in donors.items()}
    cov = [r["conformal_replicate_cal"] for r in pairs if "conformal_replicate_cal" in r]
    cov_in = [r["conformal_insample_cal"] for r in pairs if "conformal_insample_cal" in r]
    return {
        "eligible_donors": sorted(donors),
        "donor_mean_auroc": donors,
        "candidate_beats_all_by_donor": wins,
        "verdict": "PASS" if len(donors) >= 4 and all(wins.values()) else "FAIL",
        "comparator_paste2": "not run (package install not permitted in this environment)",
        "secondary_conformal_target_coverage_replicate_cal_mean": float(np.mean([c[0] for c in cov])) if cov else None,
        "secondary_conformal_target_coverage_insample_cal_mean": float(np.mean([c[0] for c in cov_in])) if cov_in else None,
        "secondary_conformal_singleton_rate_replicate_cal_mean": float(np.mean([c[1] for c in cov])) if cov else None,
        "pairs": pairs,
    }


def run():
    if OUT.exists():
        print("refusing to overwrite", OUT)
        return 1
    rows = [l.rstrip("\n").split("\t") for l in open(COVER)]
    by_donor = {}
    for _d, name, has, _n in rows:
        m = NAME_RE.match(name)
        if has == "True" and m:
            by_donor.setdefault(m.group(1), []).append((name, int(m.group(2)), int(m.group(3))))
    sigs, used, parts = load_sigs(), [], []
    for donor, items in sorted(by_donor.items()):
        sections = {}
        for name, site, rep in items:
            s = load_section(DATA / name)
            s.update(donor=donor, site=site, rep=rep)
            s["y"] = label(s.pop("hist"))
            sections[name] = s
            used.append(name)
        parts.append(evaluate(sections, sigs))
        del sections
    res = merge(parts)
    res["sections_used"] = used
    res["input_sha256"] = {n: sha(DATA / n) for n in used}
    res["code_sha256"] = sha(Path(__file__))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k not in ("pairs", "input_sha256")}, indent=2))
    return 0


def selftest():
    rng = np.random.default_rng(3)
    genes = [f"G{i}" for i in range(400)] + ["PTPRZ1", "CD44"]
    sigs = {"saadawy": [("PTPRZ1", 1)], "bastola": [("CD44", -1)]}

    def make(donor, site, rep, shift, signal=True):
        n = 400
        y = rng.integers(0, 2, n)
        base = rng.normal(size=(n, len(genes)))
        if signal:
            base[:, :40] += (0.25 * y)[:, None]
        scale = shift[None, :]
        X = np.log1p(np.exp(base) * scale)  # section-specific per-gene scale
        return dict(donor=donor, site=site, rep=rep, genes=genes, X=X, y=y,
                    cell=rng.normal(size=(n, 10)), niche=rng.normal(size=(n, 8)))

    secs = {}
    for d in range(4):
        for site in (1, 2):
            shift = np.exp(rng.normal(scale=1.0, size=len(genes)))
            for rep in (1, 2):
                secs[f"D{d}-{site}_{rep}"] = make(f"D{d}", site, rep, shift)
    r = evaluate(secs, sigs)
    print(json.dumps({"verdict": r["verdict"], "donor_mean_auroc": r["donor_mean_auroc"]}, indent=1))
    secs0 = {k: dict(v, y=rng.integers(0, 2, len(v["y"]))) for k, v in secs.items()}
    r0 = evaluate(secs0, sigs)
    print("null verdict", r0["verdict"])
    assert r["verdict"] == "PASS" and r0["verdict"] == "FAIL"
    print("selftest ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(selftest() if "--selftest" in sys.argv else run())
