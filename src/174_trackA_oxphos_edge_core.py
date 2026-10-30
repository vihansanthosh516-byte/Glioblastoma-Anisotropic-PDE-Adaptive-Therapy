"""Script 174: Track A sealed test 2 -- Hallmark OXPHOS, labeled edge vs core.

Protocol: paper/trackA_oxphos/PROTOCOL.md (sealed in commit 162e55a).
Reuses loaders, null and sign-flip from script 173 unchanged.
Writes output/trackA_oxphos/result.json.
"""
import importlib.util
import json
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "s173", PROJECT_ROOT / "src" / "173_trackA_senmayo_edge_core.py")
s173 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s173)

s173.SEAL = PROJECT_ROOT / "paper" / "trackA_oxphos"
OUTPUT_DIR = PROJECT_ROOT / "output" / "trackA_oxphos"
ALPHA = 0.025


def main():
    s173.check_seal()
    ox = s173.read_list("oxphos_genes.txt")
    hyp = s173.read_list("hallmark_hypoxia_genes.txt")
    neu = s173.read_list("neuron_genes.txt")
    rng = np.random.default_rng(s173.SEED)
    cohorts = {"ivygap": s173.load_ivy(), "angers": s173.load_angers()}
    res, E, H, N = {}, [], [], []
    for c, (x, meta) in cohorts.items():
        pts, o = s173.zone_diffs(x, meta, ox, hyp, rng)
        _, on = s173.zone_diffs(x, meta, neu, neu, rng)
        res[c] = {"patients": pts, "n_sig_genes": o["sig"]["n_genes"],
                  "n_neuron_genes": on["sig"]["n_genes"],
                  "E_sig": o["sig"]["E"].tolist(), "E_hyp": o["hyp"]["E"].tolist(),
                  "E_neuron": on["sig"]["E"].tolist(),
                  "mean_E_sig": float(o["sig"]["E"].mean())}
        E.append(o["sig"]["E"])
        H.append(o["hyp"]["E"])
        N.append(on["sig"]["E"])
    p, m = s173.sign_flip(np.concatenate([s173.scale(e) for e in E]))
    X = np.column_stack([np.concatenate(H), np.concatenate(N)])
    y = np.concatenate(E)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    split = np.cumsum([len(e) for e in E])[:-1]
    p_adj, m_adj = s173.sign_flip(np.concatenate([s173.scale(r) for r in np.split(resid, split)]))
    same = bool(np.sign(res["ivygap"]["mean_E_sig"]) == np.sign(res["angers"]["mean_E_sig"]))
    primary = {"n_patients": len(y), "alpha": ALPHA, "p": p, "mean_scaled_E": m,
               "direction": "edge > core" if m > 0 else "core > edge",
               "cohorts_same_sign": same, "coef_hypoxia": float(coef[0]),
               "coef_neuron": float(coef[1]), "p_adjusted": p_adj,
               "mean_scaled_E_adjusted": m_adj}
    primary["success"] = bool(p < ALPHA and same and p_adj < ALPHA
                              and np.sign(m_adj) == np.sign(m))
    sec = {c + "_only": dict(zip(("p", "mean"), s173.sign_flip(s173.scale(np.array(res[c]["E_sig"])))))
           for c in res}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = {"protocol_commit": "162e55a", "seed": s173.SEED, "primary": primary,
           "cohorts": res, "secondaries": sec}
    (OUTPUT_DIR / "result.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
