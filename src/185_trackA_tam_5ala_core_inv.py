"""Script 185: Track A sealed test 3 -- TAM activation, tumour core vs 5-ALA invasive region (E-MTAB-8743).

Protocol: paper/trackA_tam_5ala/PROTOCOL.md (sealed before any expression value was read).
Reuses zone_diffs, sign_flip and scale from script 173 unchanged.
Writes output/trackA_tam_5ala/result.json; refuses to overwrite it.
"""
import importlib.util
import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "s173", PROJECT_ROOT / "src" / "173_trackA_senmayo_edge_core.py")
s173 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s173)

s173.SEAL = PROJECT_ROOT / "paper" / "trackA_tam_5ala"
DATA = PROJECT_ROOT / "data" / "external" / "E-MTAB-8743_author_matrix"
OUTPUT_DIR = PROJECT_ROOT / "output" / "trackA_tam_5ala"
SEED = 20261030
ALPHA = 0.05


def load():
    counts = pd.read_csv(DATA / "Final_count_exon.csv", index_col=0).drop(columns="width")
    ann = pd.read_csv(DATA / "Annotation_col.csv", index_col=0)["Annotation"]
    ann = ann[ann.isin(["Core", "Inv"])]
    meta = pd.DataFrame({"patient": [s.split()[1].split("_")[0] for s in ann.index],
                         "zone": ann.map({"Core": "core", "Inv": "inv"}).values}, index=ann.index)
    counts.index = counts.index.astype(str)          # Entrez gene IDs (amendment 1)
    c = counts[meta.index]
    x = np.log2(c / c.sum() * 1e6 + 1)
    return x, meta


def main():
    s173.check_seal()
    fp = OUTPUT_DIR / "result.json"
    if fp.exists():
        sys.exit(f"refusing to overwrite {fp}: the sealed test is scored once")
    tam, hyp, neu, mye = (s173.read_list(f) for f in ("tam_activation_genes.txt", "hallmark_hypoxia_genes.txt",
                                                      "neuron_genes.txt", "myeloid_genes.txt"))
    ent = pd.read_csv(s173.SEAL / "symbol_to_entrez.tsv", sep="\t", dtype=str).set_index("symbol")["entrez"]
    tam, hyp, neu, mye = ([ent[g] for g in gs] for gs in (tam, hyp, neu, mye))
    x, meta = load()
    assert meta.groupby("patient").zone.nunique().eq(2).all() and meta.patient.nunique() == 10
    rng = np.random.default_rng(SEED)
    pts, o = s173.zone_diffs(x, meta, tam, hyp, rng, edge="inv")
    _, on = s173.zone_diffs(x, meta, neu, mye, rng, edge="inv")
    E, H, N, M = o["sig"]["E"], o["hyp"]["E"], on["sig"]["E"], on["hyp"]["E"]
    p, m = s173.sign_flip(s173.scale(E))
    X = np.column_stack([H, N, M])
    coef, *_ = np.linalg.lstsq(X, E, rcond=None)
    p_adj, m_adj = s173.sign_flip(s173.scale(E - X @ coef))
    primary = {"n_donors": len(pts), "alpha": ALPHA, "p": p, "mean_scaled_E": m,
               "direction": "invasive > core" if m > 0 else "core > invasive",
               "coef_hypoxia": float(coef[0]), "coef_neuron": float(coef[1]),
               "coef_myeloid": float(coef[2]), "p_adjusted": p_adj, "mean_scaled_E_adjusted": m_adj,
               "n_genes": {"tam": o["sig"]["n_genes"], "hypoxia": o["hyp"]["n_genes"],
                           "neuron": on["sig"]["n_genes"], "myeloid": on["hyp"]["n_genes"]}}
    primary["step1_pass"] = bool(p < ALPHA and m < 0)
    primary["success"] = bool(primary["step1_pass"] and p_adj < ALPHA and m_adj < 0)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fp.write_text(json.dumps({"script": "185_trackA_tam_5ala_core_inv",
                              "label": "EXPLORATORY-PROVENANCE bulk RNA-seq; sealed; scored once",
                              "seed": SEED, "patients": pts, "primary": primary,
                              "E": {"tam": E.tolist(), "hypoxia": H.tolist(), "neuron": N.tolist(),
                                    "myeloid": M.tolist()}}, indent=2))
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
