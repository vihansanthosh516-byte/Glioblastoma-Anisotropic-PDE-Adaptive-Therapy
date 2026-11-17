"""Phase F: locked primary and declared secondary analyses. Run once, only with the locked config.

    py -3 src/recurrence_reliability/run.py --config configs/recurrence_reliability.json

Refuses an unlocked config, hash drift against qc_report.json, missing annotations, or a cohort that
differs from the QC-frozen cohort. Writes output/recurrence_reliability/<run_id>/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from recurrence_reliability import inference as inf  # noqa: E402
from recurrence_reliability import scoring as sc  # noqa: E402
from recurrence_reliability import tenx_io as rio  # noqa: E402
from recurrence_reliability.metadata import OUT, PROJECT_ROOT, load_annotation, sha256  # noqa: E402
from recurrence_reliability.qc import MAX_MITO_FRAC, MIN_GENES, nucleus_metrics  # noqa: E402

CODE_FILES = ["metadata.py", "tenx_io.py", "qc.py", "scoring.py", "inference.py", "run.py"]


class Specimen:
    """Filtered nuclei of one library with author labels; counts kept as genes x cells."""

    def __init__(self, library_id: str, annotation: dict[str, str]):
        lib = rio.load_library(library_id)
        umi, genes, mito = nucleus_metrics(lib)
        keep = (genes >= MIN_GENES) & (mito <= MAX_MITO_FRAC)
        lab = np.array([annotation.get(b, "NA") for b in lib.barcodes])
        # Amendment 1 (2026-11-16, before any scoring): nuclei without an author label are dropped,
        # never assigned a label; their count is reported. A library with > 5% unlabelled nuclei stops.
        unl = keep & (lab == "NA")
        self.n_unlabelled_dropped = int(unl.sum())
        if unl.sum() > 0.05 * max(1, keep.sum()):
            raise RuntimeError(f"{library_id}: more than 5% of filtered nuclei lack an author label")
        keep = keep & (lab != "NA")
        self.library_id = library_id
        self.counts = lib.counts[:, np.flatnonzero(keep)]
        self.cell_ids = np.array([f"{library_id}:{b}" for b in lib.barcodes[keep]])
        self.tumor = lab[keep] == "Tumor"
        self.symbols = lib.symbols
        self.gene_ids = lib.gene_ids

    def subset(self, malignant_only: bool):
        return self.counts[:, np.flatnonzero(self.tumor)] if malignant_only else self.counts


def module_index(symbols: np.ndarray, genes: list[str]) -> np.ndarray:
    pos: dict[str, list[int]] = {}
    for i, s in enumerate(symbols):
        pos.setdefault(s, []).append(i)
    return np.array([pos[g][0] for g in genes if len(pos.get(g, [])) == 1], dtype=int)


def pool_stats(specimens: list[Specimen], malignant_only: bool) -> tuple[np.ndarray, np.ndarray, int]:
    """Aggregate expression log2(mean CPM + 1) and per-gene mean log2(CPM/10+1) over all pooled nuclei."""
    n_genes = specimens[0].counts.shape[0]
    cpm_sum, log_sum, n = np.zeros(n_genes), np.zeros(n_genes), 0
    for s in specimens:
        m = s.subset(malignant_only).astype(np.float64)
        lib = np.asarray(m.sum(axis=0)).ravel()
        cpm_sum += np.asarray((m @ scipy.sparse.diags(1e6 / lib)).sum(axis=1)).ravel()
        log_sum += np.asarray(sc.log_cpm(m).sum(axis=1)).ravel()
        n += m.shape[1]
    return np.log2(cpm_sum / n + 1.0), log_sum / n, n


def score_cells(counts, center, modules_idx: dict, controls: dict, which: str) -> np.ndarray:
    if which == "R":
        le = sc.log_cpm(counts)
        per = {k: sc.tirosh_score(le, center, modules_idx[k], controls[k]) for k in modules_idx}
    else:
        per = {k: sc.ucell_score(counts, modules_idx[k]) for k in modules_idx}
    return sc.mes_like(per) if set(per) == {"MESlike1", "MESlike2"} else per[next(iter(per))]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg_path = (PROJECT_ROOT / args.config).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    if not cfg.get("locked"):
        raise SystemExit("config is not locked; biological scoring refused")
    qc = json.loads((OUT / "qc_report.json").read_text(encoding="utf-8"))
    if sha256(OUT / "qc_report.json") != cfg["qc_report_sha256"]:
        raise SystemExit("qc_report.json hash differs from locked config")
    for name, h in qc["file_sha256"].items():
        if sha256(rio.GSM_DIR / name) != h:
            raise SystemExit(f"hash drift: {name}")
    for k, h in cfg["module_sha256"].items():
        if sha256(PROJECT_ROOT / k) != h:
            raise SystemExit(f"module file hash drift: {k}")
    cohort = qc["final_cohort"]
    if [c["table1_pair"] for c in cohort] != cfg["cohort_table1_pairs"]:
        raise SystemExit("cohort differs from locked config")

    t0 = time.time()
    run_dir = OUT / cfg["run_id"]
    if run_dir.exists():
        raise SystemExit(f"{run_dir} exists; the locked analysis runs once")
    run_dir.mkdir(parents=True)
    manifest = json.loads((OUT / "metadata_manifest.json").read_text(encoding="utf-8"))
    ann = load_annotation()
    key_of = {l["library_id"]: l["annotation_key"] for s in manifest["specimens"] for l in s["libraries"]}
    lib_of = {s: qc["specimens"][s]["primary_library"] for s in qc["specimens"]}

    specs: dict[str, Specimen] = {}
    for c in cohort:
        for stage in ("Primary", "Recurrent"):
            sid = c[stage]
            specs[sid] = Specimen(lib_of[sid], ann[key_of[lib_of[sid]]])
    b2_specs = {sid: Specimen(lid, ann[key_of[lid]]) for sid, lid in cfg["independent_capture_libraries"].items()}
    rows, result = analyze(specs, b2_specs, cohort, rio.neftel_mes(), rio.verhaak_mes()["Mesenchymal"], cfg)
    result["unlabelled_nuclei_dropped"] = {s.library_id: s.n_unlabelled_dropped
                                           for s in [*specs.values(), *b2_specs.values()] if s.n_unlabelled_dropped}
    tsv = ["\t".join(rows[0].keys())] + ["\t".join(str(v) for v in r.values()) for r in rows]
    (run_dir / "patient_changes.tsv").write_text("\n".join(tsv) + "\n", encoding="utf-8")
    result = {
        "run_id": cfg["run_id"], "config_sha256": sha256(cfg_path),
        "code_sha256": {f: sha256(Path(__file__).parent / f) for f in CODE_FILES},
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                        "platform": platform.platform()},
        **result, "runtime_seconds": round(time.time() - t0, 1),
    }
    (run_dir / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    primary = result["primary_endpoint"]
    print(json.dumps({"n": len(rows), "opposite": primary["opposite_direction"], "ci": primary["opposite_ci95_clopper_pearson"],
                      "indeterminate": primary["indeterminate_coverage"], "reading": primary["prespecified_reading"]}))


def analyze(specs: dict, b2_specs: dict, cohort: list[dict], neftel_genes: dict[str, list[str]],
            verhaak_genes: list[str], cfg: dict) -> tuple[list[dict], dict]:
    """All declared analyses on loaded specimens. Pure function of its inputs and cfg seeds."""
    symbols = next(iter(specs.values())).symbols
    mods = {k: module_index(symbols, v) for k, v in neftel_genes.items()}
    verh = {"Verhaak_Mesenchymal": module_index(symbols, verhaak_genes)}
    seeds = cfg["seeds"]

    # Reference pools (fixed: all filtered nuclei of the frozen cohort's primary libraries)
    agg_m, cen_m, n_m = pool_stats(list(specs.values()), True)
    agg_w, cen_w, n_w = pool_stats(list(specs.values()), False)
    bins_m, bins_w = sc.expression_bins(agg_m), sc.expression_bins(agg_w)
    rng_c = np.random.default_rng(seeds["control_genes"])
    ctrl_m = {k: sc.control_genes(v, bins_m, rng_c) for k, v in mods.items()}
    ctrl_w = {k: sc.control_genes(v, bins_w, rng_c) for k, v in mods.items()}
    ctrl_v = {k: sc.control_genes(v, bins_m, rng_c) for k, v in verh.items()}

    cells: dict[str, dict[str, np.ndarray]] = {}
    for sid, s in specs.items():
        mal, allc = s.subset(True), s.subset(False)
        cells[sid] = {
            "R": score_cells(mal, cen_m, mods, ctrl_m, "R"),
            "A": score_cells(mal, cen_m, mods, ctrl_m, "A"),
            "R_whole": score_cells(allc, cen_w, mods, ctrl_w, "R"),
            "A_whole": score_cells(allc, cen_w, mods, ctrl_w, "A"),
            "R_verhaak": score_cells(mal, cen_m, verh, ctrl_v, "R"),
        }

    rng_b = np.random.default_rng(seeds["cell_bootstrap"])
    B = cfg["n_cell_bootstrap"]
    spec_stats = {sid: {k: {"mean": float(v.mean()), "se": inf.boot_se(v, rng_b, B), "n": int(len(v))}
                        for k, v in d.items()} for sid, d in cells.items()}

    rows = []
    for c in cohort:
        p, r = c["Primary"], c["Recurrent"]
        row = {"table1_pair": c["table1_pair"], "primary": p, "recurrent": r,
               "n_malignant_primary": spec_stats[p]["R"]["n"], "n_malignant_recurrent": spec_stats[r]["R"]["n"],
               "n_all_primary": spec_stats[p]["R_whole"]["n"], "n_all_recurrent": spec_stats[r]["R_whole"]["n"]}
        for k in ("R", "A", "R_whole", "A_whole", "R_verhaak"):
            d = spec_stats[r][k]["mean"] - spec_stats[p][k]["mean"]
            se = float(np.hypot(spec_stats[r][k]["se"], spec_stats[p][k]["se"]))
            row[f"delta_{k}"], row[f"se_{k}"], row[f"dir_{k}"] = d, se, inf.direction(d, se)
        rows.append(row)
    dR, dA = [x["dir_R"] for x in rows], [x["dir_A"] for x in rows]
    primary = inf.agreement_table(dR, dA)
    primary["indeterminate_coverage"] = {"R": sum(1 for d in dR if d == 0) / len(dR), "A": sum(1 for d in dA if d == 0) / len(dA)}
    primary["zero_failure_one_sided_upper95_at_this_n"] = inf.one_sided_upper_zero(len(dR))
    tol = cfg["tolerance_opposite_fraction"]
    lo, hi = primary["opposite_ci95_clopper_pearson"]
    primary["prespecified_reading"] = ("unreliable beyond tolerance (CI lower bound > tolerance)" if lo > tol else
                                       "reliable within tolerance (CI upper bound < tolerance)" if hi < tol else
                                       "inconclusive relative to tolerance")

    rng_p = np.random.default_rng(seeds["patient_bootstrap"])
    secondary = {
        "S0_continuous_patient_change_descriptive": {k: inf.patient_bootstrap_mean(np.array([x[f"delta_{k}"] for x in rows]), rng_p, cfg["n_patient_bootstrap"])
                                                     for k in ("R", "A", "R_whole", "A_whole", "R_verhaak")},
        "S1_whole_tissue_vs_malignant_R": inf.agreement_table(dR, [x["dir_R_whole"] for x in rows]),
        "S1_whole_tissue_vs_malignant_A": inf.agreement_table(dA, [x["dir_A_whole"] for x in rows]),
        "S4_construct_Verhaak_vs_Neftel_R": inf.agreement_table(dR, [x["dir_R_verhaak"] for x in rows]),
    }
    # S2 matched cell number
    rng_s = np.random.default_rng(seeds["matched_depth"])
    m_n, reps = cfg["matched_cell_number"], cfg["matched_reps"]
    s2 = []
    for x in rows:
        same_R = same_RA = 0
        for _ in range(reps):
            ds = {}
            for k in ("R", "A"):
                vp, vr = cells[x["primary"]][k], cells[x["recurrent"]][k]
                ds[k] = rng_s.choice(vr, m_n, replace=False).mean() - rng_s.choice(vp, m_n, replace=False).mean()
            same_R += np.sign(ds["R"]) == np.sign(x["delta_R"])
            same_RA += np.sign(ds["R"]) == np.sign(ds["A"])
        s2.append({"table1_pair": x["table1_pair"], "frac_reps_sign_matches_full_R": same_R / reps, "frac_reps_R_A_same_sign": same_RA / reps})
    secondary["S2_matched_cell_number"] = {"n_cells_per_specimen": m_n, "reps": reps, "per_patient": s2,
                                           "median_frac_sign_matches_full_R": float(np.median([y["frac_reps_sign_matches_full_R"] for y in s2]))}
    # S3 library-aware (independent captures; batch2 nuclei are not in the reference pools)
    s3 = []
    for sid, lid in cfg["independent_capture_libraries"].items():
        b2 = b2_specs[sid]
        x = next(r for r in rows if sid in (r["primary"], r["recurrent"]))
        other = x["recurrent"] if sid == x["primary"] else x["primary"]
        sign = 1 if sid == x["recurrent"] else -1
        out = {"specimen": sid, "library": lid, "table1_pair": x["table1_pair"], "n_malignant_batch2": int(b2.tumor.sum())}
        if b2.tumor.sum() < 2:
            out["note"] = "fewer than 2 malignant nuclei in the independent capture; not scored"
            s3.append(out)
            continue
        for k in ("R", "A"):
            v2 = score_cells(b2.subset(True), cen_m, mods, ctrl_m, k)
            pooled = np.concatenate([cells[sid][k], v2])
            for label, vals in (("batch2_only", v2), ("pooled_libraries", pooled)):
                d = sign * (vals.mean() - spec_stats[other][k]["mean"])
                se = float(np.hypot(inf.boot_se(vals, rng_b, B), spec_stats[other][k]["se"]))
                out[f"{k}_{label}"] = {"delta": float(d), "dir": inf.direction(d, se)}
            out[f"{k}_primary_library"] = {"delta": x[f"delta_{k}"], "dir": x[f"dir_{k}"]}
        s3.append(out)
    secondary["S3_library_aware"] = s3
    # S5 equal-patient vs pooled-cell (different estimands)
    allp = np.concatenate([cells[x["primary"]]["R"] for x in rows])
    allr = np.concatenate([cells[x["recurrent"]]["R"] for x in rows])
    secondary["S5_equal_patient_vs_pooled_cell_R"] = {"equal_patient_mean_delta": float(np.mean([x["delta_R"] for x in rows])),
                                                     "pooled_cell_difference": float(allr.mean() - allp.mean())}
    # S6 split-half calibration of the indeterminate rule (no recurrence information)
    rng_h = np.random.default_rng(seeds["split_half"])
    det = {"R": [], "A": []}
    for _ in range(cfg["split_half_reps"]):
        for sid in specs:
            for k in ("R", "A"):
                v = rng_h.permutation(cells[sid][k])
                h1, h2 = v[: len(v) // 2], v[len(v) // 2:]
                se = float(np.hypot(inf.boot_se(h1, rng_h, cfg["split_half_boot"]), inf.boot_se(h2, rng_h, cfg["split_half_boot"])))
                det[k].append(inf.direction(h1.mean() - h2.mean(), se) != 0)
    secondary["S6_split_half_false_determinate_rate"] = {k: float(np.mean(v)) for k, v in det.items()}
    secondary["S6_expected_under_calibration"] = 0.05

    return rows, {
        "n_patients": len(rows), "pool_nuclei": {"malignant": n_m, "whole_tissue": n_w},
        "module_genes_used": {k: int(len(v)) for k, v in {**mods, **verh}.items()},
        "control_genes": {k: int(len(v)) for k, v in ctrl_m.items()},
        "primary_endpoint": primary, "secondary": secondary, "specimen_stats": spec_stats,
    }


if __name__ == "__main__":
    main()
