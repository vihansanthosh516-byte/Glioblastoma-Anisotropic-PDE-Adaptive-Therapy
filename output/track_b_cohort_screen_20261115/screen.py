"""Track B untouched-cohort screen. Input-side counts only.

Reads metadata tables only (scan days, file availability, clinical columns that exist before the outcome scan).
Reads NO mask, NO image, NO tumor volume at any outcome scan, NO radiomics table.
Eligibility rule of the locked question: consecutive scans, interval 14 to 365 days. The "input core > 0" rule
needs the masks and is NOT applied here, so every count is an upper limit.
"""
import csv, collections, json, statistics as st
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
WAW = ROOT / "data" / "external" / "waw"
RHUH = ROOT / "data" / "external" / "rhuh" / "clinical_data_TCIA_RHUH-GBM.csv"
CFB = ROOT / "research_notes" / "GBM patient analysis readiness 3" / "source_tables"


def waw():
    print("== WAW (Zenodo 10.5281/zenodo.22205631)")
    meta = json.load(open(ROOT / "data" / "external" / "waw.json", encoding="utf-8"))
    print("   license", meta["metadata"].get("license"), "| access", meta["metadata"].get("access_right"))
    print("   mask archives:", [(f["key"], f["size"]) for f in meta["files"] if "mask" in f["key"]])
    r = list(csv.DictReader(open(WAW / "2_2_postop_scans_coreg_meta.csv", encoding="utf-8-sig")))
    ex = collections.defaultdict(dict)
    for x in r:
        p = x["file_name"].split("_")
        ex[p[0]][p[1]] = float(x["days_postop"]) if x["days_postop"] not in ("", "nan") else None
    print("   patients with a post-operative exam", len(ex), "| exams per patient", sorted(collections.Counter(len(v) for v in ex.values()).items()))
    cl = list(csv.DictReader(open(WAW / "1_0_clinical_molecular_data.csv", encoding="utf-8-sig"), delimiter=";"))
    grade = {x["id"]: x["who_grade"].strip() for x in cl}
    pairs = []
    for pid, v in ex.items():
        d = sorted(t for t in v.values() if t is not None)
        pairs += [(pid, b - a) for a, b in zip(d[:-1], d[1:])]
    ok = [p for p in pairs if 14 <= p[1] <= 365]
    print("   consecutive post-operative pairs", len(pairs), "| patients", len({p[0] for p in pairs}))
    print("   pairs with interval 14-365 d", len(ok), "| patients", len({p[0] for p in ok}), "| grade 4 patients", len({p[0] for p in ok if grade.get(p[0]) == "4"}))
    print("   interval days min/median/max", min(p[1] for p in ok), st.median(p[1] for p in ok), max(p[1] for p in ok))


def rhuh():
    print("== RHUH-GBM (early post-operative scan -> recurrence scan)")
    r = list(csv.reader(open(RHUH, encoding="utf-8-sig")))[1:]
    f = lambda s: float(s) if s.strip() not in ("", "NA", "nan") else None
    res, pfs = [f(x[11]) for x in r], [f(x[20]) for x in r]
    print("   patients", len(r), "| residual enhancing tumor > 0 at the early post-operative scan", sum(1 for v in res if v and v > 0),
          "| progression-free time 14-365 d (proxy for the scan interval)", sum(1 for v in pfs if v and 14 <= v <= 365),
          "| both", sum(1 for a, b in zip(res, pfs) if a and a > 0 and b and 14 <= b <= 365))


def cfb():
    print("== CFB-GBM v4 (TCIA), availability tables only")
    m = pd.read_csv(CFB / "CFB-GBM_mri_availability_v02_20260129.tsv", sep="\t")
    g = pd.read_csv(CFB / "CFB-GBM_treatment_imaging_availability_v03_20260619.tsv", sep="\t")
    print("   GTV files by timepoint and type:", g.groupby(["temporality", "gtv_type"]).gtv.sum().to_dict())
    ok = m[(m.t1gd == 1) & (m.flair == 1)]
    print("   scans with T1Gd and FLAIR:", ok.groupby("temporality").size().to_dict(), "| of those with pre-contrast T1:", int(((ok.t1eg == 1) | (ok.t1tse == 1)).sum()), "| with T2:", int((ok.t2tse == 1).sum()), "| total", len(ok))
    p = m.pivot(index="id_patient", columns="temporality", values="time_diff_t0 (weeks)")
    has = {(r.id_patient, r.temporality): (r.t1gd == 1 and r.flair == 1) for r in m.itertuples()}
    elig = set()
    for a, b in [("t0", "t1"), ("t1", "t2")]:
        x = ((p[b] - p[a]) * 7.0).dropna()
        x = x.loc[[i for i in x.index if has.get((i, a), False) and has.get((i, b), False)]]
        e = x[(x >= 14) & (x <= 365)]; elig |= set(e.index)
        print(f"   {a}->{b}: pairs with T1Gd+FLAIR at both scans {len(x)} | interval 14-365 d {len(e)} | interval days at 0/10/25/50/75/90/100%:", e.quantile([0, .1, .25, .5, .75, .9, 1]).round(0).tolist())
    print("   patients with at least one eligible pair:", len(elig))
    return len(elig)


def power(n, effect=0.027, sd=0.154):
    from scipy.stats import norm
    return float(norm.cdf(effect * np.sqrt(n) / sd - 1.96))


if __name__ == "__main__":
    waw(); rhuh(); n_cfb = cfb()
    print("== power against persistence (assumes effect 0.027 from loco4_all_baselines_stdout.txt and SD 0.154 from the LUMIERE result.json; two-sided 5%)")
    for label, n in [("WAW only", 28), ("CFB only", n_cfb), ("CFB + WAW", n_cfb + 28), ("needed for 80%", 255)]:
        print(f"   {label}: n={n}, power={power(n):.2f}")
