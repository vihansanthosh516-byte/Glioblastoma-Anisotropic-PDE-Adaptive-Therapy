"""Geometry-only common support and eligible-patient count for LUMIERE (amendment G4). No tumor labels.

Inputs: registration_qc.json (rigid params/QC), transform_audit.json (native affines/shapes, ct1.tfm), and the
visit-2 atlas brain_mask (non-tumor) as forecast domain. Applies the frozen 25-patient minimum.
Usage: python -m src.fk_lumiere_eligibility
"""
from __future__ import annotations

import gzip
import hashlib
import json
import zipfile

import nibabel as nib
import numpy as np

from src.fk_lumiere_coverage_chain import domain_on_scoring_grid, visit_coverage
from src.fk_lumiere_registration_run import ATLAS_MASK, OUT_DIR as REG_DIR
from src.fk_lumiere_runner import MIN_ELIGIBLE_PATIENTS
from src.fk_lumiere_transform_audit import OUT_DIR as AUDIT_DIR, PROJECT_ROOT, ZIP_PATH, _visit_dirs, tfm_to_ras_4x4
from src.track_b_coverage_scoring import common_support

OUT_DIR = PROJECT_ROOT / "output" / "track_b_lumiere_eligibility_20261110"
VISITS = ("visit1", "visit2", "target")


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def patient_support(audit_rows, reg, atlas_mask_1mm):
    covs = []
    for k, name in enumerate(VISITS):
        v = audit_rows[k]
        leaves = [b for b in v["tfm"]["blocks"] if not b["type"].startswith("Composite")]
        params = None if k == 1 else reg["visit1" if k == 0 else "visit3"]["params"]
        covs.append(visit_coverage(name, native_shape=v["native"]["shape"], native_affine=np.array(v["native"]["affine"]),
                                   tfm_native_to_atlas_ras=tfm_to_ras_4x4(leaves),
                                   atlas_affine=np.array(v["atlas"]["affine"]), atlas_shape=tuple(v["atlas"]["shape"]),
                                   rigid_params=params))
    sup = common_support(covs, VISITS).mask
    dom = domain_on_scoring_grid(atlas_mask_1mm)
    return {"domain_voxels": int(dom.sum()), "support_voxels": int(sup.sum()),
            "domain_covered_voxels": int((dom & sup).sum()),
            "domain_covered_fraction": float((dom & sup).sum() / dom.sum()) if dom.any() else None,
            "per_visit_domain_covered_fraction": {c.visit: float((dom & c.mask).sum() / dom.sum()) for c in covs},
            "full_domain_coverage": bool(dom.any() and not (dom & ~sup).any())}


def main():
    reg_path, audit_path = REG_DIR / "registration_qc.json", AUDIT_DIR / "transform_audit.json"
    reg = {c["patient"]: c for c in json.loads(reg_path.read_text())["cases"]}
    audit = json.loads(audit_path.read_text())["patients"]
    rows = {}
    with zipfile.ZipFile(ZIP_PATH) as zf:
        names = zf.namelist()
        for pid, r in reg.items():
            row = {"registration_passed": r["passed"]}
            if r["passed"]:
                d = _visit_dirs(names, pid)[(audit[pid][1]["week"], audit[pid][1]["order"])]
                img = nib.Nifti1Image.from_bytes(gzip.decompress(zf.read(d + ATLAS_MASK)))
                try:
                    row.update(patient_support(audit[pid], r, np.asarray(img.dataobj) > 0))
                except ValueError as exc:      # empty/invalid support or domain: recorded failure, not coverage
                    row.update({"support_failure": str(exc), "full_domain_coverage": False})
                if row.get("domain_voxels") == 0:
                    row.update({"support_failure": "empty forecast domain", "full_domain_coverage": False})
            row["eligible"] = bool(r["passed"] and row.get("full_domain_coverage"))
            rows[pid] = row
    eligible = sorted(p for p, r in rows.items() if r["eligible"])
    n = len(eligible)
    res = {"rule": "paper/trackB_lumiere_geometry_amendment_2026-11-10.md G4",
           "read": "visit-2 atlas brain_mask only, plus saved transforms; no tumor labels",
           "inputs_sha256": {str(reg_path.relative_to(PROJECT_ROOT)): _sha(reg_path),
                             str(audit_path.relative_to(PROJECT_ROOT)): _sha(audit_path)},
           "n_patients": len(rows), "n_registration_failed": sum(not r["registration_passed"] for r in rows.values()),
           "n_support_failure": sum("support_failure" in r for r in rows.values()),
           "n_partial_coverage": sum(r["registration_passed"] and not r["eligible"] and "support_failure" not in r
                                     for r in rows.values()),
           "n_eligible": n, "minimum": MIN_ELIGIBLE_PATIENTS,
           "decision": ("proceed: eligible list frozen before any target read" if n >= MIN_ELIGIBLE_PATIENTS
                        else f"stop: {n} eligible < {MIN_ELIGIBLE_PATIENTS}; no scoring"),
           "patients": rows}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "eligibility.json"
    out.write_text(json.dumps(res, indent=1))
    if n >= MIN_ELIGIBLE_PATIENTS:
        frozen = {"eligible": eligible, "transforms": {p: {"visit1": reg[p]["visit1"]["params"],
                                                           "visit3": reg[p]["visit3"]["params"]} for p in eligible},
                  "eligibility_sha256": _sha(out), **res["inputs_sha256"]}
        fz = OUT_DIR / "frozen_cohort.json"
        fz.write_text(json.dumps(frozen, indent=1, sort_keys=True))
        (OUT_DIR / "frozen_cohort.sha256").write_text(_sha(fz) + "  frozen_cohort.json\n")
    print(json.dumps({k: v for k, v in res.items() if k != "patients"}, indent=1))


if __name__ == "__main__":
    main()
