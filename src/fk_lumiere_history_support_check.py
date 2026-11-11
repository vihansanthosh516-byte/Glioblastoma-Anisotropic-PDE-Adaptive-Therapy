"""Pre-target history-mask support check (amendment G4). Reads visit-1 and visit-2 atlas seg_mask ONLY.

For each geometry-qualified patient: frozen forecast domain = origin_support(visit-2 brain mask, history labels),
on the 2 mm grid (G2); visit-1 labels are first aligned to visit 2 with the saved rigid transform (nearest
neighbour). The patient stays eligible only if common support covers that whole domain. Target (visit-3)
masks are never opened. Usage: python -m src.fk_lumiere_history_support_check
"""
from __future__ import annotations

import gzip
import hashlib
import json
import zipfile

import nibabel as nib
import numpy as np

from src.fk_longitudinal_prep import apply_rigid, origin_support
from src.fk_lumiere_coverage_chain import domain_on_scoring_grid
from src.fk_lumiere_eligibility import OUT_DIR, VISITS, visit_coverages
from src.fk_lumiere_registration_run import ATLAS_MASK, OUT_DIR as REG_DIR
from src.fk_lumiere_runner import MIN_ELIGIBLE_PATIENTS
from src.fk_lumiere_transform_audit import OUT_DIR as AUDIT_DIR, ZIP_PATH, _visit_dirs
from src.track_b_coverage_scoring import common_support

SEG = "DeepBraTumIA-segmentation/atlas/segmentation/seg_mask.nii.gz"


def _img(zf, name):
    raw = zf.read(name)
    return np.asarray(nib.Nifti1Image.from_bytes(gzip.decompress(raw)).dataobj), hashlib.sha256(raw).hexdigest()


def main():
    fz = json.loads((OUT_DIR / "frozen_cohort.json").read_text())
    reg = {c["patient"]: c for c in json.loads((REG_DIR / "registration_qc.json").read_text())["cases"]}
    audit = json.loads((AUDIT_DIR / "transform_audit.json").read_text())["patients"]
    rows, read = {}, []
    with zipfile.ZipFile(ZIP_PATH) as zf:
        names = zf.namelist()
        for pid in fz["eligible"]:
            d = _visit_dirs(names, pid)
            p1 = d[(audit[pid][0]["week"], audit[pid][0]["order"])]
            p2 = d[(audit[pid][1]["week"], audit[pid][1]["order"])]
            h1, s1 = _img(zf, p1 + SEG)
            h2, s2 = _img(zf, p2 + SEG)
            bm, s3 = _img(zf, p2 + ATLAS_MASK)
            read += [(p1 + SEG, s1), (p2 + SEG, s2), (p2 + ATLAS_MASK, s3)]
            h1_on_v2 = apply_rigid((h1 > 0).astype(float), reg[pid]["visit1"]["params"], spacing_mm=(1, 1, 1), order=0) > 0
            dom1 = origin_support(bm > 0, [h1_on_v2, h2 > 0])
            dom = domain_on_scoring_grid(dom1)
            sup = common_support(visit_coverages(audit[pid], reg[pid]), VISITS).mask
            outside = int((dom & ~sup).sum())
            rows[pid] = {"domain_voxels_with_history": int(dom.sum()),
                         "history_voxels_outside_brain_mask_2mm": int((dom & ~domain_on_scoring_grid(bm > 0)).sum()),
                         "domain_voxels_outside_support": outside, "passes": outside == 0}
    assert all(n.endswith(SEG) or n.endswith(ATLAS_MASK) for n, _ in read)
    assert not any("week" in n and n.endswith(SEG) and n.startswith(tuple(
        _visit_dirs(names, pid)[(audit[pid][2]["week"], audit[pid][2]["order"])] for pid in fz["eligible"])) for n, _ in read), \
        "a target seg_mask was read"
    final = sorted(p for p, r in rows.items() if r["passes"])
    res = {"rule": "amendment G4 pre-target history-mask support check",
           "read_scope": "visit-1 and visit-2 atlas seg_mask and visit-2 atlas brain_mask of geometry-qualified patients; "
                         "no visit-3 (target) mask",
           "files_read": {n: h for n, h in read}, "n_input": len(rows), "n_pass": len(final),
           "minimum": MIN_ELIGIBLE_PATIENTS,
           "decision": ("proceed" if len(final) >= MIN_ELIGIBLE_PATIENTS else
                        f"stop: {len(final)} eligible < {MIN_ELIGIBLE_PATIENTS}; no scoring"),
           "eligible_final": final, "patients": rows}
    out = OUT_DIR / "history_support_check.json"
    out.write_text(json.dumps(res, indent=1, sort_keys=True))
    (OUT_DIR / "history_support_check.sha256").write_text(
        hashlib.sha256(out.read_bytes()).hexdigest() + "  history_support_check.json\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("patients", "files_read")}, indent=1))


if __name__ == "__main__":
    main()
