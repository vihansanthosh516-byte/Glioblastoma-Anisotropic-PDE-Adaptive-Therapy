# Addendum to the open-data readiness audit (2026-10-30, Claude)

Continues the two open checks at the end of Codex's `paper/open_data_readiness_2026-10-30.md` (commit 4d8b112, branch `codex/duke-access-request-confirmed`). Metadata only. No MRI, mask or expression arrays were opened. Nothing was downloaded except the small public SDRF text file (32.8 KB, read in a scratch folder).

## 1. CFB crosswalk hash mismatch: RESOLVED (explained, not a data change)
- Current file `research_notes/GBM CFB scoring readiness/CFB_v4_candidate_crosswalk_20261028.csv`: SHA-256 `05a4ead1...` (164,573 bytes, 2026-10-28 18:54).
- Backup `CFB_v4_candidate_crosswalk_20261028.csv.pre_t0_reference.bak`: SHA-256 `b1106bb2...` (164,681 bytes, 17:46). This is the hash in the older preflight summaries (`output/cfb_origin_preflight_20261028/`, `_hardened/`, `_verified/`, `_final/`).
- Difference: same 107 rows; 0 changed values in shared columns. Only the registration-reference columns were renamed from t1-reference (`t0_to_t1_*`, `t2_to_t1_*`) to t0-reference (`t1_to_t0_*`, `t2_to_t0_*`).
- The current hash is bound by `output/cfb_origin_preflight_20261028_t0_reference/summary.json` (19:01).
- All five preflight runs: 0 pass, 107 fail, `scoring_authorized: false`. So no preflight validates any CFB patient for scoring, under either version.

## 2. E-MTAB-8743 count-matrix provenance: PARTLY RESOLVED
- ArrayExpress record (BioStudies API): released 2020-12-15; 50 assays, 10 individuals x 5 samples; sampling site: whole tumour 20, tumour core 10, tumour rim 10, invasive region 10; phenotype 5ALA-positive 10, 5ALA-negative 10. Disease: glioblastoma multiforme (all 50).
- The SDRF lists no processed (derived) data file. ArrayExpress holds raw reads only.
- The authors' `Data/Final_count_exon.csv` (2,580,293 bytes) was uploaded 2022-07-18 with the commit message "Upload data files". The repository README and `Data/README.md` do not say how the matrix was made (aligner, reference, or which runs). Provenance of the matrix: **UNKNOWN**.
- Region labels: the invasive samples were picked by 5-ALA fluorescence during surgery (SDRF phenotype column), not by expression. This is an independent label, which is good for the Track A question.
- Open-access rule (your definition: anyone can download, no email): the raw reads and the SDRF are public with no request step. The GitHub matrix is public too. Reuse terms: the repository has only an MIT licence (written for code); no separate data licence was found. The ArrayExpress API record shows no licence field.

## What this means
- CFB: the hash question is closed. CFB is still not scoring-ready (0/107 preflight pass).
- E-MTAB-8743: the strongest open Track A lead. Two gaps remain:
  1. Matrix provenance is undocumented. A clean fix is to rebuild counts from the public raw reads with a stated pipeline. This needs a large download and your OK (size not yet checked).
  2. It is bulk RNA-seq (10 donors), not single cells. The question must be stated for bulk samples.
