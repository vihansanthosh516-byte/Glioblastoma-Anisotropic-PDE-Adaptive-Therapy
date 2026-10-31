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

## 3. Further metadata checks (same day)
- **E-MTAB-8743 raw size (ENA filereport API, study PRJEB36495):** core 10 runs 27.24 GB; invasive 10 runs 26.59 GB; rim 10 runs 27.15 GB; whole tumour 20 runs 47.43 GB; total 128.41 GB. All 10 individuals have both a core and an invasive sample. Rebuilding counts for core + invasive needs 53.83 GB of FASTQ. Not downloaded.
- **E-MTAB-8743 prior project use:** no file in `src/`, `paper/` or `START_HERE.md` (other than this addendum) names E-MTAB-8743, the 5ALA repository or 5-ALA. No prior exposure found.
- **Zenodo 20709378 author metadata table** (`Figure_2/spatial_data_metadata.txt`, 781 bytes): 25 sections. Six sample-name prefixes have both an explicit "core" and an explicit "edge" label: B185, B186, ST073021, WU1220, WU1221, WU1227 (Codex counted five; B185 core and edge sections both have damaged or low-res H&E). Labels are per section, some are hedged ("margin, probably core"), many H&E entries are "damaged". Prefix = donor is not verified. Six is below the 10-donor gate.
- **CFB GTV vs MU labels 1+3:** existing notes (`research_notes/GBM CFB scoring readiness/`) already record that CFB GTV is a binary, mixed-provenance contour (clinical at some visits, nnU-Net-completed at others; only 44 generated masks reviewed). It is not equivalent to MU labels 1+3. Status FAIL for direct substitution; a CFB-specific endpoint would be a new protocol.

## 4. Source table (PASS / FAIL / UNKNOWN)
| Item | E-MTAB-8743 (Track A) | Zenodo 20709378 (Track A) | CFB-GBM v4 (Track B) | MU-Glioma-Post (Track B) |
|---|---|---|---|---|
| Open download, no request | PASS | PASS (Zenodo CC BY 4.0) | PASS (CC BY 4.0, local copy) | PASS (CC BY 4.0) |
| Reuse terms for exact files | UNKNOWN (matrix: MIT code licence only) | PASS for Zenodo files; GEO separate | PASS | PASS |
| Label independent of outcome | PASS (5-ALA at surgery) | UNKNOWN (per section; some hedged) | FAIL vs MU target (binary, mixed provenance) | PASS (source labels) |
| Donor / patient count | PASS (10 paired) | FAIL (6 prefixes, unverified) | UNKNOWN (38 metadata candidates, 0/107 preflight pass) | UNKNOWN (18 candidates) |
| Data provenance | UNKNOWN (matrix pipeline undocumented) | UNKNOWN (spot crosswalk) | PASS (hash difference explained) | PASS |
| Quality / registration | n/a (bulk) | FAIL in part (damaged H&E) | UNKNOWN (no registration QC) | UNKNOWN (registration, label-zero support) |
| Prior project exposure | PASS (none found) | PASS (none found) | UNKNOWN (screened before) | FAIL (examined before; exploratory only) |
| Ready for a sealed test | NO | NO | NO | NO |

**Closest to ready:** E-MTAB-8743 for a bulk-tissue Track A question. Two routes: (a) use the authors' 2.6 MB matrix and label the provenance gap; or (b) rebuild counts from 53.83 GB of public FASTQ (needs the user's OK). Track B has no dataset that passes; the remaining gates need image arrays and registration QC, which this audit does not open.
