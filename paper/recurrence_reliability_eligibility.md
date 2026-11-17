# Recurrence reliability: metadata, annotation, and batch eligibility (Phase B)

Date: 2026-11-16. Metadata only. No expression matrix, gene list, program score, or survival value was opened. Machine-readable record: `output/recurrence_reliability/metadata_manifest.json`; table: `output/recurrence_reliability/candidate_table.tsv`; code: `src/recurrence_reliability/metadata.py`.

## Result

**23 of the 27 conservative candidate pairs (46 specimens) pass the metadata, annotation, and batch gates.** Four pairs are excluded by rules fixed in code before counts were computed. No pair was added: the eight Table 1 snRNA specimens without an exact GEO title match (for example SF10022, SF10127, SF4400, SF6996) have author annotations but no GEO matrix under that ID, so adding them would need a guessed join.

| Table 1 pair | Specimens (P / R) | Exclusion reason |
|---|---|---|
| 15 | SF3996 / SF4449 | SF4449 has 22 author-labelled tumor nuclei (< 50) |
| 23 | SF7307 / SF7388 | SF7307 primary library: 197 of 223 barcodes annotated (join 0.88 < 0.95) |
| 28 | SF9871 / SF10108 | SF9871 has 6 author-labelled tumor nuclei (< 50) |
| 35 | SF11587 / SF11981 | SF11981 has no author annotation sample |

## Gate states

| Gate | State | Evidence |
|---|---|---|
| Specimen identity and pairing | Pass (conservative) | 27 pairs from the Table 1 crosswalk with unique exact SF-ID GSM matches and agreeing stage labels (`recurrence_specimen_crosswalk_2026-11-16.md`) |
| Recurrence order and treatment | Pass, author-described | Wang et al. Results: "For 36 patients we profiled paired specimens from the primary untreated tumor and matched first recurrence ... all patients had been treated only with standard-of-care therapy (TMZ, IR and surgical resection)." Table 1 has no per-row first-recurrence field, so first recurrence is taken from the authors' cohort description, not verified per patient. Pair 36 (two recurrences) is not a candidate. |
| Malignant annotation | Pass for 23 pairs | GEO `GSE174554_Tumor_normal_metadata.txt.gz` (1,422,162 bytes, SHA-256 `5136964f...65ac6`): 254,288 rows, 100 annotation sample keys, labels Tumor/Normal only, 8,347 exact duplicate rows with identical labels, no conflicting labels. Wang Methods: labels follow CONICSmat v1.0 CNV calls (LRT < 0.05, BIC difference > 50, posterior > 0.5). These are operational CNV-based labels, independent of any mesenchymal score, not perfect ground truth. |
| Annotation join | Pass for 23 pairs | Join by specimen and 16-base barcode after removing the GEM suffix (`-1` in barcodes.tsv, `.1` in the annotation). The matching key must be unique among the specimen ID and specimen ID + `v2`, its barcodes must be a subset of the library's barcodes, and it must cover at least 95% of them. Eight candidates match under the `v2` key (SF3996, SF7025, SF7388, SF9372, SF4324, SF11916, SF11815, SF10857); every one is a full set match. Join rate is 1.0 for every eligible primary library except SF4810 (0.9965) and SF10433 (0.9991). |
| Library relationships | Resolved | See below. |
| Access and terms | Pass for local analysis | Processed files are public on GEO (no login). GEO states that NCBI places no use restrictions but submitters may hold rights. The Wang paper is open access. Raw reads are at EGA EGAS00001004909 (controlled; not used). Plan: keep matrices local and uncommitted; publish only code and aggregate outputs. |
| Minimum malignant support | Pass for 23 pairs | Floor of 50 author-labelled tumor nuclei per specimen, fixed in code before counts were computed. The smallest eligible counts are SF6098 (55) and SF12115 (59). Phase C re-checks this after technical filters. |

## Batch and library interpretation

Seven candidate specimens have a second (`batch2`) library. Barcode overlap with the primary library decides the relationship. A 10x v3 barcode whitelist has millions of barcodes, so large overlap cannot be chance.

- **Same capture (barcode subset), never added:** SF10099 (117 of 117 batch2 barcodes in the primary library), SF4449 (100%), SF4810 (92.7%), SF7307 (100%). These are re-sequencing or re-processing of the same nuclei, not new material.
- **Independent capture:** SF12243 (overlap 0.15%), SF10433 (0.17%), SF4209 (0%). These are separate captures of the same specimen. They are used only in a library-aware secondary sensitivity analysis, with library provenance kept in cell IDs. They never add patients.
- **Primary-library rule:** each specimen uses the library whose file prefix is `<GSM>_<specimen>` (no batch suffix). The rule was fixed before any expression was read.

## Incidental exposure

Computing the annotation join exposed per-specimen counts of author-labelled tumor and normal nuclei (composition). These counts relate to the secondary whole-tissue vs malignant comparison. They are not the program score. No mesenchymal or other program values were seen.

## Remaining uncertainty

- First recurrence is author-described, not verified per patient.
- CNV-based labels can miss tumor nuclei without detectable CNVs and can mislabel normal nuclei; the error rate is unknown.
- The GEO barcodes are already filtered by the authors (Cell Ranger 3.0.2, pre-mRNA reference, > 2.5% mitochondrial reads and < 200 genes removed, DoubletFinder). Phase C must verify that the matrices are integer counts.
- 23 patients give low precision. Even zero failures gives an exact one-sided 95% upper bound of `1 - 0.05**(1/23)` = 0.122 (computed in Phase D).
