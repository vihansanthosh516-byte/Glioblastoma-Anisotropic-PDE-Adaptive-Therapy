# Track A HGCC bridge: metadata gate check (2026-11-12)

Protocol: `paper/trackA_hgcc_tmz_response_bridge_2026-11-12.md` (committed in 1d112ce). This check used file headers, sample IDs, gene names, and TMZ **presence** only. **No TMZ response magnitude was printed, summarized, or analyzed.** No expression value was read.

## Files (official HGCC portal, https://hgcc.se/portal/, downloaded 2026-11-12, stored in git-ignored `data/hgcc/`)

| File | Bytes | SHA-256 |
|---|---|---|
| HGCC_RNA_expression.txt | 46,628,606 | 31199414be5763c4f665a28f8b9d79ade923a704d808171c765dd4209b952e08 |
| HGCC_drug_response_AUC.txt | 369,637 | 1e6ce054239896926d7157e12d4ecfe8b702cbaf6705467bb70a28a6045bf54d |
| HGCC_clinical.txt | 2,573 | d0ec5b7bf94f8cca6164a771c69669d64f1de30d228fd7bb23d8e8c1fe19e7fc |

## Gates

| Gate | Result | Evidence |
|---|---|---|
| Release and hashes | Pass | Table above |
| Exact join, at least 25 origins | **Pass** | Drug file: 261 drugs x 110 columns; one row named "Temozolomide", present in 110/110 columns. 85 columns have an exact `U####` ID; 21 suffix columns (`b`, `replicated`, `U3065_<n>` sublines) are excluded rather than prefix-matched. Expression: 23,833 genes, 114 profiles, 101 distinct IDs (11 IDs have two profiles, to be averaged per protocol). Exact ID joins: **83**, matching Ju et al.'s 83. The clinical table lists one row per `U####MG` ID with its own sex and age, which supports one patient origin per ID; HGCC does not publish an explicit origin crosswalk, so this stays an assumption. |
| No overlap with Ntafoulis score-source cultures | Pass (by institution) | HGCC lines are Uppsala `U####MG`; Ntafoulis cultures are Erasmus MC `GS.###`. Different biobanks and IDs. Not proven at patient level. |
| Exact unique mapping of all 119 frozen genes | **FAIL** | 102 of 119 map exactly; none duplicated. Missing (17): AMY2B, CFAP44, EXTL3-AS1, KANTR, LENG8-AS1, LINC00907, LINCR-0001, MRNIP, NBPF11, PRKCQ-AS1, PSMA3-AS1, PTOV1-AS1, PTOV1-AS2, RSRP1, SLC25A25-AS1, TWF2, ZNF547. Most are long non-coding RNAs or renamed symbols absent from the HTA 2.0 gene-level annotation. |
| Data-use terms | Recorded, unclear | Portal and resource page describe open data; no explicit licence for the response table was found. Only aggregate statistics may be shared. |

## Decision

By the frozen protocol (step 2: "Require exact, unique mapping of all 119 frozen genes; otherwise stop"), **the analysis stops here.** No outcome was accessed.

## Possible amendment (not adopted; needs the user's choice)

A dated amendment made **before any TMZ value is read** could define a 102-of-119 exact-mapped score, following the project's own 119-of-121 precedent: no aliases, no substitutions, original directions, matched controls rebuilt with the same direction counts. It would be labelled as a second, smaller proxy. Weakness: 14% of the panel is lost, and the lost genes are mostly non-coding, so the score is a different, narrower signature. Alias mapping (for example old symbols for renamed genes) is not allowed under the project rules.

Track A's A5 failure, Track B's outcomes, and the Track C stop are unchanged.
