# Track A myeloid-origin analysis: result (2026-11-10, exploratory)

Protocol frozen and committed before scoring: `paper/trackA_myeloid_origin_protocol_2026-11-10.md` (commit 992fabf). Script `src/222_trackA_myeloid_origin.py` (run once from the worktree at 992fabf with the script not yet committed; committed unchanged afterwards). Output `output/trackA_myeloid_origin_20261110/result.json`. Exposed cohort (A5 used it); exploratory only. Not a drug-response measurement.

**Verdict under the frozen rule: INCONCLUSIVE.** Two of three positive conditions held in 4 of 5 donors; the contamination control condition held in only 3 of 5.

| Donor | Core / Periphery Myeloid cells | Program AUC | Matched-null 95th pct | NEO control AUC | Logistic AUC (LODO) |
|---|---|---|---|---|---|
| GS1 | 21 / 39 | 0.922 | 0.678 | 0.498 | 0.873 |
| GS2 | 71 / 57 | **0.305** | 0.616 | 0.474 | 0.368 |
| GS3 | 60 / 11 | 0.736 | 0.676 | 0.541 | 0.694 |
| GS5 | 23 / 11 | 0.830 | 0.731 | 0.731 | 0.933 |
| GS13 | 70 / 98 | 0.878 | 0.704 | **0.992** | 0.914 |
| Equal-donor mean | | 0.734 | 0.681 | 0.647 | 0.757 |

All 24 genes were found. Rule counts: program AUC > 0.5 in 4/5; above its matched null in 4/5; above NEO in 3/5 (needed 4). Program above logistic in 2/5 (descriptive).

## Reading, with limits

- In four donors, Myeloid cells from Core points carried a lower microglia-minus-macrophage score than Periphery cells, above size- and expression-matched random gene sets. This matches the published core-macrophage / periphery-microglia pattern.
- GS2 went the opposite way (0.305); the trained logistic baseline was also below 0.5 there (0.368), so GS2's labels separate poorly in both directions.
- In GS13 the neoplastic-marker score separated Core from Periphery almost perfectly inside "Myeloid" cells (0.992). That points to ambient tumor RNA, doublets, or annotation mixing in that donor's Core Myeloid cells. The program result in GS13 cannot be separated from that effect. GS5 also showed a NEO AUC of 0.731.
- The untrained 16-gene program was close to the trained logistic baseline on average (0.734 vs 0.757) and beat it in 2 donors.
- Weakest points: 11 Periphery Myeloid cells in GS3 and GS5; one exposed cohort; sampling-point labels, not cell origin; no ambient-RNA correction; gene lists are literature-derived and not matched to one paper's exact list.
- Novelty: low. The direction was published before; this is a donor-level check in another cohort. It does not meet the project's "new method" bar.
