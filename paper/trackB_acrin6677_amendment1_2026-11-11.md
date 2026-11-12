# ACRIN 6677 amendment 1 (2026-11-11)

No reviewed source established the stored unit of S0 `s0e33`. The [Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3688018/fullTextXML) failed. [NCBI efetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=3688018) returned no verified unit. Europe PMC search identified [Boxerman et al.](https://europepmc.org/articles/PMC6620020) and [Schmainda et al.](https://europepmc.org/articles/PMC4490871), whose XML endpoints failed. The [TCIA collection description](https://www.cancerimagingarchive.net/collection/acrin-dsc-mr-brain/) did not establish the stored unit.

The primary measure switches to `s0e23`, TOTAL TUMOR AREA, for documentation reasons before any measurement value was read. The archive data dictionary labels its components “TUMOR AREA (CM2)”. Use `c = 0.1 cm²` for `y = log(value + c)`. [Stage 0](trackB_acrin6677_stage0_2026-11-12.md) counted 83 patients with at least three measured area visits. The [earlier no-go](trackB_acrin6677_nogo_2026-11-11.md) applies only to `s0e33`, which is secondary and exploratory in arbitrary units.

Use patient `cn`, day `S0e2d`, reader `s0e4deid`, and visit label `s0e1`. A measured visit is a distinct day with at least one numeric value after applying the rules below. Average reader values on the same day, sort visits by day, and use the first three measured visits only. Require day1 < day2 < day3. Exclude uninterpretable rows (`s0e5 = 1`). For other rows, count missing area as zero when lesions are absent (`s0e13 = 1`).

The [protocol](trackB_acrin6677_protocol_2026-11-12.md) has one history increment for two model parameters. Set `g = 0` and fit `k` from that increment if `y2 < y1`. Otherwise set `k = 0` and fit `g`. Keep `tau = 56 days`. M equals B2 on flat or rising histories and differs from B2 only when the measurement fell. The frozen comparison uses all eligible patients against the best of B1, B2, and B3. The falling-history statement concerns M versus B2 only.

No patient measurement value, outcome, or treatment arm was read before this freeze. The affected analysis is exploratory under the protocol's amendment rule.
