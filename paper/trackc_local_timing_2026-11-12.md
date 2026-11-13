# Track C v2: local timing of one development world (2026-11-12)

One development-pilot world (`dev_pilot-single-00000`, primary cell, 2 patients per cohort, 135 controller models, 500 particles, 1 worker) ran on the Windows PC with the real benchmark. Source: `output/trackc_v2_local_timing/dev_pilot/primary/m02/manifest.json` and its checkpoint (`elapsed_seconds`).

- Time: **4,079 seconds (1.13 hours)** for one world.
- Two other single-core jobs (precision screen, interval cross-check) ran on the PC at the same time, so this may be a little slow.
- Earlier record: 0.793 hours per world at 2 patients and 81 models (`paper/gap_status_2026-11-09.md` C6). Ratio 1.43, below the 1.67 model-count ratio.
- Only the elapsed time and the names of the contrasts were read. No policy utility was examined. This world is development data, never evaluation data.

## Rough Spark estimate (assumption: one Spark core about as fast as one PC core; not measured)

At m = 2 and the screen floor of 276 worlds per family (552 worlds): 552 x 1.13 = about 624 core-hours. With 18 workers: **about 35 hours** for the primary cell. Each of the 6 sensitivity cells costs about the same (scan28 more, scan112 a little more because of the longer forecast). All 7 cells: about 10 days. The pilot at m = 12 (40 worlds) and the timing run add about 1 to 2 days.

These numbers are estimates. Spark runbook step 2 measures the real ones.
