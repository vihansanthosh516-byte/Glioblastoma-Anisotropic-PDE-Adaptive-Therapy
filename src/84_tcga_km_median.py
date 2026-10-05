"""Script 84: Kaplan-Meier median overall survival for the TCGA-GBM cohort used in scripts 36-37.

Script 36 reports np.median(time[event == 1]), the median among deaths only. That is not a
survival median because it ignores censored patients. This script reads the same cohort file
and reports the Kaplan-Meier median with a 95% CI, plus the old statistic for comparison.

Input:  output/clinical_mapped_cohort.csv (columns os_time_days, os_event)
Output: output/tcga_km_summary.json
"""
import json
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
IN_CSV = OUTPUT_DIR / "clinical_mapped_cohort.csv"
OUT_JSON = OUTPUT_DIR / "tcga_km_summary.json"


def main():
    df = pd.read_csv(IN_CSV)
    n_total = len(df)
    df = df.dropna(subset=["os_time_days", "os_event"])
    time = df["os_time_days"].to_numpy(float)
    event = df["os_event"].to_numpy(int)

    km = KaplanMeierFitter().fit(time, event)
    ci = median_survival_times(km.confidence_interval_)
    res = {
        "script": "84_tcga_km_median",
        "input": str(IN_CSV.relative_to(PROJECT_ROOT)),
        "n_rows_in_file": int(n_total),
        "n_analysed": int(len(df)),
        "n_events": int(event.sum()),
        "n_censored": int((event == 0).sum()),
        "km_median_days": float(km.median_survival_time_),
        "km_median_ci95_days": [float(ci.iloc[0, 0]), float(ci.iloc[0, 1])],
        "median_of_all_times_days_NOT_a_survival_median": float(np.median(time)),
        "median_of_death_times_only_days_script36_statistic": float(np.median(time[event == 1])),
    }
    OUT_JSON.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
