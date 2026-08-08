"""Pipeline step 34 — how overdispersed are these counts, measured rather than asserted.

The paper's account of why the observation model matters rests on one empirical
step: weekly dengue counts are far more variable than a Poisson process allows,
so a Poisson likelihood treats ordinary noise as a large residual and rewards any
model flexible enough to chase it. Climate covariates are exactly that flexibility.

That step had been asserted. Two versions of the assertion appeared in earlier
drafts --- "far larger" and "two orders of magnitude larger" --- and the second was
wrong: it is true only of the extreme tail.

The measure here is deliberately model-free. Residuals are taken about a five-week
centred moving average rather than about a fitted transmission model, so the
number does not depend on any of the structures this paper compares being correct.
A Poisson process gives an index of one.

Writes ``41_dispersion.csv``, one row per window.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SMOOTH_WEEKS = 5
MIN_WEEKS = 20
MIN_USABLE = 10

cfg = load_config()
raw = resolve(cfg, "raw")
tables = resolve(cfg, "tables")


def dispersion_index(counts: np.ndarray) -> float:
    """Mean of residual squared over the local mean, about a moving average.

    A Poisson process gives one. Values above one mean the counts vary more from
    week to week than a Poisson likelihood is willing to believe, and the size of
    the excess is what decides how hard that likelihood pushes a model to explain
    ordinary noise.
    """
    y = np.asarray(counts, dtype=float)
    smooth = (pd.Series(y).rolling(SMOOTH_WEEKS, center=True, min_periods=1)
              .mean().to_numpy())
    keep = smooth > 0
    if keep.sum() < MIN_USABLE:
        return float("nan")
    return float(np.mean((y - smooth)[keep] ** 2 / smooth[keep]))


def main() -> None:
    cases = pd.read_csv(
        raw / cfg["data"]["opendengue"]["csv_name"],
        usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                 "calendar_start_date", "calendar_end_date", "dengue_total"],
        low_memory=False)
    cases["start"] = pd.to_datetime(cases["calendar_start_date"], errors="coerce")
    cases["end"] = pd.to_datetime(cases["calendar_end_date"], errors="coerce")
    cases = cases[(cases["end"] - cases["start"]).dt.days + 1 == 7]

    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])

    rows = []
    for _, w in inv.iterrows():
        sub = cases[cases["adm_0_name"].astype(str).str.upper()
                    == str(w["country"]).upper()]
        if w["level"] == "national":
            sub = sub[sub["adm_1_name"].isna()]
        else:
            sub = sub[sub["adm_1_name"].astype(str).str.upper()
                      == str(w["unit"]).upper()]
            sub = sub[sub["adm_2_name"].isna()]
        sub = sub[(sub["start"] >= w["start"]) & (sub["start"] <= w["end"])]
        s = sub.groupby("start")["dengue_total"].sum().sort_index()
        if len(s) < MIN_WEEKS:
            continue
        idx = dispersion_index(s.to_numpy(float))
        if np.isfinite(idx):
            rows.append(dict(country=w["country"], unit=w["unit"],
                             window_start=str(w["start"].date()),
                             weeks=len(s), mean_weekly=float(s.mean()),
                             dispersion_index=round(idx, 3)))

    d = pd.DataFrame(rows)
    d.to_csv(tables / "41_dispersion.csv", index=False)
    r = d["dispersion_index"].to_numpy()

    print("=" * 74)
    print("HOW OVERDISPERSED ARE WEEKLY DENGUE COUNTS?")
    print("=" * 74)
    print(f"  {len(r)} outbreaks; residuals about a {SMOOTH_WEEKS}-week centred "
          f"moving average.")
    print("  A Poisson process gives an index of 1.\n")
    print(f"  median    {np.median(r):8.1f}")
    print(f"  quartiles {np.percentile(r, 25):8.1f} to {np.percentile(r, 75):.1f}")
    print(f"  90th pct  {np.percentile(r, 90):8.1f}")
    print(f"  maximum   {r.max():8.0f}")
    print(f"\n  above  5: {(r > 5).mean() * 100:5.0f}% of outbreaks")
    print(f"  above 20: {(r > 20).mean() * 100:5.0f}% of outbreaks")
    print(f"\n  A week of 2,000 cases is treated by Poisson as accurate to "
          f"+/- {np.sqrt(2000):.0f}.")
    print(f"  At the median dispersion the honest figure is "
          f"+/- {np.sqrt(2000 * np.median(r)):.0f}.")
    print("\n  'Two orders of magnitude', which an earlier draft claimed, is true")
    print("  only of the extreme tail. Six times at the median is the number.")
    print(f"\nTable: {tables / '41_dispersion.csv'}")


if __name__ == "__main__":
    main()
