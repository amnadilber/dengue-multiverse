"""
Pipeline step 1 — build analysis-ready weekly series.

For each study window this script extracts the case counts, aligns the daily
climate record to the same weeks, applies the mosquito development lag to
rainfall, and writes one tidy CSV per window.

Two checks are enforced rather than assumed:

* the weekly grid must be unbroken, since a gap would silently shift every
  subsequent week's alignment with the climate covariates;
* climate coverage must extend far enough before the first modelled week to
  supply the lagged rainfall term without back-filling into the fitting period.

A window failing either check is reported and skipped rather than quietly
producing a series that looks usable.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
# Enforced by tests/test_environment.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.climate import lagged_smoothed_rain  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def load_cases(csv_path: Path, window: dict) -> pd.DataFrame:
    """Weekly reported cases for one window, at the configured admin level."""
    usecols = ["adm_0_name", "adm_1_name", "adm_2_name", "calendar_start_date",
               "calendar_end_date", "dengue_total", "S_res", "T_res"]
    df = pd.read_csv(csv_path, usecols=usecols, low_memory=False)
    df = df[df["adm_0_name"].astype(str).str.upper() == "PAKISTAN"].copy()

    df["start"] = pd.to_datetime(df["calendar_start_date"], errors="coerce")
    df["end"] = pd.to_datetime(df["calendar_end_date"], errors="coerce")
    df = df[(df["end"] - df["start"]).dt.days + 1 == 7]          # weekly only

    if window["admin_level"] == 0:
        sub = df[df["adm_1_name"].isna()]
    else:
        sub = df[df["adm_1_name"].astype(str).str.upper() == window["region"]]

    sub = sub[(sub["start"] >= window["start"]) & (sub["start"] <= window["end"])]
    weekly = (sub.groupby("start")["dengue_total"].sum().sort_index()
                 .rename("cases").to_frame())
    return weekly


def load_climate(path: Path) -> pd.DataFrame:
    """Read a NASA POWER CSV, skipping its header block and marking gaps."""
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("YEAR,"))
    clim = pd.read_csv(path, skiprows=start)
    clim = clim.replace(-999.0, np.nan)
    clim["date"] = (pd.to_datetime(clim["YEAR"].astype(str), format="%Y")
                    + pd.to_timedelta(clim["DOY"] - 1, unit="D"))
    return clim.set_index("date").sort_index()


def build(window_name: str, window: dict, cases_csv: Path, raw_dir: Path,
          out_dir: Path) -> bool:
    print(f"\n--- {window['label']} ---")
    weekly = load_cases(cases_csv, window)
    if weekly.empty:
        print("  no weekly case records in this window; skipped")
        return False

    # Continuity check.
    gaps = weekly.index.to_series().diff().dt.days.dropna()
    if not (gaps == 7).all():
        bad = int((gaps != 7).sum())
        print(f"  {len(weekly)} weeks, but {bad} break(s) in the weekly grid; skipped")
        return False
    print(f"  {len(weekly)} unbroken weeks, "
          f"{weekly.index.min().date()} to {weekly.index.max().date()}")
    print(f"  {weekly['cases'].sum():,.0f} cases, peak {weekly['cases'].max():,.0f}/week")

    pt = window["climate_point"]
    clim_path = raw_dir / f"climate_{window_name}_{pt['name'].lower()}.csv"
    clim = load_climate(clim_path)

    cfg_forcing = load_config()["model"]["climate_forcing"]
    lag, smooth = cfg_forcing["rain_lag_weeks"], cfg_forcing["rain_smooth_weeks"]

    needed_before = (lag + smooth) * 7
    available_before = (weekly.index.min() - clim.index.min()).days
    if available_before < needed_before:
        print(f"  climate record starts only {available_before} days before the "
              f"first week; {needed_before} needed for the rainfall lag; skipped")
        return False

    clim["rain_lagged"] = lagged_smoothed_rain(clim["PRECTOTCORR"], lag, smooth)

    # Weekly means of the daily climate, labelled by week start.
    weekly_clim = (clim.resample("7D", origin=weekly.index.min())
                       .mean(numeric_only=True))
    joined = weekly.join(weekly_clim[["T2M", "T2M_MAX", "T2M_MIN", "RH2M",
                                      "PRECTOTCORR", "rain_lagged"]], how="left")

    if joined[["T2M", "rain_lagged"]].isna().any().any():
        print("  climate covariates incomplete after alignment; skipped")
        return False

    joined.insert(0, "week_index", np.arange(len(joined)))
    joined.insert(1, "days_from_start",
                  (joined.index - joined.index.min()).days.astype(float))
    joined["population"] = window["population"]

    out = out_dir / f"{window_name}.csv"
    joined.to_csv(out, index_label="week_start")
    print(f"  written: {out.name}  ({len(joined)} rows x {joined.shape[1]} cols)")
    print(f"  mean weekly temperature {joined['T2M'].mean():.1f} C, "
          f"lagged rainfall {joined['rain_lagged'].mean():.2f} mm/day")
    return True


def main() -> None:
    cfg = load_config()
    raw = resolve(cfg, "raw")
    out = resolve(cfg, "processed")
    cases_csv = raw / cfg["data"]["opendengue"]["csv_name"]

    built = [name for name, w in cfg["windows"].items()
             if build(name, w, cases_csv, raw, out)]

    print(f"\nBuilt {len(built)} of {len(cfg['windows'])} windows: {', '.join(built)}")
    if not built:
        raise SystemExit("No usable window; the analysis cannot proceed.")


if __name__ == "__main__":
    main()
