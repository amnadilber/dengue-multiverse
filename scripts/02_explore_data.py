"""
Pipeline step 2 — inspect the assembled series before modelling anything.

Looking at the data first is not a formality. The figure produced here decides
whether the modelling assumptions are defensible: whether each window contains a
complete wave, whether the climate covariates vary enough to be identifiable
from a single season, and whether the seasonal relationship the model presumes
is visible at all.

Anything the figure contradicts must change the model, not be explained away.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
processed = resolve(cfg, "processed")
figures = resolve(cfg, "figures")
tables = resolve(cfg, "tables")

windows = list(cfg["windows"].items())

fig, axes = plt.subplots(len(windows), 1, figsize=(11, 3.6 * len(windows)))
summary = []

for ax, (name, w) in zip(np.atleast_1d(axes), windows):
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])

    ax.bar(df["week_start"], df["cases"], width=6, color="#8c1c13",
           label="reported cases")
    ax.set_ylabel("cases / week")
    ax.set_title(f"{w['label']} — {len(df)} weeks")
    ax.grid(alpha=0.25, axis="y")

    ax_t = ax.twinx()
    ax_t.plot(df["week_start"], df["T2M"], color="#e08a1e", lw=1.6,
              label="temperature")
    ax_t.plot(df["week_start"], df["rain_lagged"] * 5, color="#1f6f8b", lw=1.6,
              ls="--", label=f"lagged rain x5 ({cfg['model']['climate_forcing']['rain_lag_weeks']}w)")
    ax_t.set_ylabel("°C   /   mm day⁻¹ (×5)")

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax_t.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")

    # How far behind the climate driver does the epidemic follow?
    #
    # Comparing the single largest case week against the single largest climate
    # week is unreliable: the 77-week national window spans two summers, so the
    # hottest week in it need not belong to the epidemic season at all. Lagged
    # cross-correlation uses the whole series instead and returns the delay that
    # actually maximises association, which is the quantity the model needs.
    best = {}
    for driver in ("T2M", "PRECTOTCORR"):
        cors = {L: df["cases"].corr(df[driver].shift(L))
                for L in range(0, min(13, len(df) // 3))}
        lag = max(cors, key=lambda k: (cors[k] if pd.notna(cors[k]) else -np.inf))
        best[driver] = (lag, round(float(cors[lag]), 3))

    summary.append(dict(
        window=name, weeks=len(df), total_cases=int(df["cases"].sum()),
        peak_cases=int(df["cases"].max()),
        peak_week=df.loc[df["cases"].idxmax(), "week_start"].date(),
        best_temp_lag_wks=best["T2M"][0], temp_corr=best["T2M"][1],
        best_rain_lag_wks=best["PRECTOTCORR"][0], rain_corr=best["PRECTOTCORR"][1],
        temp_range_C=round(df["T2M"].max() - df["T2M"].min(), 1),
        attack_rate_per_100k=round(df["cases"].sum() / df["population"].iloc[0] * 1e5, 2),
    ))

fig.tight_layout()
fig.savefig(figures / "01_data_overview.png", dpi=150)
print(f"Figure: {figures / '01_data_overview.png'}")

tab = pd.DataFrame(summary)
tab.to_csv(tables / "01_data_summary.csv", index=False)
print(f"Table:  {tables / '01_data_summary.csv'}\n")
print(tab.to_string(index=False))

configured = cfg["model"]["climate_forcing"]["rain_lag_weeks"]
print("\nReading:")
for r in summary:
    print(f"  {r['window']}: cases correlate best with rainfall lagged "
          f"{r['best_rain_lag_wks']} weeks (r = {r['rain_corr']}) and with "
          f"temperature lagged {r['best_temp_lag_wks']} weeks (r = {r['temp_corr']}); "
          f"temperature spans {r['temp_range_C']} C")

rain_lags = [r["best_rain_lag_wks"] for r in summary]
print(f"\nConfigured rainfall lag: {configured} weeks. "
      f"Empirical best lags: {rain_lags}.")
if configured not in rain_lags:
    print("  The configured lag is not the empirical optimum in any window. "
          "Either justify it from the vector-biology literature or refit with the "
          "lag estimated rather than assumed — do not leave the discrepancy "
          "unexplained.")
