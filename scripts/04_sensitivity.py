"""
Pipeline step 4 — sensitivity of the conclusions to choices that were made
rather than estimated.

Three quantities were fixed by judgement, and a reader is entitled to know how
much each one is holding up the result:

* **The reporting fraction.** Structurally unidentifiable from the population at
  risk, so it had to be fixed. The prediction is invariant to trading one
  against the other, which makes a strong prediction: the estimated population
  should scale exactly inversely with rho while every transmission parameter
  stays put. If it does not, something is wrong with the reformulation.
* **The rainfall lag.** Set to five weeks from the biological chain and from
  cross-correlation, but the empirical optimum differed by window (6, 5, 4).
* **The climate location.** A single city stands in for a province or a country.

Each is varied while everything else is held fixed. The quantity tracked is R0,
because that is what the conclusions rest on.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.climate import lagged_smoothed_rain, standardise  # noqa: E402
from dengue_pk.inference import Dataset, fit  # noqa: E402
from dengue_pk.models import FixedParams, basic_reproduction_number  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OBS = "nb"          # the observation model adopted after step 03

cfg = load_config()
fixed = FixedParams.from_config(cfg)
processed = resolve(cfg, "processed")
raw = resolve(cfg, "raw")
figures = resolve(cfg, "figures")
tables = resolve(cfg, "tables")

rows = []
TABLE = tables / "03_sensitivity.csv"

# Baseline fits, computed once per window with the full multi-start search. Every
# perturbed configuration warm-starts from the corresponding baseline and adds
# three random restarts.
#
# This is an efficiency choice with a cost worth stating. A sensitivity analysis
# asks how an estimate moves under a small change of assumption, so beginning
# near the unperturbed optimum is appropriate — but if a change of assumption
# introduced a new and distant optimum, three random restarts might miss it. The
# full search takes about four hours for this sweep against forty minutes;
# spot-checking a handful of settings at full strength is the guard against that,
# and any disagreement would invalidate the shortcut.
BASELINE: dict[tuple[str, str], dict] = {}


def fit_and_record(data, cfg_local, axis, value, window):
    for model in ("climate", "constant"):
        key = (window, model)
        try:
            if key not in BASELINE:
                res = fit(data, cfg_local, fixed, model=model, observation=OBS)
                BASELINE[key] = res.theta
            else:
                res = fit(data, cfg_local, fixed, model=model, observation=OBS,
                          start_from=BASELINE[key], n_starts_override=4)
        except Exception as exc:                      # a setting that will not fit
            rows.append(dict(window=window, axis=axis, value=value, model=model,
                             R0=np.nan, pop_at_risk_M=np.nan, aic=np.nan,
                             note=str(exc)[:60]))
            pd.DataFrame(rows).to_csv(TABLE, index=False)
            continue
        rows.append(dict(
            window=window, axis=axis, value=value, model=model,
            R0=round(basic_reproduction_number(res.theta["beta_0"], fixed), 3),
            pop_at_risk_M=round(res.theta["pop_frac"] * data.population / 1e6, 4),
            aic=round(res.aic, 1),
            a_temp=round(res.theta.get("a_temp", np.nan), 3),
            a_rain=round(res.theta.get("a_rain", np.nan), 3),
            nb_k=round(res.nb_k, 2) if res.nb_k else np.nan,
            starts=f"{res.n_converged}/{res.n_starts}", note=""))
        # Written after every fit: the sweep takes tens of minutes and a failure
        # in a later setting should not discard the earlier ones.
        pd.DataFrame(rows).to_csv(TABLE, index=False)


# ---------------------------------------------------------------------------
# 1. Reporting fraction
# ---------------------------------------------------------------------------
print("=== Sensitivity to the assumed reporting fraction ===")
for name, w in cfg["windows"].items():
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])
    data = Dataset.from_frame(df, label=w["label"])
    for rho in cfg["model"]["fixed"]["rho_sensitivity"]:
        cfg_local = copy.deepcopy(cfg)
        cfg_local["model"]["fixed"]["rho_fixed"] = rho
        print(f"  {name}, rho = {rho}")
        fit_and_record(data, cfg_local, "rho", rho, name)


# ---------------------------------------------------------------------------
# 2. Rainfall lag — rebuild the covariate rather than reusing the stored one
# ---------------------------------------------------------------------------
def load_climate(path: Path) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("YEAR,"))
    clim = pd.read_csv(path, skiprows=start).replace(-999.0, np.nan)
    clim["date"] = (pd.to_datetime(clim["YEAR"].astype(str), format="%Y")
                    + pd.to_timedelta(clim["DOY"] - 1, unit="D"))
    return clim.set_index("date").sort_index()


print("\n=== Sensitivity to the rainfall lag ===")
for name, w in cfg["windows"].items():
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])
    pt = w["climate_point"]
    clim = load_climate(raw / f"climate_{name}_{pt['name'].lower()}.csv")
    smooth = cfg["model"]["climate_forcing"]["rain_smooth_weeks"]

    for lag in cfg["model"]["climate_forcing"]["rain_lag_sensitivity"]:
        rain = lagged_smoothed_rain(clim["PRECTOTCORR"], lag, smooth)
        weekly = rain.resample("7D", origin=df["week_start"].min()).mean()
        aligned = weekly.reindex(df["week_start"]).to_numpy()
        if np.isnan(aligned).any():
            print(f"  {name}, lag {lag}w: incomplete coverage; skipped")
            continue
        z_r, *_ = standardise(aligned)
        data = Dataset(df["days_from_start"].to_numpy(float),
                       df["cases"].to_numpy(float), df["T2M"].to_numpy(float),
                       z_r, float(df["population"].iloc[0]), w["label"])
        print(f"  {name}, lag {lag}w")
        fit_and_record(data, cfg, "rain_lag_weeks", lag, name)


# ---------------------------------------------------------------------------
# 3. Climate location — does the choice of city matter?
# ---------------------------------------------------------------------------
ALTERNATIVES = {
    "national_2013": [("Lahore", 31.55, 74.35), ("Karachi", 24.86, 67.01),
                      ("Peshawar", 34.01, 71.58), ("Multan", 30.20, 71.45)],
    "sindh_2021": [("Karachi", 24.86, 67.01), ("Hyderabad", 25.40, 68.37),
                   ("Sukkur", 27.71, 68.83)],
    "kp_2021": [("Peshawar", 34.01, 71.58), ("Abbottabad", 34.15, 73.21),
                ("Bannu", 32.99, 70.60)],
}

print("\n=== Sensitivity to the climate location ===")
print("  (requires downloading additional NASA POWER series; skipped if absent)")
for name, w in cfg["windows"].items():
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])
    for city, lat, lon in ALTERNATIVES.get(name, []):
        path = raw / f"climate_alt_{name}_{city.lower()}.csv"
        if not path.exists():
            continue
        clim = load_climate(path)
        smooth = cfg["model"]["climate_forcing"]["rain_smooth_weeks"]
        lag = cfg["model"]["climate_forcing"]["rain_lag_weeks"]
        clim["rain_lagged"] = lagged_smoothed_rain(clim["PRECTOTCORR"], lag, smooth)
        weekly = clim.resample("7D", origin=df["week_start"].min()).mean(
            numeric_only=True).reindex(df["week_start"])
        if weekly[["T2M", "rain_lagged"]].isna().any().any():
            continue
        z_r, *_ = standardise(weekly["rain_lagged"].to_numpy())
        data = Dataset(df["days_from_start"].to_numpy(float),
                       df["cases"].to_numpy(float),
                       weekly["T2M"].to_numpy(float), z_r,
                       float(df["population"].iloc[0]), w["label"])
        print(f"  {name}, climate from {city}")
        fit_and_record(data, cfg, "climate_location", city, name)


# ---------------------------------------------------------------------------
tab = pd.DataFrame(rows)
tab.to_csv(TABLE, index=False)
print(f"\nTable: {TABLE}  ({len(tab)} fits)")

print("\n--- R0 across every setting, climate model ---")
clim_only = tab[(tab["model"] == "climate") & tab["R0"].notna()]
for (window, axis), grp in clim_only.groupby(["window", "axis"]):
    r0 = grp["R0"]
    print(f"  {window:15s} {axis:17s} R0 {r0.min():.2f}–{r0.max():.2f} "
          f"(spread {r0.max() - r0.min():.2f} over {len(grp)} settings)")

# The rho check is a prediction, not an observation: the estimated population
# must scale exactly inversely with the assumed reporting fraction.
print("\n--- Does the population at risk scale as 1/rho, as the algebra requires? ---")
for window, grp in clim_only[clim_only["axis"] == "rho"].groupby("window"):
    g = grp.sort_values("value")
    product = g["pop_at_risk_M"] * g["value"]
    print(f"  {window:15s} rho x population = "
          f"{', '.join(f'{p:.4f}' for p in product)}"
          f"   (constant to {product.std() / product.mean() * 100:.1f}%)")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
for ax, axis in zip(axes, ["rho", "rain_lag_weeks", "climate_location"]):
    sub = clim_only[clim_only["axis"] == axis]
    if sub.empty:
        ax.set_axis_off()
        ax.set_title(f"{axis}: not run")
        continue
    for window, grp in sub.groupby("window"):
        ax.plot(grp["value"].astype(str), grp["R0"], "o-", label=window)
    ax.set_title(f"R₀ against {axis}")
    ax.set_ylabel("R₀")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(figures / "03_sensitivity.png", dpi=150)
print(f"\nFigure: {figures / '03_sensitivity.png'}")
