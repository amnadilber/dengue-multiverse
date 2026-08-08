"""
Pipeline step 15 — the study this project exists to make possible.

On the Pakistani windows, three routine analysis choices each reversed the
in-sample verdict on whether transmission is climate-driven. Three windows make
that an anecdote. This runs the same factorial across every usable outbreak in
the global inventory and asks how often it happens.

Each window is fitted under every combination of four choices, none of which a
reader of a finished paper would normally see stated:

* **observation model** — Poisson or negative binomial
* **temperature term** — a unimodal response from vector biology, or the
  conventional free log-linear coefficient
* **rainfall lag** — 3, 5 or 7 weeks
* **how much of the series is fitted** — all of it, or the first 75%

Twenty-four combinations, two models each. For every window the outcome recorded
is which model AIC prefers, and separately how the two compare on weeks withheld
from fitting.

The quantity of interest is not how often climate forcing wins. It is how often
the answer **changes** within a single window as these choices vary, because a
conclusion that moves with the analyst's unremarked decisions is not evidence
about the world.

Run with a positional argument to limit the number of windows, for timing.
"""

from __future__ import annotations

import copy
import itertools
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.climate import lagged_smoothed_rain, standardise  # noqa: E402
from dengue_pk.inference import (Dataset, fit, nb_deviance_residuals,  # noqa: E402
                                 poisson_deviance_residuals, predict,
                                 set_temperature_form)
from dengue_pk.locations import point_for  # noqa: E402
from dengue_pk.models import (FixedParams, IntegrationFailure,  # noqa: E402
                              basic_reproduction_number)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OBSERVATIONS = ("nb", "poisson")
TEMP_FORMS = ("briere", "loglinear")
RAIN_LAGS = (3, 5, 7)
TRAIN_FRACS = (1.0, 0.75)
# Model structure is the fifth factor and the one a referee is most likely to
# raise: a result obtained under a single transmission model might be a property
# of that model. `seir` is the directly transmitted formulation much of the
# applied literature actually fits to case counts.
STRUCTURES = ("hostvector", "seir")

BASE = load_config()
fixed = FixedParams.from_config(BASE)
raw = resolve(BASE, "raw")
tables = resolve(BASE, "tables")
clim_dir = raw / "climate_global"

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else None
OUT = tables / ("20_global_robustness_5factor.csv" if LIMIT is None
                else "20_global_robustness_5factor_pilot.csv")

inv = pd.read_csv(tables / "12_global_windows.csv", parse_dates=["start", "end"])
if LIMIT:
    # Spread the pilot across the inventory rather than taking the first N,
    # which would be one country.
    inv = inv.iloc[:: max(len(inv) // LIMIT, 1)].head(LIMIT)

print(f"{len(inv)} windows x {len(OBSERVATIONS) * len(TEMP_FORMS) * len(RAIN_LAGS) * len(TRAIN_FRACS)} "
      f"combinations x 2 models\n")

# --- case series, loaded once ----------------------------------------------
print("Loading case data...")
cases_all = pd.read_csv(raw / BASE["data"]["opendengue"]["csv_name"],
                        usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                                 "calendar_start_date", "calendar_end_date",
                                 "dengue_total"], low_memory=False)
cases_all["start"] = pd.to_datetime(cases_all["calendar_start_date"], errors="coerce")
cases_all["end"] = pd.to_datetime(cases_all["calendar_end_date"], errors="coerce")
cases_all = cases_all[(cases_all["end"] - cases_all["start"]).dt.days + 1 == 7]

_climate_cache: dict[str, pd.DataFrame] = {}


def load_climate(slug: str):
    if slug not in _climate_cache:
        path = clim_dir / f"{slug}.csv"
        if not path.exists():
            return None
        with open(path, "r", encoding="utf-8") as fh:
            lines = fh.readlines()
        try:
            head = next(i for i, ln in enumerate(lines) if ln.startswith("YEAR,"))
        except StopIteration:
            return None
        c = pd.read_csv(path, skiprows=head).replace(-999.0, np.nan)
        c["date"] = (pd.to_datetime(c["YEAR"].astype(str), format="%Y")
                     + pd.to_timedelta(c["DOY"] - 1, unit="D"))
        _climate_cache[slug] = c.set_index("date").sort_index()
    return _climate_cache[slug]


def case_series(w) -> pd.Series:
    sub = cases_all[cases_all["adm_0_name"].astype(str).str.upper()
                    == str(w["country"]).upper()]
    if w["level"] == "national":
        sub = sub[sub["adm_1_name"].isna()]
    else:
        sub = sub[sub["adm_1_name"].astype(str).str.upper() == str(w["unit"]).upper()]
        sub = sub[sub["adm_2_name"].isna()]
    sub = sub[(sub["start"] >= w["start"]) & (sub["start"] <= w["end"])]
    return sub.groupby("start")["dengue_total"].sum().sort_index()


def build(w, lag_weeks):
    """Dataset for one window at one rainfall lag, or None if unbuildable."""
    pt = point_for(w["country"], w["unit"], w["level"])
    if pt is None:
        return None, "no climate point"
    name, lat, lon = pt
    slug = f"{name.lower().replace(' ', '_').replace('/', '_')}_{round(lat, 4)}_{round(lon, 4)}"
    clim = load_climate(slug)
    if clim is None:
        return None, "climate file missing"

    s = case_series(w)
    if len(s) < 20:
        return None, f"only {len(s)} weekly records"
    gaps = s.index.to_series().diff().dt.days.dropna()
    if not (gaps == 7).all():
        return None, "weekly grid broken"

    smooth = BASE["model"]["climate_forcing"]["rain_smooth_weeks"]
    c = clim.copy()
    c["rain_lagged"] = lagged_smoothed_rain(c["PRECTOTCORR"], lag_weeks, smooth)
    weekly = (c.resample("7D", origin=s.index.min())
               .mean(numeric_only=True).reindex(s.index))
    if weekly[["T2M", "rain_lagged"]].isna().any().any():
        return None, "climate coverage incomplete"

    z_rain, *_ = standardise(weekly["rain_lagged"].to_numpy())
    days = (s.index - s.index.min()).days.to_numpy(float)
    # Population is a scale only: R0 is invariant to it, and pop_frac absorbs it.
    return Dataset(days, s.to_numpy(float), weekly["T2M"].to_numpy(float),
                   z_rain, 1e7, f"{w['country']}/{w['unit']}"), None


def cfg_for(observation, temp_form, lag, structure="hostvector"):
    c = copy.deepcopy(BASE)
    c["model"]["temperature_form"] = temp_form
    c["model"]["structure"] = structure
    c["model"]["climate_forcing"]["rain_lag_weeks"] = lag
    c["inference"]["observation"] = observation
    c["inference"]["classical"]["n_restarts"] = 3
    return c


rows = []
t_start = time.time()
skipped = []

for idx, (_, w) in enumerate(inv.iterrows(), 1):
    label = f"{w['country']}/{w['unit']} {w['start'].date()}"
    data_by_lag = {}
    for lag in RAIN_LAGS:
        d, why = build(w, lag)
        if d is not None:
            data_by_lag[lag] = d
    if not data_by_lag:
        skipped.append((label, why))
        continue

    # One careful fit per window sets the dispersion and the warm start. Every
    # combination then reuses them, which is what makes 24 x 2 fits per window
    # affordable; the alternative is a multi-start search per combination and a
    # runtime measured in days.
    base_data = data_by_lag[BASE["model"]["climate_forcing"]["rain_lag_weeks"]
                            if BASE["model"]["climate_forcing"]["rain_lag_weeks"]
                            in data_by_lag else min(data_by_lag)]
    set_temperature_form("briere")
    try:
        anchor = fit(base_data, cfg_for("nb", "briere", 5), fixed,
                     model="constant", observation="nb")
    except Exception as exc:
        skipped.append((label, f"anchor fit failed: {str(exc)[:40]}"))
        continue
    k_hat = anchor.nb_k or 1e6
    warm = dict(anchor.theta)
    warm.setdefault("a_temp", 0.5)
    warm.setdefault("a_rain", 0.0)

    n_ok = 0
    for observation, temp_form, lag, frac, structure in itertools.product(
            OBSERVATIONS, TEMP_FORMS, RAIN_LAGS, TRAIN_FRACS, STRUCTURES):
        if lag not in data_by_lag:
            continue
        full = data_by_lag[lag]
        k_train = int(round(len(full.days) * frac))
        train = full if frac >= 1.0 else full.head(k_train)
        cfg = cfg_for(observation, temp_form, lag, structure)
        set_temperature_form(temp_form)

        res = {}
        for model in ("climate", "constant"):
            try:
                res[model] = fit(train, cfg, fixed, model=model,
                                 observation=observation,
                                 start_from=warm, n_starts_override=3,
                                 fixed_nb_k=k_hat if observation == "nb" else None)
            except (RuntimeError, IntegrationFailure, ValueError):
                res[model] = None
        if res["climate"] is None or res["constant"] is None:
            continue

        d_aic = res["climate"].aic - res["constant"].aic
        held = np.nan
        if frac < 1.0:
            try:
                tail = slice(k_train, None)
                dev = {}
                for model in ("climate", "constant"):
                    mu = predict(res[model].theta, full, fixed, model, cfg)
                    dev[model] = float(np.sum(
                        (nb_deviance_residuals(full.cases[tail], mu[tail], k_hat)
                         if observation == "nb"
                         else poisson_deviance_residuals(full.cases[tail], mu[tail]))
                        ** 2))
                held = dev["climate"] - dev["constant"]
            except (IntegrationFailure, ValueError):
                held = np.nan

        rows.append(dict(
            country=w["country"], unit=w["unit"], level=w["level"],
            window_start=w["start"].date(), weeks=int(w["weeks"]),
            cases=int(w["cases"]),
            observation=observation, temp_form=temp_form, rain_lag=lag,
            train_frac=frac, structure=structure,
            delta_aic=round(d_aic, 3),
            climate_wins=bool(d_aic < 0),
            heldout_delta=round(held, 3) if np.isfinite(held) else np.nan,
            R0_climate=round(basic_reproduction_number(
                res["climate"].theta["beta_0"], fixed, structure), 4),
            R0_constant=round(basic_reproduction_number(
                res["constant"].theta["beta_0"], fixed, structure), 4),
            a_temp=round(res["climate"].theta.get("a_temp", np.nan), 4),
            a_rain=round(res["climate"].theta.get("a_rain", np.nan), 4)))
        n_ok += 1

    pd.DataFrame(rows).to_csv(OUT, index=False)
    rate = (time.time() - t_start) / idx
    print(f"[{idx:3d}/{len(inv)}] {label:45s} {n_ok:2d} combos  "
          f"{rate:5.1f} s/window  eta {rate * (len(inv) - idx) / 60:5.1f} min")

print(f"\n{len(rows)} fits recorded in {(time.time() - t_start) / 60:.1f} min")
if skipped:
    print(f"{len(skipped)} windows skipped:")
    for label, why in skipped[:15]:
        print(f"  {label:45s} {why}")
print(f"Table: {OUT}")
