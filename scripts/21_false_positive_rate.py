"""
Pipeline step 21 — what is the false-positive rate of this analysis?

Everything so far measures disagreement: how often the verdict changes when the
analyst's choices change. Disagreement is evidence that something is wrong, but
it does not say which verdict is wrong, because on real data there is no truth to
compare against.

Here there is. Counts are simulated from a model whose answer is known, the full
factorial is run on them, and the verdict is scored against the truth:

* **Null truth.** Data generated from the constant-transmission model, with the
  real climate covariates present but exerting no influence. Any endorsement of
  climate forcing is a false positive, and its rate under each choice is the
  quantity that matters most.
* **Climate truth.** Data generated from the climate-forced model with a genuine
  effect. Failures to detect it are false negatives, and a rule that avoids false
  positives by never detecting anything must be caught here.

Simulation uses each window's own fitted parameters, its own climate series and
its own estimated dispersion, so the synthetic data resemble the real thing in
scale, shape and noise rather than being idealised.

This turns "the verdict is unstable" into "the verdict is wrong this often, in
this direction, under this choice", which is the form in which a methodological
finding is actionable.
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
from dengue_pk.inference import (Dataset, fit, make_forcing,  # noqa: E402
                                 predict, set_temperature_form)
from dengue_pk.locations import point_for  # noqa: E402
from dengue_pk.models import FixedParams, IntegrationFailure  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OBSERVATIONS = ("nb", "poisson")
TEMP_FORMS = ("briere", "loglinear")
RAIN_LAGS = (3, 5, 7)
TRAIN_FRACS = (1.0, 0.75)
STRUCTURES = ("hostvector", "seir")

# A genuine but moderate climate effect for the alternative truth: the rainfall
# coefficient is set to the median magnitude actually estimated on real data, so
# the power calculation is against an effect of realistic size rather than one
# chosen to be easy to detect.
TRUE_A_RAIN = 0.30
TRUE_A_TEMP = 1.00

BASE = load_config()
fixed = FixedParams.from_config(BASE)
raw = resolve(BASE, "raw")
tables = resolve(BASE, "tables")
clim_dir = raw / "climate_global"

N_WINDOWS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
OUT = tables / "21_false_positive_rate.csv"

inv = pd.read_csv(tables / "12_global_windows.csv", parse_dates=["start", "end"])
inv = inv.iloc[:: max(len(inv) // N_WINDOWS, 1)].head(N_WINDOWS)
rng = np.random.default_rng(BASE["seed"])

n_combos = (len(OBSERVATIONS) * len(TEMP_FORMS) * len(RAIN_LAGS)
            * len(TRAIN_FRACS) * len(STRUCTURES))
print(f"{len(inv)} windows x 2 truths x {n_combos} combinations x 2 models\n")

print("Loading case data (for the real series' scale and timing)...")
cases_all = pd.read_csv(raw / BASE["data"]["opendengue"]["csv_name"],
                        usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                                 "calendar_start_date", "calendar_end_date",
                                 "dengue_total"], low_memory=False)
cases_all["start"] = pd.to_datetime(cases_all["calendar_start_date"], errors="coerce")
cases_all["end"] = pd.to_datetime(cases_all["calendar_end_date"], errors="coerce")
cases_all = cases_all[(cases_all["end"] - cases_all["start"]).dt.days + 1 == 7]

_clim: dict[str, pd.DataFrame] = {}


def climate_for(w):
    pt = point_for(w["country"], w["unit"], w["level"])
    if pt is None:
        return None
    name, lat, lon = pt
    slug = f"{name.lower().replace(' ', '_').replace('/', '_')}_{round(lat, 4)}_{round(lon, 4)}"
    if slug not in _clim:
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
        _clim[slug] = c.set_index("date").sort_index()
    return _clim[slug]


def build(w, lag):
    clim = climate_for(w)
    if clim is None:
        return None
    sub = cases_all[cases_all["adm_0_name"].astype(str).str.upper()
                    == str(w["country"]).upper()]
    if w["level"] == "national":
        sub = sub[sub["adm_1_name"].isna()]
    else:
        sub = sub[sub["adm_1_name"].astype(str).str.upper() == str(w["unit"]).upper()]
        sub = sub[sub["adm_2_name"].isna()]
    sub = sub[(sub["start"] >= w["start"]) & (sub["start"] <= w["end"])]
    s = sub.groupby("start")["dengue_total"].sum().sort_index()
    if len(s) < 20:
        return None
    gaps = s.index.to_series().diff().dt.days.dropna()
    if not (gaps == 7).all():
        return None
    c = clim.copy()
    c["rain_lagged"] = lagged_smoothed_rain(
        c["PRECTOTCORR"], lag, BASE["model"]["climate_forcing"]["rain_smooth_weeks"])
    weekly = (c.resample("7D", origin=s.index.min())
               .mean(numeric_only=True).reindex(s.index))
    if weekly[["T2M", "rain_lagged"]].isna().any().any():
        return None
    z_rain, *_ = standardise(weekly["rain_lagged"].to_numpy())
    days = (s.index - s.index.min()).days.to_numpy(float)
    return Dataset(days, s.to_numpy(float), weekly["T2M"].to_numpy(float),
                   z_rain, 1e7, f"{w['country']}/{w['unit']}")


def cfg_for(observation, temp_form, lag, structure):
    c = copy.deepcopy(BASE)
    c["model"]["temperature_form"] = temp_form
    c["model"]["structure"] = structure
    c["model"]["climate_forcing"]["rain_lag_weeks"] = lag
    c["inference"]["observation"] = observation
    c["inference"]["classical"]["n_restarts"] = 3
    return c


rows = []
t0 = time.time()
done = 0

for idx, (_, w) in enumerate(inv.iterrows(), 1):
    label = f"{w['country']}/{w['unit']} {w['start'].date()}"
    data_by_lag = {lag: build(w, lag) for lag in RAIN_LAGS}
    data_by_lag = {k: v for k, v in data_by_lag.items() if v is not None}
    if not data_by_lag:
        continue
    ref = data_by_lag[min(data_by_lag)]

    # Fit the real series once to obtain realistic generating parameters.
    set_temperature_form("briere")
    base_cfg = cfg_for("nb", "briere", min(data_by_lag), "hostvector")
    try:
        anchor = fit(ref, base_cfg, fixed, model="constant", observation="nb")
    except Exception:
        continue
    k_hat = float(anchor.nb_k or 20.0)
    # Labelled "no_effect" rather than "null": pandas reads the string
    # "null" back as NaN, which silently deleted an entire arm of this
    # study once already (see docs/ANALYSIS_LOG.md).
    truth_null = dict(anchor.theta)
    truth_null.update(a_temp=1e-9, a_rain=0.0)
    truth_climate = dict(anchor.theta)
    truth_climate.update(a_temp=TRUE_A_TEMP, a_rain=TRUE_A_RAIN)
    # Match the two truths on mean transmission.
    #
    # Switching the climate coefficients on without this multiplies beta(t) by a
    # factor whose average is not 1, so the alternative differs from the null in
    # overall transmission intensity as well as in whether beta varies. That
    # confounds the thing being measured, and it also broke the study: for 21 of
    # 59 windows the raised R_0 made the generating integration diverge, so power
    # was estimated only on windows where a strong effect happened not to break
    # the model — a subsample selected on something related to the outcome.
    #
    # Rescaling beta_0 by the mean of the climate multiplier leaves both truths
    # with the same average transmission, differing only in whether it varies
    # with climate. That is the contrast this study is meant to test.
    mult = make_forcing({**truth_climate, "beta_0": 1.0}, ref,
                        "climate", base_cfg).on_grid(ref.days)
    truth_climate["beta_0"] = float(anchor.theta["beta_0"]
                                    / max(float(np.mean(mult)), 1e-12))

    for truth_name, truth in (("no_effect", truth_null), ("climate", truth_climate)):
        try:
            mu_true = predict(truth, ref, fixed, "climate", base_cfg)
        except IntegrationFailure:
            continue
        # Gamma-Poisson mixture reproduces NB2 with the window's own dispersion.
        lam = rng.gamma(shape=k_hat, scale=np.maximum(mu_true, 1e-12) / k_hat)
        synth_counts = rng.poisson(lam).astype(float)
        if synth_counts.sum() < 100:
            continue

        for observation, temp_form, lag, frac, structure in itertools.product(
                OBSERVATIONS, TEMP_FORMS, RAIN_LAGS, TRAIN_FRACS, STRUCTURES):
            if lag not in data_by_lag:
                continue
            d0 = data_by_lag[lag]
            full = Dataset(d0.days, synth_counts, d0.temp_c, d0.z_rain,
                           d0.population, d0.label)
            train = full if frac >= 1.0 else full.head(
                int(round(len(full.days) * frac)))
            cfg = cfg_for(observation, temp_form, lag, structure)
            set_temperature_form(temp_form)
            res = {}
            for model in ("climate", "constant"):
                try:
                    res[model] = fit(train, cfg, fixed, model=model,
                                     observation=observation,
                                     start_from=truth, n_starts_override=3,
                                     fixed_nb_k=k_hat if observation == "nb" else None)
                except (RuntimeError, IntegrationFailure, ValueError):
                    res[model] = None
            if res["climate"] is None or res["constant"] is None:
                continue
            d_aic = res["climate"].aic - res["constant"].aic
            rows.append(dict(
                country=w["country"], unit=w["unit"],
                window_start=w["start"].date(), truth=truth_name,
                observation=observation, temp_form=temp_form, rain_lag=lag,
                train_frac=frac, structure=structure,
                delta_aic=round(d_aic, 3),
                climate_wins=bool(d_aic < 0),
                climate_wins_margin4=bool(d_aic < -4)))
    done += 1
    pd.DataFrame(rows).to_csv(OUT, index=False)
    rate = (time.time() - t0) / idx
    print(f"[{idx:3d}/{len(inv)}] {label:42s} {rate:5.1f} s/window  "
          f"eta {rate * (len(inv) - idx) / 60:5.1f} min")

d = pd.DataFrame(rows)
print(f"\n{len(d):,} fits from {done} windows in {(time.time() - t0) / 60:.1f} min\n")

if len(d):
    print("=" * 78)
    print("FALSE POSITIVES: data generated with NO climate effect")
    print("=" * 78)
    null = d[d["truth"] == "no_effect"]
    print(f"  climate endorsed by the sign of dAIC : "
          f"{null['climate_wins'].mean() * 100:5.1f}% of fits")
    print(f"  climate endorsed with a margin of 4  : "
          f"{null['climate_wins_margin4'].mean() * 100:5.1f}% of fits\n")
    for col in ("observation", "temp_form", "structure", "train_frac"):
        g = null.groupby(col)[["climate_wins", "climate_wins_margin4"]].mean() * 100
        print(f"  by {col}:")
        for lvl, r in g.iterrows():
            print(f"    {str(lvl):12s} sign {r['climate_wins']:5.1f}%   "
                  f"margin 4 {r['climate_wins_margin4']:5.1f}%")

    print("\n" + "=" * 78)
    print("POWER: data generated WITH a real climate effect")
    print("=" * 78)
    clim = d[d["truth"] == "climate"]
    print(f"  climate detected by the sign of dAIC : "
          f"{clim['climate_wins'].mean() * 100:5.1f}% of fits")
    print(f"  climate detected with a margin of 4  : "
          f"{clim['climate_wins_margin4'].mean() * 100:5.1f}% of fits")
    for col in ("observation", "structure"):
        g = clim.groupby(col)[["climate_wins", "climate_wins_margin4"]].mean() * 100
        print(f"  by {col}:")
        for lvl, r in g.iterrows():
            print(f"    {str(lvl):12s} sign {r['climate_wins']:5.1f}%   "
                  f"margin 4 {r['climate_wins_margin4']:5.1f}%")

print(f"\nTable: {OUT}")
