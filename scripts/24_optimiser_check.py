"""Pipeline step 24 — is the instability optimiser noise?

The sharpest remaining objection to the robustness result is that it may not be
about analysis choices at all. The global run fitted each model with three
multi-starts from a warm start supplied by an anchor fit, which is a compromise
made for runtime: 11,326 fits at ten restarts would have taken most of a week.
If three starts sometimes miss the optimum, and miss it more often for the
climate model — which has two extra parameters and therefore a rougher surface —
then Delta AIC would move for reasons that have nothing to do with the
observation model or the rainfall lag, and what this study calls instability
would partly be the optimiser failing to converge.

The objection is testable, so it is tested rather than argued with. A sample of
windows is refitted twice under the same 48 combinations:

* **as-run** — three restarts, warm-started from the anchor fit, exactly as in
  step 16;
* **thorough** — ten restarts, cold, no warm start, which is the configured
  default and roughly three times the work.

Three things are then measured, and they answer different questions:

1. **Does the thorough arm find better optima at all?** If it does not, the
   objection is empty because there was nothing left to find. Measured as the
   share of fits where the thorough arm reaches a strictly lower AIC, and by how
   much.
2. **Does the verdict change?** The share of individual comparisons whose sign
   of Delta AIC differs between the arms. This is the direct analogue of the
   flip rates reported for the five analysis choices, and can be read on the
   same scale: if the optimiser flips fewer verdicts than the weakest genuine
   factor, it is not what is driving the result.
3. **Does the headline change?** The share of windows whose verdict is unstable,
   computed within each arm separately. This is the number the paper reports,
   and it is the one that has to survive.

A finding of "the thorough arm finds better optima but the instability is
unchanged" would be the strongest possible answer: it would mean the optimiser
does have room to improve and that improving it does not rescue the conclusion.
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
from dengue_pk.inference import Dataset, fit, set_temperature_form  # noqa: E402
from dengue_pk.locations import point_for  # noqa: E402
from dengue_pk.models import FixedParams, IntegrationFailure  # noqa: E402
from dengue_pk.robustness import complete_windows, window_verdicts  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OBSERVATIONS = ("nb", "poisson")
TEMP_FORMS = ("briere", "loglinear")
RAIN_LAGS = (3, 5, 7)
TRAIN_FRACS = (1.0, 0.75)
STRUCTURES = ("hostvector", "seir")

#: (label, restarts, warm-started) for the two arms.
ARMS = (("as_run", 3, True), ("thorough", 10, False))

BASE = load_config()
fixed = FixedParams.from_config(BASE)
raw = resolve(BASE, "raw")
tables = resolve(BASE, "tables")
clim_dir = raw / "climate_global"

OUT = tables / "23_optimiser_check.csv"

# `analyse` re-reads the stored table and reports without refitting. The fitting
# arm writes after every window, so a run stopped early — this one competes for
# the machine with step 21 and was stopped once it had enough windows to answer
# the question — still leaves a complete, analysable table.
ANALYSE_ONLY = len(sys.argv) > 1 and sys.argv[1] == "analyse"
N_WINDOWS = (0 if ANALYSE_ONLY
             else int(sys.argv[1]) if len(sys.argv) > 1 else 20)

if not ANALYSE_ONLY:
    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])
    # Even spacing through the inventory rather than a random draw: the
    # inventory is ordered by country, so this samples across countries by
    # construction and is reproducible without depending on a seed.
    inv = inv.iloc[:: max(len(inv) // N_WINDOWS, 1)].head(N_WINDOWS)

n_combos = (len(OBSERVATIONS) * len(TEMP_FORMS) * len(RAIN_LAGS)
            * len(TRAIN_FRACS) * len(STRUCTURES))

if not ANALYSE_ONLY:
    print(f"{len(inv)} windows x {n_combos} combinations x 2 models x "
          f"{len(ARMS)} optimiser settings\n")
    print("Loading case data...")
    cases_all = pd.read_csv(raw / BASE["data"]["opendengue"]["csv_name"],
                            usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                                     "calendar_start_date", "calendar_end_date",
                                     "dengue_total"], low_memory=False)
    cases_all["start"] = pd.to_datetime(cases_all["calendar_start_date"],
                                        errors="coerce")
    cases_all["end"] = pd.to_datetime(cases_all["calendar_end_date"],
                                      errors="coerce")
    cases_all = cases_all[
        (cases_all["end"] - cases_all["start"]).dt.days + 1 == 7]

_clim: dict[str, pd.DataFrame] = {}


def load_climate(slug):
    if slug not in _clim:
        path = clim_dir / f"{slug}.csv"
        if not path.exists():
            _clim[slug] = None
            return None
        with open(path, "r", encoding="utf-8") as fh:
            lines = fh.readlines()
        try:
            head = next(i for i, ln in enumerate(lines) if ln.startswith("YEAR,"))
        except StopIteration:
            _clim[slug] = None
            return None
        c = pd.read_csv(path, skiprows=head).replace(-999.0, np.nan)
        c["date"] = (pd.to_datetime(c["YEAR"].astype(str), format="%Y")
                     + pd.to_timedelta(c["DOY"] - 1, unit="D"))
        _clim[slug] = c.set_index("date").sort_index()
    return _clim[slug]


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
    pt = point_for(w["country"], w["unit"], w["level"])
    if pt is None:
        return None
    name, lat, lon = pt
    slug = f"{name.lower().replace(' ', '_').replace('/', '_')}_{round(lat, 4)}_{round(lon, 4)}"
    clim = load_climate(slug)
    if clim is None:
        return None
    s = case_series(w)
    if len(s) < 20:
        return None
    gaps = s.index.to_series().diff().dt.days.dropna()
    if not (gaps == 7).all():
        return None
    smooth = BASE["model"]["climate_forcing"]["rain_smooth_weeks"]
    c = clim.copy()
    c["rain_lagged"] = lagged_smoothed_rain(c["PRECTOTCORR"], lag_weeks, smooth)
    weekly = (c.resample("7D", origin=s.index.min())
               .mean(numeric_only=True).reindex(s.index))
    if weekly[["T2M", "rain_lagged"]].isna().any().any():
        return None
    z_rain, *_ = standardise(weekly["rain_lagged"].to_numpy())
    days = (s.index - s.index.min()).days.to_numpy(float)
    return Dataset(days, s.to_numpy(float), weekly["T2M"].to_numpy(float),
                   z_rain, 1e7, f"{w['country']}/{w['unit']}")


def cfg_for(observation, temp_form, lag, structure, restarts):
    c = copy.deepcopy(BASE)
    c["model"]["temperature_form"] = temp_form
    c["model"]["structure"] = structure
    c["model"]["climate_forcing"]["rain_lag_weeks"] = lag
    c["inference"]["observation"] = observation
    c["inference"]["classical"]["n_restarts"] = restarts
    return c


if not ANALYSE_ONLY:
    rows = []
    t0 = time.time()

    for idx, (_, w) in enumerate(inv.iterrows(), 1):
        label = f"{w['country']}/{w['unit']} {w['start'].date()}"
        data_by_lag = {lag: build(w, lag) for lag in RAIN_LAGS}
        data_by_lag = {k: v for k, v in data_by_lag.items() if v is not None}
        if not data_by_lag:
            continue
        ref = data_by_lag[min(data_by_lag)]

        # The warm start the main run used: one anchor fit per window.
        set_temperature_form("briere")
        anchor_cfg = cfg_for("nb", "briere", min(data_by_lag), "hostvector", 3)
        try:
            anchor = fit(ref, anchor_cfg, fixed, model="climate", observation="nb")
        except Exception as exc:
            print(f"[{idx:3d}/{len(inv)}] {label:40s} anchor failed: {exc}")
            continue
        warm, k_hat = dict(anchor.theta), float(anchor.nb_k or 20.0)

        for observation, temp_form, lag, frac, structure in itertools.product(
                OBSERVATIONS, TEMP_FORMS, RAIN_LAGS, TRAIN_FRACS, STRUCTURES):
            if lag not in data_by_lag:
                continue
            full = data_by_lag[lag]
            train = full if frac >= 1.0 else full.head(
                int(round(len(full.days) * frac)))
            set_temperature_form(temp_form)

            for arm, restarts, warm_start in ARMS:
                cfg = cfg_for(observation, temp_form, lag, structure, restarts)
                res = {}
                for model in ("climate", "constant"):
                    try:
                        res[model] = fit(
                            train, cfg, fixed, model=model, observation=observation,
                            start_from=warm if warm_start else None,
                            n_starts_override=restarts,
                            fixed_nb_k=k_hat if observation == "nb" else None)
                    except (RuntimeError, IntegrationFailure, ValueError):
                        res[model] = None
                if res["climate"] is None or res["constant"] is None:
                    continue
                d_aic = res["climate"].aic - res["constant"].aic
                rows.append(dict(
                    country=w["country"], unit=w["unit"],
                    window_start=w["start"].date(), arm=arm,
                    observation=observation, temp_form=temp_form, rain_lag=lag,
                    train_frac=frac, structure=structure,
                    aic_climate=round(float(res["climate"].aic), 4),
                    aic_constant=round(float(res["constant"].aic), 4),
                    delta_aic=round(float(d_aic), 4),
                    climate_wins=int(d_aic < 0),
                    weeks=len(train.days), cases=float(train.cases.sum())))

        pd.DataFrame(rows).to_csv(OUT, index=False)
        rate = (time.time() - t0) / idx
        print(f"[{idx:3d}/{len(inv)}] {label:40s} {rate:6.1f} s/window  "
              f"eta {rate * (len(inv) - idx) / 60:5.1f} min")


    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
else:
    d = pd.read_csv(OUT)
    print(f"analysing the stored table: {OUT.name}\n")

# ---------------------------------------------------------------------------
KEY = ["country", "unit", "window_start"]
FIT = KEY + ["observation", "temp_form", "rain_lag", "train_frac", "structure"]
wide = d.pivot_table(index=FIT, columns="arm",
                     values=["aic_climate", "aic_constant", "delta_aic",
                             "climate_wins"]).dropna()
print(f"\n{len(d):,} fits, {len(wide):,} comparisons present in both arms\n")

print("=" * 76)
print("1. DOES THE THOROUGH ARM FIND BETTER OPTIMA?")
print("=" * 76)
TOL = 1e-6  # below this an AIC difference is numerical, not a better optimum
for model in ("climate", "constant"):
    a = wide[(f"aic_{model}", "as_run")]
    t = wide[(f"aic_{model}", "thorough")]
    better = (t < a - TOL)
    gain = (a - t)[better]
    print(f"  {model:9s} thorough is better in {better.mean() * 100:5.1f}% of fits"
          f"   median gain {gain.median() if len(gain) else 0:8.3f} AIC"
          f"   max {gain.max() if len(gain) else 0:9.2f}")
print("\n  If these are near zero the objection is empty: there was nothing")
print("  left for a longer search to find.")

print("\n" + "=" * 76)
print("2. DOES THE VERDICT CHANGE? (same scale as the five analysis choices)")
print("=" * 76)
flip = (wide[("climate_wins", "as_run")]
        != wide[("climate_wins", "thorough")])
p_a = wide[("climate_wins", "as_run")].mean()
p_t = wide[("climate_wins", "thorough")].mean()
print(f"  optimiser setting 3 warm -> 10 cold: flipped {flip.mean() * 100:5.1f}%"
      f"   P(climate wins) {p_a:.3f} -> {p_t:.3f}")
print(f"  weakest genuine factor (model structure, step 17):    8.3%")
print(f"  strongest genuine factor (observation model, step 17): 38.1%")
shift = (wide[("delta_aic", "thorough")] - wide[("delta_aic", "as_run")])
print(f"\n  shift in delta AIC: median {shift.median():+.3f}, "
      f"IQR {shift.quantile(.25):+.3f} to {shift.quantile(.75):+.3f}")

print("\n" + "=" * 76)
print("3. DOES THE HEADLINE CHANGE?")
print("=" * 76)
for arm, _, _ in ARMS:
    sub = d[d["arm"] == arm].copy()
    sub, n_combos = complete_windows(sub, KEY)
    per = window_verdicts(sub, n_combos, KEY)
    print(f"  {arm:9s} unstable in {per['unstable'].mean() * 100:5.1f}% of "
          f"{len(per)} windows completing all {n_combos} combinations")
print("\n  This is the number the paper reports. It is the one that has to")
print("  survive a more thorough optimiser.")

print(f"\nTable: {OUT}")
