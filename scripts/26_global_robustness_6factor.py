"""Pipeline step 26 — the six-factor factorial, run in parallel.

Two objections to the five-factor study were left standing in its own
limitations section, and both are the kind a reviewer raises first because the
paper raised them itself.

**"Your model is too simple."** Every fit held mosquito lifespan constant across
seasons whose temperature swings by fifteen degrees. Applied dengue modelling
does not do that. A third structure is added here: the same host--vector model
with adult mosquito mortality following a thermal response fixed from vector
biology. It estimates no extra parameter, so the comparison remains between
structures rather than between one model given more freedom than another.

**"You fixed the parameters that matter."** Mosquito lifespan, the incubation
periods and the vector-to-host ratio were fixed from the literature and never
varied, while the study criticised others for exactly this class of unexamined
choice. Published ranges are wide — viraemia 4-7 days, incubation 4-7 days,
lifespan 8-15 days, vector-to-host ratio from below one to above ten — so
choosing within them is a degree of freedom like the rainfall lag. A second
parameter set, at the other end of the same ranges, is now a factor.

The design is therefore:

    observation  x2   negative binomial, Poisson
    temp_form    x2   Briere, log-linear
    rain_lag     x3   3, 5, 7 weeks
    train_frac   x2   all, first 75%
    structure    x3   host-vector, human-only SEIR, host-vector + thermal mortality
    params       x2   central literature set, alternative literature set
    ----------------
                144 combinations, two models each: 288 fits per window

At roughly 1.2 s per fit that is 22 hours single-threaded for the full
inventory, so the windows are distributed across processes. The parent does all
of the I/O — reading OpenDengue once, building each window's datasets — and
workers receive only prepared Datasets and return rows. Workers therefore hold
no copy of the 2.2-million-row case table, which is what makes running twelve of
them possible on an ordinary laptop.

Results are written after every window completes, so an interrupted run leaves a
usable table rather than nothing.

**What the null model means in the third structure, stated plainly.** Under
``hostvector_tempmort`` temperature enters *both* compared models, because it
drives mosquito mortality in the structure itself. The constant model there is
therefore not climate-free: it is "climate affects mosquito survival but not the
transmission coefficient". The comparison being made is the same one in all
three structures --- does letting the transmission coefficient depend on climate
improve the fit, holding the rest of the mechanism fixed --- but the baseline it
is measured against is richer in the third.

That is deliberate and it reflects practice: mechanistic dengue models routinely
carry temperature in their vector components while still testing whether a
climate term on transmission earns its parameters. It also makes the third
structure the hardest test of the three, since its null already contains thermal
biology. If climate forcing is still endorsed there at least as often, the
endorsement is not coming from temperature being absent from the alternative.
A reader must be told this, because a structure whose null is not the same null
would otherwise look like a like-for-like comparison and is not.
"""

from __future__ import annotations

import copy
import itertools
import multiprocessing as mp
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.climate import (ThermalMortality, lagged_smoothed_rain,  # noqa: E402
                               standardise)
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
STRUCTURES = ("hostvector", "seir", "hostvector_tempmort")
PARAM_SETS = ("central", "alt")

BASE = load_config()
raw = resolve(BASE, "raw")
tables = resolve(BASE, "tables")
clim_dir = raw / "climate_global"
OUT = tables / "25_global_robustness_6factor.csv"

FIXED = {"central": FixedParams.from_config(BASE, "fixed"),
         "alt": FixedParams.from_config(BASE, "fixed_alt")}

N_COMBOS = (len(OBSERVATIONS) * len(TEMP_FORMS) * len(RAIN_LAGS)
            * len(TRAIN_FRACS) * len(STRUCTURES) * len(PARAM_SETS))


def cfg_for(observation, temp_form, lag, structure):
    c = copy.deepcopy(BASE)
    c["model"]["temperature_form"] = temp_form
    c["model"]["structure"] = structure
    c["model"]["climate_forcing"]["rain_lag_weeks"] = lag
    c["inference"]["observation"] = observation
    c["inference"]["classical"]["n_restarts"] = 3
    return c


# ---------------------------------------------------------------------------
# Worker: fits one window under every combination. Receives prepared Datasets,
# never touches the filesystem.
# ---------------------------------------------------------------------------
def fit_window(package):
    meta, data_by_lag = package
    label = f"{meta['country']}/{meta['unit']} {meta['window_start']}"
    if not data_by_lag:
        return label, [], "no usable lag"

    ref = data_by_lag[min(data_by_lag)]
    # One anchor fit per window supplies the dispersion and the warm start, as
    # in step 16. It uses the central parameter set; using a different anchor
    # per parameter set would confound the factor with its own starting point.
    set_temperature_form("briere")
    anchor_cfg = cfg_for("nb", "briere", min(data_by_lag), "hostvector")
    try:
        anchor = fit(ref, anchor_cfg, FIXED["central"], model="climate",
                     observation="nb")
    except Exception as exc:                       # noqa: BLE001
        return label, [], f"anchor failed: {exc}"
    warm, k_hat = dict(anchor.theta), float(anchor.nb_k or 20.0)

    rows = []
    for observation, temp_form, lag, frac, structure, pset in itertools.product(
            OBSERVATIONS, TEMP_FORMS, RAIN_LAGS, TRAIN_FRACS, STRUCTURES,
            PARAM_SETS):
        if lag not in data_by_lag:
            continue
        full = data_by_lag[lag]
        train = full if frac >= 1.0 else full.head(
            int(round(len(full.days) * frac)))
        fixed = FIXED[pset]
        cfg = cfg_for(observation, temp_form, lag, structure)
        set_temperature_form(temp_form)

        res = {}
        for model in ("climate", "constant"):
            try:
                res[model] = fit(train, cfg, fixed, model=model,
                                 observation=observation, start_from=warm,
                                 n_starts_override=3,
                                 fixed_nb_k=k_hat if observation == "nb" else None)
            except (RuntimeError, IntegrationFailure, ValueError):
                res[model] = None
        if res["climate"] is None or res["constant"] is None:
            continue

        # R0 for the thermal-mortality structure uses the time-averaged
        # mortality the model actually experienced, not the constant it was
        # scaled from.
        mu_eff = (ThermalMortality(train.days, train.temp_c, fixed.mu_v).mean()
                  if structure == "hostvector_tempmort" else None)
        r0_struct = "seir" if structure == "seir" else "hostvector"
        d_aic = res["climate"].aic - res["constant"].aic

        # Held-out deviance on the withheld tail, where a tail was withheld.
        # Without this the out-of-sample comparison — an entire results section,
        # and the one that shows validation halves endorsement without restoring
        # stability — cannot be recomputed from the table.
        held = np.nan
        if frac < 1.0:
            tail = slice(len(train.days), None)
            try:
                dev = {}
                for model in ("climate", "constant"):
                    mu = predict(res[model].theta, full, fixed, model, cfg)
                    resid = (nb_deviance_residuals(full.cases[tail], mu[tail], k_hat)
                             if observation == "nb"
                             else poisson_deviance_residuals(full.cases[tail],
                                                             mu[tail]))
                    dev[model] = float(np.sum(resid ** 2))
                held = dev["climate"] - dev["constant"]
            except (IntegrationFailure, ValueError):
                held = np.nan
        rows.append(dict(
            country=meta["country"], unit=meta["unit"], level=meta["level"],
            window_start=meta["window_start"], weeks=len(train.days),
            cases=float(train.cases.sum()),
            observation=observation, temp_form=temp_form, rain_lag=lag,
            train_frac=frac, structure=structure, params=pset,
            delta_aic=round(float(d_aic), 3),
            climate_wins=bool(d_aic < 0),
            heldout_delta=round(float(held), 3) if np.isfinite(held) else np.nan,
            R0_climate=round(basic_reproduction_number(
                res["climate"].theta["beta_0"], fixed, r0_struct, mu_eff), 4),
            R0_constant=round(basic_reproduction_number(
                res["constant"].theta["beta_0"], fixed, r0_struct, mu_eff), 4),
            a_temp=round(float(res["climate"].theta["a_temp"]), 4),
            a_rain=round(float(res["climate"].theta["a_rain"]), 4)))
    return label, rows, None


# ---------------------------------------------------------------------------
# Parent: all I/O happens here.
# ---------------------------------------------------------------------------
def build_packages():
    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    if limit:
        inv = inv.iloc[:: max(len(inv) // limit, 1)].head(limit)

    print(f"{len(inv)} windows x {N_COMBOS} combinations x 2 models "
          f"= {len(inv) * N_COMBOS * 2:,} fits\n")
    print("Loading case data...")
    cases = pd.read_csv(raw / BASE["data"]["opendengue"]["csv_name"],
                        usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                                 "calendar_start_date", "calendar_end_date",
                                 "dengue_total"], low_memory=False)
    cases["start"] = pd.to_datetime(cases["calendar_start_date"], errors="coerce")
    cases["end"] = pd.to_datetime(cases["calendar_end_date"], errors="coerce")
    cases = cases[(cases["end"] - cases["start"]).dt.days + 1 == 7]

    clim_cache: dict[str, pd.DataFrame | None] = {}

    def climate(slug):
        if slug not in clim_cache:
            path = clim_dir / f"{slug}.csv"
            if not path.exists():
                clim_cache[slug] = None
                return None
            with open(path, "r", encoding="utf-8") as fh:
                lines = fh.readlines()
            try:
                head = next(i for i, ln in enumerate(lines)
                            if ln.startswith("YEAR,"))
            except StopIteration:
                clim_cache[slug] = None
                return None
            c = pd.read_csv(path, skiprows=head).replace(-999.0, np.nan)
            c["date"] = (pd.to_datetime(c["YEAR"].astype(str), format="%Y")
                         + pd.to_timedelta(c["DOY"] - 1, unit="D"))
            clim_cache[slug] = c.set_index("date").sort_index()
        return clim_cache[slug]

    def series(w):
        sub = cases[cases["adm_0_name"].astype(str).str.upper()
                    == str(w["country"]).upper()]
        if w["level"] == "national":
            sub = sub[sub["adm_1_name"].isna()]
        else:
            sub = sub[sub["adm_1_name"].astype(str).str.upper()
                      == str(w["unit"]).upper()]
            sub = sub[sub["adm_2_name"].isna()]
        sub = sub[(sub["start"] >= w["start"]) & (sub["start"] <= w["end"])]
        return sub.groupby("start")["dengue_total"].sum().sort_index()

    smooth = BASE["model"]["climate_forcing"]["rain_smooth_weeks"]
    packages = []
    for _, w in inv.iterrows():
        pt = point_for(w["country"], w["unit"], w["level"])
        if pt is None:
            continue
        name, lat, lon = pt
        slug = (f"{name.lower().replace(' ', '_').replace('/', '_')}"
                f"_{round(lat, 4)}_{round(lon, 4)}")
        clim = climate(slug)
        if clim is None:
            continue
        s = series(w)
        if len(s) < 20:
            continue
        gaps = s.index.to_series().diff().dt.days.dropna()
        if not (gaps == 7).all():
            continue

        by_lag = {}
        for lag in RAIN_LAGS:
            c = clim.copy()
            c["rain_lagged"] = lagged_smoothed_rain(c["PRECTOTCORR"], lag, smooth)
            weekly = (c.resample("7D", origin=s.index.min())
                       .mean(numeric_only=True).reindex(s.index))
            if weekly[["T2M", "rain_lagged"]].isna().any().any():
                continue
            z_rain, *_ = standardise(weekly["rain_lagged"].to_numpy())
            days = (s.index - s.index.min()).days.to_numpy(float)
            by_lag[lag] = Dataset(days, s.to_numpy(float),
                                  weekly["T2M"].to_numpy(float), z_rain, 1e7,
                                  f"{w['country']}/{w['unit']}")
        if not by_lag:
            continue
        packages.append(({"country": w["country"], "unit": w["unit"],
                          "level": w["level"],
                          "window_start": str(w["start"].date())}, by_lag))
    return packages


def main():
    packages = build_packages()
    n_workers = max(1, min(12, (mp.cpu_count() or 2) - 2))
    print(f"{len(packages)} windows prepared; {n_workers} worker processes\n")

    rows, skipped, t0 = [], [], time.time()
    with mp.Pool(n_workers) as pool:
        for i, (label, got, why) in enumerate(
                pool.imap_unordered(fit_window, packages), 1):
            if why:
                skipped.append((label, why))
                print(f"[{i:3d}/{len(packages)}] {label:42s} SKIPPED: {why}")
                continue
            rows.extend(got)
            pd.DataFrame(rows).to_csv(OUT, index=False)
            rate = (time.time() - t0) / i
            print(f"[{i:3d}/{len(packages)}] {label:42s} {len(got):3d} combos  "
                  f"{rate:5.1f} s/window  eta "
                  f"{rate * (len(packages) - i) / 60:5.1f} min", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    print(f"\n{len(d):,} fits from {d.groupby(['country', 'unit', 'window_start']).ngroups} "
          f"windows in {(time.time() - t0) / 60:.1f} min")
    if skipped:
        print(f"\n{len(skipped)} windows skipped:")
        for label, why in skipped:
            print(f"  {label}: {why}")
    print(f"\nTable: {OUT}")


if __name__ == "__main__":
    main()
