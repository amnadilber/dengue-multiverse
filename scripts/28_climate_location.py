"""Pipeline step 28 — does the choice of climate grid cell change the verdict?

This is the degree of freedom the study named in its own limitations and did not
measure: one representative point supplies the climate for a whole country or
province. Step 27 downloaded a second point one degree of latitude away, about
110 km, well inside the area any of these case series aggregates over. A verdict
that changes at that distance is changing for a reason no analyst could defend
as substantive.

The design mirrors the factorial exactly rather than inventing a new comparison,
so the answer lands on the same scale as the other six factors: a sample of
windows is fitted under all 144 combinations at each of the two locations, and
the reported quantity is the share of paired comparisons whose verdict differs.
That number can be read directly against the 38% for the observation model and
the 8% for model structure.

Run on a subsample rather than the full inventory. The cost is 288 combinations
per window against 144, and the question is whether location belongs in the same
league as the other factors, which a few dozen windows settle.

**How to read the result, and a caveat that cuts against it.** The offset is
mechanical: one degree of latitude, applied identically everywhere, with no
regard for terrain. In flat units it lands somewhere an analyst might genuinely
have chosen — the two series differ by a few tenths of a degree. In others it
crosses high ground: for Costa Rica the two points differ by 4.3 C in mean
temperature, which is not a choice between two defensible stations for the same
population so much as a choice between two climates. The flip rate should
therefore be read as an **upper bound** on what a sensible alternative would do,
not as an estimate of it. Reported that way, because the alternative — hand-picking
a second city per unit — would substitute a gazetteer of our own choices for a
rule anyone can reproduce, and this study is in no position to do that quietly.
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
from dengue_pk.climate import lagged_smoothed_rain, standardise  # noqa: E402
from dengue_pk.inference import Dataset, fit, set_temperature_form  # noqa: E402
from dengue_pk.locations import offset_point, point_for  # noqa: E402
from dengue_pk.models import FixedParams, IntegrationFailure  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OBSERVATIONS = ("nb", "poisson")
TEMP_FORMS = ("briere", "loglinear")
RAIN_LAGS = (3, 5, 7)
TRAIN_FRACS = (1.0, 0.75)
STRUCTURES = ("hostvector", "seir", "hostvector_tempmort")
PARAM_SETS = ("central", "alt")
LOCATIONS = ("primary", "offset")

BASE = load_config()
raw = resolve(BASE, "raw")
tables = resolve(BASE, "tables")
clim_dir = raw / "climate_global"
OUT = tables / "26_climate_location.csv"

FIXED = {"central": FixedParams.from_config(BASE, "fixed"),
         "alt": FixedParams.from_config(BASE, "fixed_alt")}


def cfg_for(observation, temp_form, lag, structure):
    c = copy.deepcopy(BASE)
    c["model"]["temperature_form"] = temp_form
    c["model"]["structure"] = structure
    c["model"]["climate_forcing"]["rain_lag_weeks"] = lag
    c["inference"]["observation"] = observation
    c["inference"]["classical"]["n_restarts"] = 3
    return c


def fit_window(package):
    """All 144 combinations at both locations, for one window."""
    meta, by_loc = package
    label = f"{meta['country']}/{meta['unit']} {meta['window_start']}"
    if len(by_loc) < 2:
        return label, [], "second location unavailable"

    ref = by_loc["primary"][min(by_loc["primary"])]
    set_temperature_form("briere")
    anchor_cfg = cfg_for("nb", "briere", min(by_loc["primary"]), "hostvector")
    try:
        anchor = fit(ref, anchor_cfg, FIXED["central"], model="climate",
                     observation="nb")
    except Exception as exc:                        # noqa: BLE001
        return label, [], f"anchor failed: {exc}"
    warm, k_hat = dict(anchor.theta), float(anchor.nb_k or 20.0)

    rows = []
    for loc, observation, temp_form, lag, frac, structure, pset in \
            itertools.product(LOCATIONS, OBSERVATIONS, TEMP_FORMS, RAIN_LAGS,
                              TRAIN_FRACS, STRUCTURES, PARAM_SETS):
        data_by_lag = by_loc[loc]
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
        d_aic = res["climate"].aic - res["constant"].aic
        rows.append(dict(
            country=meta["country"], unit=meta["unit"],
            window_start=meta["window_start"], location=loc,
            observation=observation, temp_form=temp_form, rain_lag=lag,
            train_frac=frac, structure=structure, params=pset,
            delta_aic=round(float(d_aic), 3), climate_wins=int(d_aic < 0),
            weeks=len(train.days), cases=float(train.cases.sum())))
    return label, rows, None


def read_climate(slug):
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
    return c.set_index("date").sort_index()


def build_packages(n_windows: int):
    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])
    inv = inv.iloc[:: max(len(inv) // n_windows, 1)].head(n_windows)

    print("Loading case data...", flush=True)
    cases = pd.read_csv(raw / BASE["data"]["opendengue"]["csv_name"],
                        usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                                 "calendar_start_date", "calendar_end_date",
                                 "dengue_total"], low_memory=False)
    cases["start"] = pd.to_datetime(cases["calendar_start_date"], errors="coerce")
    cases["end"] = pd.to_datetime(cases["calendar_end_date"], errors="coerce")
    cases = cases[(cases["end"] - cases["start"]).dt.days + 1 == 7]

    smooth = BASE["model"]["climate_forcing"]["rain_smooth_weeks"]
    packages = []
    for _, w in inv.iterrows():
        pt = point_for(w["country"], w["unit"], w["level"])
        if pt is None:
            continue
        name, lat, lon = pt
        stem = name.lower().replace(" ", "_").replace("/", "_")
        alt_lat, alt_lon = offset_point(lat, lon)
        clim = {"primary": read_climate(f"{stem}_{round(lat, 4)}_{round(lon, 4)}"),
                "offset": read_climate(f"{stem}_alt_{alt_lat}_{alt_lon}")}
        if clim["primary"] is None or clim["offset"] is None:
            continue

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
        if len(s) < 20:
            continue
        gaps = s.index.to_series().diff().dt.days.dropna()
        if not (gaps == 7).all():
            continue

        by_loc = {}
        for loc, c0 in clim.items():
            by_lag = {}
            for lag in RAIN_LAGS:
                c = c0.copy()
                c["rain_lagged"] = lagged_smoothed_rain(c["PRECTOTCORR"],
                                                        lag, smooth)
                weekly = (c.resample("7D", origin=s.index.min())
                           .mean(numeric_only=True).reindex(s.index))
                if weekly[["T2M", "rain_lagged"]].isna().any().any():
                    continue
                z_rain, *_ = standardise(weekly["rain_lagged"].to_numpy())
                days = (s.index - s.index.min()).days.to_numpy(float)
                by_lag[lag] = Dataset(days, s.to_numpy(float),
                                      weekly["T2M"].to_numpy(float), z_rain,
                                      1e7, f"{w['country']}/{w['unit']}")
            if by_lag:
                by_loc[loc] = by_lag
        if len(by_loc) == 2:
            packages.append(({"country": w["country"], "unit": w["unit"],
                              "window_start": str(w["start"].date())}, by_loc))
    return packages


def report(d: pd.DataFrame) -> None:
    KEY = ["country", "unit", "window_start"]
    CELL = ["observation", "temp_form", "rain_lag", "train_frac", "structure",
            "params"]
    wide = d.pivot_table(index=KEY + CELL, columns="location",
                         values=["climate_wins", "delta_aic"]).dropna()
    if wide.empty:
        print("no paired comparisons")
        return

    flip = (wide[("climate_wins", "primary")]
            != wide[("climate_wins", "offset")])
    p_pri = wide[("climate_wins", "primary")].mean()
    p_off = wide[("climate_wins", "offset")].mean()
    shift = wide[("delta_aic", "offset")] - wide[("delta_aic", "primary")]

    print("=" * 74)
    print("DOES THE CLIMATE GRID CELL CHANGE THE VERDICT?")
    print("=" * 74)
    print(f"  {len(wide):,} paired comparisons across "
          f"{d.groupby(KEY).ngroups} outbreaks\n")
    print(f"  primary -> offset (110 km): flipped {flip.mean() * 100:5.1f}%"
          f"   P(climate wins) {p_pri:.3f} -> {p_off:.3f}")
    print(f"  shift in delta AIC: median {shift.median():+.2f}, "
          f"IQR {shift.quantile(.25):+.2f} to {shift.quantile(.75):+.2f}")
    # Read the comparison from the stored sensitivity table rather than writing
    # the numbers here: hard-coded ones went stale the moment the factorial grew,
    # and this is the line a reader uses to place the location factor.
    sens_path = tables / "15_choice_sensitivity.csv"
    if sens_path.exists():
        sens = pd.read_csv(sens_path).sort_values("flipped", ascending=False)
        print("\n  For comparison, every factor in the current factorial:")
        for _, r in sens.iterrows():
            mark = "  <-- climate location sits here" if (
                r["flipped"] < flip.mean()
                and (sens["flipped"] > flip.mean()).sum()
                == sens.index.get_loc(r.name)) else ""
            print(f"    {r['factor']:12s} {r['comparison']:38s} "
                  f"{r['flipped'] * 100:5.1f}%{mark}")

    print("\n" + "=" * 74)
    print("PER-WINDOW: DOES EITHER LOCATION GIVE A STABLE ANSWER?")
    print("=" * 74)
    for loc in LOCATIONS:
        sub = d[d["location"] == loc]
        counts = sub.groupby(KEY)["climate_wins"].agg(["sum", "size"])
        counts = counts[counts["size"] == counts["size"].max()]
        unstable = ((counts["sum"] >= 1)
                    & (counts["sum"] <= counts["size"] - 1)).mean()
        print(f"  {loc:8s} unstable in {unstable * 100:5.1f}% of "
              f"{len(counts)} windows")

    out = pd.DataFrame([dict(
        factor="climate_location", comparison="primary -> offset (110 km)",
        pairs=len(wide), p_climate_base=round(float(p_pri), 3),
        p_climate_alt=round(float(p_off), 3),
        flipped=round(float(flip.mean()), 3),
        net_change=round(float(p_off - p_pri), 3))])
    out.to_csv(tables / "27_location_sensitivity.csv", index=False)
    print(f"\nTable: {tables / '27_location_sensitivity.csv'}")


def main():
    n_windows = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    packages = build_packages(n_windows)
    n_workers = max(1, min(12, (mp.cpu_count() or 2) - 2))
    print(f"{len(packages)} windows with both locations; "
          f"{n_workers} workers\n", flush=True)

    rows, t0 = [], time.time()
    with mp.Pool(n_workers) as pool:
        for i, (label, got, why) in enumerate(
                pool.imap_unordered(fit_window, packages), 1):
            if why:
                print(f"[{i:3d}/{len(packages)}] {label:40s} SKIPPED: {why}",
                      flush=True)
                continue
            rows.extend(got)
            pd.DataFrame(rows).to_csv(OUT, index=False)
            rate = (time.time() - t0) / i
            print(f"[{i:3d}/{len(packages)}] {label:40s} {len(got):3d} rows  "
                  f"eta {rate * (len(packages) - i) / 60:5.1f} min", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    print(f"\n{len(d):,} fits in {(time.time() - t0) / 60:.1f} min\n")
    report(d)


if __name__ == "__main__":
    main()
