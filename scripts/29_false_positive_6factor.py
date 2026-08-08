"""Pipeline step 29 — error rates against known truths, on the six-factor design.

Supersedes step 21, which ran the five-factor design. The paper says the
simulation applies "the full factorial", and once the factorial grew to 144
combinations that stopped being true. Rather than weaken the sentence, the
simulation is brought back into step with the study.

The design is otherwise unchanged from step 21, including the correction that
made its power estimates usable: the two truths are matched on **mean
transmission**. Switching the climate coefficients on multiplies beta(t) by a
factor whose average is not one, so an unadjusted alternative differs from the
null in overall transmission as well as in whether transmission varies — and in
the first version that raised R0 enough to diverge the generating integration
for 21 of 59 windows, leaving power estimated only where a strong effect
happened not to break the model. Rescaling beta_0 by the mean of the climate
multiplier isolates the contrast the study is about.

What it measures, per observation model and decision rule:

* **False positives** — data generated with no climate effect, verdict says there
  is one.
* **Power** — data generated with an effect of realistic size, verdict finds it.
* **Separation** — power minus false positives. Reported because either rate
  alone is misleading: a rule with 95% power and 83% false positives is not
  detecting anything, it is agreeing with everything.
* **Whether the instability is manufactured by the size of the factorial** —
  the same design is applied to both truths, so any difference between them is
  about the data rather than about how many analyses were run.

Parallel across windows, on the same pattern as step 26: the parent does the
I/O, workers receive prepared Datasets. At 288 fits per window per truth this is
about an hour on twelve processes, against most of a day serially.
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
STRUCTURES = ("hostvector", "seir", "hostvector_tempmort")
PARAM_SETS = ("central", "alt")

TRUE_A_TEMP = 1.00
TRUE_A_RAIN = 0.30

BASE = load_config()
raw = resolve(BASE, "raw")
tables = resolve(BASE, "tables")
clim_dir = raw / "climate_global"
OUT = tables / "28_false_positive_6factor.csv"

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


def fit_window(package):
    """Both truths, all 144 combinations, for one window."""
    meta, data_by_lag, seed = package
    label = f"{meta['country']}/{meta['unit']} {meta['window_start']}"
    if not data_by_lag:
        return label, [], "no usable lag"
    rng = np.random.default_rng(seed)
    ref = data_by_lag[min(data_by_lag)]

    set_temperature_form("briere")
    base_cfg = cfg_for("nb", "briere", min(data_by_lag), "hostvector")
    try:
        anchor = fit(ref, base_cfg, FIXED["central"], model="constant",
                     observation="nb")
    except Exception as exc:                        # noqa: BLE001
        return label, [], f"anchor failed: {exc}"
    k_hat = float(anchor.nb_k or 20.0)

    # Labelled "no_effect" rather than "null": pandas reads the string "null"
    # back as NaN, which silently deleted an entire arm of this study once.
    truth_null = dict(anchor.theta)
    truth_null.update(a_temp=1e-9, a_rain=0.0)
    truth_climate = dict(anchor.theta)
    truth_climate.update(a_temp=TRUE_A_TEMP, a_rain=TRUE_A_RAIN)
    mult = make_forcing({**truth_climate, "beta_0": 1.0}, ref, "climate",
                        base_cfg).on_grid(ref.days)
    truth_climate["beta_0"] = float(anchor.theta["beta_0"]
                                    / max(float(np.mean(mult)), 1e-12))

    rows = []
    for truth_name, truth in (("no_effect", truth_null),
                              ("climate", truth_climate)):
        try:
            mu_true = predict(truth, ref, FIXED["central"], "climate", base_cfg)
        except IntegrationFailure:
            continue
        # Gamma-Poisson mixture reproduces NB2 with the window's own dispersion.
        lam = rng.gamma(shape=k_hat, scale=np.maximum(mu_true, 1e-12) / k_hat)
        synth = rng.poisson(lam).astype(float)
        if synth.sum() < 100:
            continue

        for observation, temp_form, lag, frac, structure, pset in \
                itertools.product(OBSERVATIONS, TEMP_FORMS, RAIN_LAGS,
                                  TRAIN_FRACS, STRUCTURES, PARAM_SETS):
            if lag not in data_by_lag:
                continue
            d0 = data_by_lag[lag]
            full = Dataset(d0.days, synth, d0.temp_c, d0.z_rain,
                           d0.population, d0.label)
            train = full if frac >= 1.0 else full.head(
                int(round(len(full.days) * frac)))
            fixed = FIXED[pset]
            cfg = cfg_for(observation, temp_form, lag, structure)
            set_temperature_form(temp_form)

            res = {}
            for model in ("climate", "constant"):
                try:
                    res[model] = fit(train, cfg, fixed, model=model,
                                     observation=observation, start_from=truth,
                                     n_starts_override=3,
                                     fixed_nb_k=k_hat if observation == "nb" else None)
                except (RuntimeError, IntegrationFailure, ValueError):
                    res[model] = None
            if res["climate"] is None or res["constant"] is None:
                continue
            d_aic = res["climate"].aic - res["constant"].aic
            rows.append(dict(
                country=meta["country"], unit=meta["unit"],
                window_start=meta["window_start"], truth=truth_name,
                observation=observation, temp_form=temp_form, rain_lag=lag,
                train_frac=frac, structure=structure, params=pset,
                delta_aic=round(float(d_aic), 3),
                climate_wins=int(d_aic < 0),
                climate_wins_margin4=int(d_aic < -4)))
    return label, rows, None


def build_packages(n_windows: int):
    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])
    inv = inv.iloc[:: max(len(inv) // n_windows, 1)].head(n_windows)

    print(f"{len(inv)} windows x 2 truths x {N_COMBOS} combinations x 2 models"
          f" = {len(inv) * 2 * N_COMBOS * 2:,} fits\n", flush=True)
    print("Loading case data...", flush=True)
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

    smooth = BASE["model"]["climate_forcing"]["rain_smooth_weeks"]
    packages = []
    for idx, (_, w) in enumerate(inv.iterrows()):
        pt = point_for(w["country"], w["unit"], w["level"])
        if pt is None:
            continue
        name, lat, lon = pt
        slug = (f"{name.lower().replace(' ', '_').replace('/', '_')}"
                f"_{round(lat, 4)}_{round(lon, 4)}")
        clim = climate(slug)
        if clim is None:
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
        # A distinct seed per window, derived from the configured one, so the
        # synthetic data are reproducible and independent across processes.
        packages.append(({"country": w["country"], "unit": w["unit"],
                          "window_start": str(w["start"].date())},
                         by_lag, BASE["seed"] + idx))
    return packages


def report(d: pd.DataFrame) -> None:
    KEY = ["country", "unit", "window_start"]
    MARGINS = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24]
    null = d[d["truth"] == "no_effect"]
    alt = d[d["truth"] == "climate"]

    rows = []
    for obs in ("poisson", "nb"):
        n_sub, a_sub = null[null.observation == obs], alt[alt.observation == obs]
        for m in MARGINS:
            fp = float((n_sub["delta_aic"] < -m).mean())
            pw = float((a_sub["delta_aic"] < -m).mean())
            rows.append(dict(observation=obs, margin=m,
                             false_positive_pct=round(fp * 100, 1),
                             power_pct=round(pw * 100, 1),
                             youden=round((pw - fp) * 100, 1)))
    oc = pd.DataFrame(rows)
    oc.to_csv(tables / "29_operating_characteristics_6factor.csv", index=False)

    print("=" * 74)
    print("OPERATING CHARACTERISTICS")
    print("=" * 74)
    for obs in ("poisson", "nb"):
        sub = oc[oc.observation == obs]
        best = sub.loc[sub["youden"].idxmax()]
        print(f"\n  {obs}:")
        for _, r in sub.iterrows():
            mark = "  <-- best" if r["margin"] == best["margin"] else ""
            print(f"    margin {r['margin']:2.0f}  FP {r['false_positive_pct']:5.1f}%"
                  f"  power {r['power_pct']:5.1f}%"
                  f"  separation {r['youden']:5.1f}{mark}")

    print("\n" + "=" * 74)
    print("IS THE INSTABILITY AN ARTEFACT OF RUNNING THE FACTORIAL?")
    print("=" * 74)
    for truth, label in (("climate", "a real climate effect"),
                         ("no_effect", "no climate effect")):
        sub = d[d["truth"] == truth]
        if sub.empty:
            continue
        counts = sub.groupby(KEY).size()
        n = int(counts.max())
        full = sub.set_index(KEY).loc[counts[counts == n].index].reset_index()
        for margin, name in ((0, "sign of dAIC"), (4, "abstain within 4")):
            uns = dec = 0
            for _, g in full.groupby(KEY):
                v = np.where(g.delta_aic < -margin, 1,
                             np.where(g.delta_aic > margin, -1, 0))
                spoken = v[v != 0]
                if len(spoken) == 0:
                    continue
                dec += 1
                if len(set(spoken.tolist())) > 1:
                    uns += 1
            print(f"  {label:24s} {name:18s} unstable {uns / dec * 100:5.1f}%"
                  f"   answers {dec / (counts == n).sum() * 100:5.1f}%")
    print("\n  Same factorial in both rows; only the truth differs.")


def main():
    n_windows = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    packages = build_packages(n_windows)
    n_workers = max(1, min(12, (mp.cpu_count() or 2) - 2))
    print(f"{len(packages)} windows prepared; {n_workers} workers\n", flush=True)

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
            print(f"[{i:3d}/{len(packages)}] {label:40s} {len(got):4d} rows  "
                  f"eta {rate * (len(packages) - i) / 60:5.1f} min", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    print(f"\n{len(d):,} fits in {(time.time() - t0) / 60:.1f} min\n")
    report(d)
    print(f"\nTable: {OUT}")


if __name__ == "__main__":
    main()
