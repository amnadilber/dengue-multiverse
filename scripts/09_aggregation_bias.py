"""
Pipeline step 9 — does spatial aggregation bias the recovered transmission rate?

This is the project's second research question and it has not been answered until
now. Every fit so far has treated an administrative unit as one well-mixed
population, but a province is not well mixed: its districts have separate
epidemics that peak at different times, and summing them produces a curve broader
and flatter than any single outbreak. Fitting a homogeneous model to that sum must
distort something. The question is what, and by how much.

The 2021 Sindh season permits a direct test. The province and several of its
districts report over the same nineteen weeks, so the same epidemic can be fitted
at two spatial resolutions and the results compared. If aggregation were harmless,
the province-level R0 would sit among the district-level values.

One convenience makes this test unusually clean: R0 does not depend on the assumed
population at all. Section `docs/MODEL.md` shows the population enters only as a
product with the reporting fraction, confirmed numerically to four decimal places.
District census populations are therefore only needed to report a catchment size,
and any error in them cannot touch the quantity being compared.
"""

from __future__ import annotations

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

cfg = load_config()
OBS = cfg["inference"]["observation"]
fixed = FixedParams.from_config(cfg)
raw = resolve(cfg, "raw")
figures = resolve(cfg, "figures")
tables = resolve(cfg, "tables")

# 2017 census populations, in millions. Approximate, and it does not matter:
# R0 is invariant to the population, so these affect only the reported catchment.
DISTRICT_POP_M = {
    "THARPARKAR": 1.65, "KARACHI MALIR": 2.01, "KARACHI EAST": 2.91,
    "HYDERABAD": 2.20, "SHIKARPUR": 1.23, "LARKANA": 1.52,
    "HARIPUR": 1.00, "LAKKI MARAWAT": 0.88, "KOHAT": 0.99, "ABBOTTABAD": 1.33,
}

PROVINCES = {"SINDH": ("sindh_2021", cfg["windows"]["sindh_2021"]),
             "KHYBER PAKHTUNKHWA": ("kp_2021", cfg["windows"]["kp_2021"])}
MIN_TOTAL_CASES = 100          # below this a wave cannot constrain a fit
MIN_WEEKS = 15                 # a district reporting briefly cannot be compared

# District name variants in the OpenDengue Pakistan subset. These are the same
# places entered under different spellings, and they matter more than they look:
# "HYDERBAD" appears in exactly one week carrying 107 cases, and because the
# comparison below requires a window shared by every unit, that single stray
# record collapsed the shared window to one week and silently produced an empty
# analysis. Data-entry variants are not a cosmetic problem in aggregation work.
NAME_VARIANTS = {
    "HYDERBAD": "HYDERABAD",
    "LAKKIMARWAT": "LAKKI MARAWAT", "LAKIMARWAT": "LAKKI MARAWAT",
    "QAMBER": "KAMBER",
    "KARWEST": "KARACHI WEST",
}


def canonical(name: str) -> str:
    n = " ".join(str(name).upper().split())
    return NAME_VARIANTS.get(n, n)


def load_climate(path: Path) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("YEAR,"))
    clim = pd.read_csv(path, skiprows=start).replace(-999.0, np.nan)
    clim["date"] = (pd.to_datetime(clim["YEAR"].astype(str), format="%Y")
                    + pd.to_timedelta(clim["DOY"] - 1, unit="D"))
    return clim.set_index("date").sort_index()


def build_dataset(weeks: pd.Series, cases: np.ndarray, clim: pd.DataFrame,
                  population_m: float, label: str) -> Dataset:
    fc = cfg["model"]["climate_forcing"]
    clim = clim.copy()
    clim["rain_lagged"] = lagged_smoothed_rain(
        clim["PRECTOTCORR"], fc["rain_lag_weeks"], fc["rain_smooth_weeks"])
    weekly = (clim.resample("7D", origin=weeks.min())
                  .mean(numeric_only=True).reindex(weeks))
    z_rain, *_ = standardise(weekly["rain_lagged"].to_numpy())
    days = (weeks - weeks.min()).dt.days.to_numpy(float)
    return Dataset(days, cases.astype(float), weekly["T2M"].to_numpy(float),
                   z_rain, population_m * 1e6, label)


# --- Extract every district series for the 2021 season ----------------------
pk = pd.read_csv(raw / cfg["data"]["opendengue"]["csv_name"],
                 usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                          "calendar_start_date", "calendar_end_date",
                          "dengue_total"], low_memory=False)
pk = pk[pk["adm_0_name"].astype(str).str.upper() == "PAKISTAN"].copy()
pk["start"] = pd.to_datetime(pk["calendar_start_date"], errors="coerce")
pk["end"] = pd.to_datetime(pk["calendar_end_date"], errors="coerce")
pk = pk[(pk["end"] - pk["start"]).dt.days + 1 == 7]

rows, panels = [], []

for province, (window_key, window) in PROVINCES.items():
    print(f"\n{'=' * 74}\n{province} 2021\n{'=' * 74}")
    pt = window["climate_point"]
    clim = load_climate(raw / f"climate_{window_key}_{pt['name'].lower()}.csv")

    # The whole provincial record, not only the configured window: districts
    # report over different spans and restricting first would discard overlap.
    prov = pk[pk["adm_1_name"].astype(str).str.upper() == province].copy()
    prov["district"] = prov["adm_2_name"].map(canonical)

    districts = {}
    rejected = []
    for name, grp in prov.dropna(subset=["adm_2_name"]).groupby("district"):
        s = grp.groupby("start")["dengue_total"].sum().sort_index()
        if s.sum() >= MIN_TOTAL_CASES and len(s) >= MIN_WEEKS:
            districts[name] = s
        elif s.sum() >= MIN_TOTAL_CASES:
            rejected.append(f"{name} ({len(s)} wks)")

    if rejected:
        print(f"  excluded for reporting fewer than {MIN_WEEKS} weeks: "
              f"{', '.join(rejected)}")
    if len(districts) < 2:
        print(f"  only {len(districts)} usable district(s); "
              f"cannot test aggregation here")
        continue

    # Restrict to weeks every district shares, so the comparison is like for
    # like rather than over different periods.
    common = set.intersection(*(set(s.index) for s in districts.values()))
    weeks = pd.Series(sorted(common))
    if len(weeks) < 12:
        print(f"  only {len(weeks)} shared weeks; too few to fit")
        continue
    # Climate must span the shared weeks; the series was downloaded with a
    # year of lead-in, but the districts can report past the configured end.
    if weeks.max() > clim.index.max():
        weeks = weeks[weeks <= clim.index.max()]
        print(f"  trimmed to {len(weeks)} weeks covered by the climate record")
        if len(weeks) < 12:
            continue
    print(f"  {len(districts)} districts, {len(weeks)} shared weeks: "
          f"{weeks.min().date()} to {weeks.max().date()}")

    # --- fit each district separately -------------------------------------
    district_results = {}
    for name, s in sorted(districts.items()):
        cases = s.reindex(weeks).to_numpy()
        pop = DISTRICT_POP_M.get(name, 1.0)
        data = build_dataset(weeks, cases, clim, pop, name)
        try:
            res = fit(data, cfg, fixed, model="constant", observation=OBS)
        except Exception as exc:
            print(f"    {name:18s} fit failed: {str(exc)[:50]}")
            continue
        r0 = basic_reproduction_number(res.theta["beta_0"], fixed)
        district_results[name] = r0
        print(f"    {name:18s} {int(cases.sum()):6d} cases   R0 = {r0:.3f}")
        rows.append(dict(province=province, level="district", unit=name,
                         weeks=len(weeks), cases=int(cases.sum()), R0=round(r0, 4),
                         pop_at_risk_M=round(res.theta["pop_frac"] * pop, 4),
                         starts=f"{res.n_converged}/{res.n_starts}"))

    # --- fit the aggregate of exactly those districts ----------------------
    agg_cases = sum(s.reindex(weeks).to_numpy() for s in districts.values())
    agg_pop = sum(DISTRICT_POP_M.get(n, 1.0) for n in districts)
    agg_data = build_dataset(weeks, agg_cases, clim, agg_pop,
                             f"{province} aggregate")
    res_agg = fit(agg_data, cfg, fixed, model="constant", observation=OBS)
    r0_agg = basic_reproduction_number(res_agg.theta["beta_0"], fixed)
    print(f"    {'AGGREGATE':18s} {int(agg_cases.sum()):6d} cases   "
          f"R0 = {r0_agg:.3f}")
    rows.append(dict(province=province, level="aggregate",
                     unit=f"{province} (sum of {len(districts)})",
                     weeks=len(weeks), cases=int(agg_cases.sum()),
                     R0=round(r0_agg, 4),
                     pop_at_risk_M=round(res_agg.theta["pop_frac"] * agg_pop, 4),
                     starts=f"{res_agg.n_converged}/{res_agg.n_starts}"))

    # --- the comparison ----------------------------------------------------
    vals = np.array(list(district_results.values()))
    if len(vals) >= 2:
        inside = vals.min() <= r0_agg <= vals.max()
        print(f"\n    district R0 range: {vals.min():.3f} to {vals.max():.3f} "
              f"(mean {vals.mean():.3f})")
        print(f"    aggregate R0     : {r0_agg:.3f}  "
              f"{'inside' if inside else 'OUTSIDE'} the district range")
        print(f"    bias vs district mean: {r0_agg - vals.mean():+.3f} "
              f"({(r0_agg - vals.mean()) / vals.mean() * 100:+.1f}%)")
        panels.append((province, weeks, districts, agg_cases,
                       district_results, r0_agg))

# --- output -----------------------------------------------------------------
tab = pd.DataFrame(rows)
tab.to_csv(tables / "07_aggregation_bias.csv", index=False)
print(f"\n{'=' * 74}")
print(tab.to_string(index=False))
print(f"\nTable: {tables / '07_aggregation_bias.csv'}")

if panels:
    fig, axes = plt.subplots(len(panels), 2, figsize=(13, 4.2 * len(panels)),
                             squeeze=False)
    for row, (province, weeks, districts, agg_cases, dres, r0_agg) in enumerate(panels):
        ax = axes[row][0]
        # The aggregate is drawn first and thinner: drawn last and thick it sits
        # exactly on top of whichever district dominates, which made the earlier
        # version of this figure look as though the sum exceeded its parts.
        ax.plot(weeks, agg_cases, color="k", lw=3.0, alpha=0.28,
                label="aggregate", zorder=1)
        for name, s in sorted(districts.items()):
            ax.plot(weeks, s.reindex(weeks).to_numpy(), lw=1.4, zorder=3,
                    label=name.title())
        ax.set_title(f"{province.title()} 2021 — districts and their sum")
        ax.set_ylabel("cases / week")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25)

        ax = axes[row][1]
        names = list(dres)
        ax.barh([n.title() for n in names], [dres[n] for n in names],
                color="#c9ccd1")
        ax.axvline(r0_agg, color="#8c1c13", lw=2,
                   label=f"aggregate $R_0$ = {r0_agg:.2f}")
        ax.axvline(np.mean(list(dres.values())), color="#1f6f8b", ls="--", lw=1.5,
                   label=f"district mean = {np.mean(list(dres.values())):.2f}")
        ax.set_xlabel("$R_0$")
        ax.set_title("Does aggregation change the estimate?")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25, axis="x")
    fig.tight_layout()
    fig.savefig(figures / "07_aggregation_bias.png", dpi=150)
    print(f"Figure: {figures / '07_aggregation_bias.png'}")
