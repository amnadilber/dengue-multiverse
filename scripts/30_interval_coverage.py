"""Pipeline step 30 — do the intervals published with these models cover the truth?

Everything before this measures disagreement. Disagreement is a diagnosis, and a
diagnosis is worth less than a working instrument. This step asks the
constructive question instead: **when a paper reports a transmission parameter
with a 95% confidence interval, how often does that interval actually contain the
right answer — and can we build one that does?**

The question is answerable here and almost nowhere else, because the truth is
known. Counts are simulated from each window's own fitted parameters with the
temperature exponent set to a known value, and every analysis in the factorial is
then asked to recover it. A conventional analysis reports one estimate and one
interval from one set of choices. We can see whether that interval covers.

Two intervals are compared.

**Conventional.** One analysis, chosen as an analyst would choose it, reporting
``theta_hat +/- 1.96 * SE`` from the asymptotic standard error. This is what
appears in papers. Its coverage is computed over every cell of the factorial, so
the figure is not a statement about one arbitrary choice but about what an
analyst drawing any defensible set of choices would obtain.

**Multiverse.** All analyses of the same outbreak combined by Rubin's rules,
which were built for multiple imputation and apply unchanged here because the
structure is identical --- several analyses of one dataset, each with its own
estimate and its own within-analysis variance:

    theta_bar = mean of the estimates
    W         = mean of the within-analysis variances        (sampling uncertainty)
    B         = variance between the estimates               (analytical uncertainty)
    T         = W + (1 + 1/m) B                              (total)

The interval is ``theta_bar +/- t * sqrt(T)`` on Rubin's degrees of freedom. The
one line that matters: a conventional interval reports ``sqrt(W)`` and behaves as
though ``B`` were zero.

Coverage is the metric because it is the property a confidence interval claims to
have. An interval that is narrow and wrong is worse than a wide one that is
honest, and only coverage distinguishes them.

Restricted to fits using the Brière temperature parameterisation, since that is
the form the data were generated under and the only one in which the target
parameter means the same thing. That is 72 of the 144 combinations. Analyses
using a different model structure or a different fixed parameter set are kept:
their estimates are biased by that mis-specification, and coverage should be
charged for that bias, because a real analyst does not know they have chosen the
wrong structure either.
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
from scipy import stats  # noqa: E402

OBSERVATIONS = ("nb", "poisson")
RAIN_LAGS = (3, 5, 7)
TRAIN_FRACS = (1.0, 0.75)
STRUCTURES = ("hostvector", "seir", "hostvector_tempmort")
PARAM_SETS = ("central", "alt")

#: The generating truth. a_temp is the exponent on the Brière temperature
#: response: 0 means temperature does not matter, 1 means the literature
#: response applies in full.
TRUE_A_TEMP = 1.00
TRUE_A_RAIN = 0.30

BASE = load_config()
raw = resolve(BASE, "raw")
tables = resolve(BASE, "tables")
clim_dir = raw / "climate_global"
OUT = tables / "30_interval_coverage.csv"

FIXED = {"central": FixedParams.from_config(BASE, "fixed"),
         "alt": FixedParams.from_config(BASE, "fixed_alt")}

N_COMBOS = (len(OBSERVATIONS) * len(RAIN_LAGS) * len(TRAIN_FRACS)
            * len(STRUCTURES) * len(PARAM_SETS))


def cfg_for(observation, lag, structure):
    c = copy.deepcopy(BASE)
    c["model"]["temperature_form"] = "briere"
    c["model"]["structure"] = structure
    c["model"]["climate_forcing"]["rain_lag_weeks"] = lag
    c["inference"]["observation"] = observation
    c["inference"]["classical"]["n_restarts"] = 3
    return c


def fit_window(package):
    """Every Brière combination on one window of simulated data."""
    meta, data_by_lag, seed = package
    label = f"{meta['country']}/{meta['unit']} {meta['window_start']}"
    if not data_by_lag:
        return label, [], "no usable lag"
    rng = np.random.default_rng(seed)
    ref = data_by_lag[min(data_by_lag)]

    set_temperature_form("briere")
    base_cfg = cfg_for("nb", min(data_by_lag), "hostvector")
    try:
        anchor = fit(ref, base_cfg, FIXED["central"], model="constant",
                     observation="nb")
    except Exception as exc:                        # noqa: BLE001
        return label, [], f"anchor failed: {exc}"
    k_hat = float(anchor.nb_k or 20.0)

    truth = dict(anchor.theta)
    truth.update(a_temp=TRUE_A_TEMP, a_rain=TRUE_A_RAIN)
    # Match mean transmission to the anchor, as in step 29, so the generated
    # epidemic is the size the real one was.
    mult = make_forcing({**truth, "beta_0": 1.0}, ref, "climate",
                        base_cfg).on_grid(ref.days)
    truth["beta_0"] = float(anchor.theta["beta_0"]
                            / max(float(np.mean(mult)), 1e-12))
    try:
        mu_true = predict(truth, ref, FIXED["central"], "climate", base_cfg)
    except IntegrationFailure:
        return label, [], "generating integration diverged"
    lam = rng.gamma(shape=k_hat, scale=np.maximum(mu_true, 1e-12) / k_hat)
    synth = rng.poisson(lam).astype(float)
    if synth.sum() < 100:
        return label, [], "simulated epidemic too small"

    rows = []
    for observation, lag, frac, structure, pset in itertools.product(
            OBSERVATIONS, RAIN_LAGS, TRAIN_FRACS, STRUCTURES, PARAM_SETS):
        if lag not in data_by_lag:
            continue
        d0 = data_by_lag[lag]
        full = Dataset(d0.days, synth, d0.temp_c, d0.z_rain, d0.population,
                       d0.label)
        train = full if frac >= 1.0 else full.head(
            int(round(len(full.days) * frac)))
        set_temperature_form("briere")
        cfg = cfg_for(observation, lag, structure)
        try:
            res = fit(train, cfg, FIXED[pset], model="climate",
                      observation=observation, start_from=truth,
                      n_starts_override=3,
                      fixed_nb_k=k_hat if observation == "nb" else None)
        except (RuntimeError, IntegrationFailure, ValueError):
            continue

        row = dict(country=meta["country"], unit=meta["unit"],
                   window_start=meta["window_start"],
                   observation=observation, rain_lag=lag, train_frac=frac,
                   structure=structure, params=pset,
                   true_a_temp=TRUE_A_TEMP, true_a_rain=TRUE_A_RAIN,
                   # The generating transmission coefficient, so that coverage
                   # can be computed for R0 — the quantity this field actually
                   # reports. R0 is a fixed multiple of beta_0 within a
                   # structure, so its coverage equals beta_0's.
                   true_beta_0=round(float(truth["beta_0"]), 6))
        for name in ("a_temp", "a_rain", "beta_0"):
            est = float(res.theta[name])
            se = float(res.stderr.get(name, np.nan))
            row[f"{name}_hat"] = round(est, 6)
            row[f"{name}_se"] = round(se, 6) if np.isfinite(se) else np.nan
        rows.append(row)
    return label, rows, None


def build_packages(n_windows: int):
    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])
    inv = inv.iloc[:: max(len(inv) // n_windows, 1)].head(n_windows)
    print(f"{len(inv)} windows x {N_COMBOS} Brière combinations "
          f"= {len(inv) * N_COMBOS:,} fits\n", flush=True)
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
        if by_lag:
            packages.append(({"country": w["country"], "unit": w["unit"],
                              "window_start": str(w["start"].date())},
                             by_lag, BASE["seed"] + 7919 * idx))
    return packages


# ---------------------------------------------------------------------------
def rubin(estimates: np.ndarray, ses: np.ndarray, level: float = 0.95):
    """Combine analyses of one dataset by Rubin's rules.

    Returns (point, lower, upper, within_sd, between_sd). ``m`` is the number of
    analyses; degrees of freedom follow Rubin (1987).
    """
    ok = np.isfinite(estimates) & np.isfinite(ses)
    est, se = estimates[ok], ses[ok]
    m = len(est)
    if m == 0:
        return (np.nan,) * 5
    if m == 1:
        z = stats.norm.ppf(0.5 + level / 2)
        return est[0], est[0] - z * se[0], est[0] + z * se[0], se[0], 0.0
    q_bar = float(np.mean(est))
    w = float(np.mean(se ** 2))                    # within-analysis variance
    b = float(np.var(est, ddof=1))                 # between-analysis variance
    t = w + (1.0 + 1.0 / m) * b                    # total variance
    if t <= 0 or not np.isfinite(t):
        return (np.nan,) * 5
    # Rubin's degrees of freedom; guarded because b = 0 sends it to infinity.
    lam = (1.0 + 1.0 / m) * b / t
    df = (m - 1) / max(lam ** 2, 1e-12) if lam > 0 else np.inf
    crit = stats.t.ppf(0.5 + level / 2, max(df, 1.0)) if np.isfinite(df) \
        else stats.norm.ppf(0.5 + level / 2)
    half = crit * np.sqrt(t)
    return q_bar, q_bar - half, q_bar + half, np.sqrt(w), np.sqrt(b)


def report(d: pd.DataFrame) -> None:
    KEY = ["country", "unit", "window_start"]
    z = stats.norm.ppf(0.975)

    print("=" * 78)
    print("1. DOES A CONVENTIONAL 95% INTERVAL COVER THE TRUTH?")
    print("=" * 78)
    for name, truth_col in (("a_temp", "true_a_temp"), ("a_rain", "true_a_rain")):
        est, se = d[f"{name}_hat"], d[f"{name}_se"]
        ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
        lo, hi = est[ok] - z * se[ok], est[ok] + z * se[ok]
        covers = (lo <= d.loc[ok, truth_col]) & (d.loc[ok, truth_col] <= hi)
        print(f"  {name:8s} nominal 95%, actual "
              f"\033[1m{covers.mean() * 100:5.1f}%\033[0m "
              f"over {ok.sum():,} fits "
              f"({(~np.isfinite(se)).sum():,} had no usable standard error)")
        print(f"           median half-width {(z * se[ok]).median():.4f}, "
              f"median |bias| {np.abs(est[ok] - d.loc[ok, truth_col]).median():.4f}")

    print("\n" + "=" * 78)
    print("2. DOES COMBINING THE ANALYSES FIX IT?")
    print("=" * 78)
    rows = []
    for key, g in d.groupby(KEY):
        for name, truth_col in (("a_temp", "true_a_temp"),
                                ("a_rain", "true_a_rain")):
            q, lo, hi, w_sd, b_sd = rubin(g[f"{name}_hat"].to_numpy(),
                                          g[f"{name}_se"].to_numpy())
            if not np.isfinite(q):
                continue
            truth = float(g[truth_col].iloc[0])
            rows.append(dict(country=key[0], unit=key[1], window_start=key[2],
                             parameter=name, m=len(g), point=q, lower=lo,
                             upper=hi, within_sd=w_sd, between_sd=b_sd,
                             truth=truth, covers=bool(lo <= truth <= hi),
                             half_width=(hi - lo) / 2))
    mv = pd.DataFrame(rows)
    for name in ("a_temp", "a_rain"):
        sub = mv[mv["parameter"] == name]
        if sub.empty:
            continue
        print(f"  {name:8s} nominal 95%, actual "
              f"\033[1m{sub['covers'].mean() * 100:5.1f}%\033[0m "
              f"over {len(sub)} outbreaks")
        print(f"           median half-width {sub['half_width'].median():.4f}; "
              f"analytical SD is {(sub['between_sd'] / sub['within_sd'].replace(0, np.nan)).median():.1f}x "
              f"the sampling SD")
    mv.to_csv(tables / "31_multiverse_intervals.csv", index=False)

    print("\n" + "=" * 78)
    print("3. HOW MUCH TOO NARROW IS THE CONVENTIONAL INTERVAL?")
    print("=" * 78)
    for name in ("a_temp", "a_rain"):
        sub = mv[mv["parameter"] == name]
        conv = d.groupby(KEY)[f"{name}_se"].median() * z
        ratio = (sub.set_index(KEY)["half_width"] / conv).replace(
            [np.inf, -np.inf], np.nan).dropna()
        if ratio.empty:
            continue
        print(f"  {name:8s} the honest interval is "
              f"\033[1m{ratio.median():.1f}x\033[0m wider than the one a single "
              f"analysis reports")
        print(f"           (quartiles {ratio.quantile(.25):.1f}x to "
              f"{ratio.quantile(.75):.1f}x)")

    print("\n" + "=" * 78)
    print("4. THE PARAMETER THE FIELD ACTUALLY REPORTS")
    print("=" * 78)
    print("  R0 is a fixed multiple of the transmission coefficient within a")
    print("  structure, so the two have identical coverage. This is the number")
    print("  a dengue paper puts in its abstract.\n")
    if "true_beta_0" in d.columns:
        est, se, tr = d["beta_0_hat"], d["beta_0_se"], d["true_beta_0"]
        ok = np.isfinite(est) & np.isfinite(se) & (se > 0) & np.isfinite(tr)
        conv = float(((est[ok] - z * se[ok] <= tr[ok])
                      & (tr[ok] <= est[ok] + z * se[ok])).mean())
        cov = []
        for _, g in d.groupby(KEY):
            q, lo_i, hi_i, _, _ = rubin(g["beta_0_hat"].to_numpy(),
                                        g["beta_0_se"].to_numpy())
            if np.isfinite(q):
                cov.append(lo_i <= float(g["true_beta_0"].iloc[0]) <= hi_i)
        print(f"  conventional 95% interval: {conv * 100:5.1f}%   "
              f"({int(ok.sum()):,} fits)")
        print(f"  multiverse (Rubin):        {np.mean(cov) * 100:5.1f}%   "
              f"({len(cov)} outbreaks)")
        ratio = (est[ok] / tr[ok]).replace([np.inf, -np.inf], np.nan).dropna()
        print(f"\n  estimate divided by truth: median {ratio.median():.2f}, "
              f"IQR {ratio.quantile(.25):.2f} to {ratio.quantile(.75):.2f}")
        print("  In the middle half of analyses the reported R0 is between those")
        print("  multiples of the value that generated the data.")
    else:
        print("  (rerun this step: true_beta_0 is not in the stored table)")

    print("\n" + "=" * 78)
    print("5. HOW MANY ANALYSES DOES THE INTERVAL NEED?")
    print("=" * 78)
    print("  A recommendation costing the whole factorial would be admired and")
    print("  not used. Drawing m analyses at random from each outbreak's set and")
    print("  combining only those, averaged over 20 draws.\n")
    rng2 = np.random.default_rng(BASE["seed"])
    groups = [g for _, g in d.groupby(KEY)]
    m_rows = []
    targets = [("a_rain", "true_a_rain"), ("a_temp", "true_a_temp")]
    if "true_beta_0" in d.columns:
        targets.append(("beta_0", "true_beta_0"))
    header = "  " + f"{'m':>4}" + "".join(f"{n:>14}" for n, _ in targets) \
        + f"{'width (a_rain)':>17}"
    print(header)
    for m in (2, 4, 8, 16, 24, 48, 72):
        line, row = f"  {m:4d}", dict(m=m)
        width = np.nan
        for name, truth_col in targets:
            cov, wid = [], []
            for _ in range(20):
                for g in groups:
                    if len(g) < m:
                        continue
                    sub = g.iloc[rng2.choice(len(g), m, replace=False)]
                    q, lo_i, hi_i, _, _ = rubin(sub[f"{name}_hat"].to_numpy(),
                                                sub[f"{name}_se"].to_numpy())
                    if np.isfinite(q):
                        cov.append(lo_i <= float(sub[truth_col].iloc[0]) <= hi_i)
                        wid.append((hi_i - lo_i) / 2)
            if not cov:
                continue
            line += f"{np.mean(cov) * 100:13.1f}%"
            row[name] = round(float(np.mean(cov)) * 100, 1)
            if name == "a_rain":
                width = float(np.median(wid))
                row["half_width_a_rain"] = round(width, 3)
        print(line + f"{width:17.3f}")
        m_rows.append(row)
    pd.DataFrame(m_rows).to_csv(tables / "37_analyses_needed.csv", index=False)
    print("\n  Table: " + str(tables / "37_analyses_needed.csv"))

    print("\n" + "=" * 78)
    print("6. DOES IT MATTER WHICH EIGHT?")
    print("=" * 78)
    print("  Eight analyses restore coverage, but a practitioner has to pick")
    print("  eight. Comparing a spanning draw against draws that hold one factor")
    print("  fixed shows which choices must be varied and which need not be.\n")
    groups2 = [g for _, g in d.groupby(KEY)]
    rng3 = np.random.default_rng(BASE["seed"] + 1)

    def draw_random(g, m):
        return g.iloc[rng3.choice(len(g), min(m, len(g)), replace=False)]

    def draw_fixing(col):
        def f(g, m):
            v = rng3.choice(g[col].unique())
            sub = g[g[col] == v]
            if len(sub) < 2:
                return None
            return sub.iloc[rng3.choice(len(sub), min(m, len(sub)),
                                        replace=False)]
        return f

    def coverage_of(selector, name, truth_col, reps=25, m=8):
        got = []
        for _ in range(reps):
            for g in groups2:
                sub = selector(g, m)
                if sub is None or len(sub) < 2:
                    continue
                q, lo_i, hi_i, _, _ = rubin(sub[f"{name}_hat"].to_numpy(),
                                            sub[f"{name}_se"].to_numpy())
                if np.isfinite(q):
                    got.append(lo_i <= float(sub[truth_col].iloc[0]) <= hi_i)
        return float(np.mean(got)) * 100 if got else float("nan")

    which_targets = [("a_rain", "true_a_rain"), ("a_temp", "true_a_temp")]
    if "true_beta_0" in d.columns:
        which_targets.append(("beta_0", "true_beta_0"))
    print("  " + f"{'eight analyses sharing...':32s}"
          + "".join(f"{n:>10}" for n, _ in which_targets))
    which_rows = []
    options = [("nothing (a spanning draw)", draw_random)]
    options += [(f"one {c}", draw_fixing(c))
                for c in ("observation", "structure", "params", "train_frac",
                          "rain_lag") if c in d.columns]
    for label, sel in options:
        vals = [coverage_of(sel, n, t) for n, t in which_targets]
        flag = "   <-- coverage lost" if np.nanmin(vals) < 93 else ""
        print(f"  {label:32s}" + "".join(f"{v:9.1f}%" for v in vals) + flag)
        row = dict(held_fixed=label)
        row.update({n: round(v, 1) for (n, _), v in zip(which_targets, vals)})
        which_rows.append(row)
    pd.DataFrame(which_rows).to_csv(tables / "39_which_analyses.csv", index=False)
    print()
    print("  The choices that define a covariate are the ones that must be")
    print("  varied to cover its coefficient. Holding the rainfall lag fixed")
    print("  loses the rainfall coefficient; holding the observation model or")
    print("  the structure fixed costs nothing. This is not what we expected —")
    print("  the observation model is the largest lever on the *verdict* and")
    print("  almost irrelevant to the *interval*.")

    print("\n" + "=" * 78)
    print("7. WHICH ANALYSES ARE WORST?")
    print("=" * 78)
    for factor in ("observation", "structure", "params", "train_frac"):
        est, se = d["a_temp_hat"], d["a_temp_se"]
        ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
        sub = d[ok].copy()
        sub["covers"] = ((est[ok] - z * se[ok] <= sub["true_a_temp"])
                         & (sub["true_a_temp"] <= est[ok] + z * se[ok]))
        cov = sub.groupby(factor)["covers"].mean() * 100
        print(f"  by {factor:12s} " +
              "   ".join(f"{k}: {v:.1f}%" for k, v in cov.items()))


def main():
    n_windows = int(sys.argv[1]) if len(sys.argv) > 1 else 80
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
            print(f"[{i:3d}/{len(packages)}] {label:40s} {len(got):3d} fits  "
                  f"eta {rate * (len(packages) - i) / 60:5.1f} min", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    print(f"\n{len(d):,} fits in {(time.time() - t0) / 60:.1f} min\n")
    report(d)
    print(f"\nTables: {OUT}, {tables / '31_multiverse_intervals.csv'}")


if __name__ == "__main__":
    main()
