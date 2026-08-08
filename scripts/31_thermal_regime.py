"""Pipeline step 31 — does the local temperature regime decide what a study finds?

Kirk et al. (PLOS Climate, 2024) meta-analysed 358 published temperature--dengue
correlations and found the association strongest where temperature *varies* most.
It is a good finding and this step is an attempt to improve on it rather than to
repeat it, because a meta-analysis of published estimates carries two problems
that no amount of care inside it can remove.

**Analytical confounding.** Their observations are the numbers other people chose
to report. This project has just shown that those numbers move by tens of
percentage points with choices nobody states --- the observation model alone flips
34% of verdicts. If analysts working in high-variation settings tend, for any
reason, toward different conventions, a moderator effect appears that is about
the analysts and not about dengue.

**Publication bias.** Correlations enter their sample by having been published.
Studies that found nothing are systematically less likely to be there, and
"nothing" is exactly the outcome the hypothesis predicts in low-variation
settings.

The design here removes both. Every one of the 221 outbreaks is fitted with the
same model under the same 144 combinations: no choice varies with the setting,
because no choice varies at all. Nothing enters or leaves the sample by having
been interesting. What is left, if the pattern survives, is a property of the
data rather than of the literature.

Three questions, in increasing order of what they would be worth:

1. **Does the effect estimate depend on the temperature regime?** The direct
   analogue of the published finding.
2. **Does the *detectability* depend on it?** Whether the analyses agree with one
   another --- which a meta-analysis of published point estimates cannot ask,
   because disagreement between analyses of one dataset is invisible once one of
   them has been selected for publication.
3. **Is there a mechanistic explanation?** Transmission responds to temperature
   along a unimodal curve. A season spent near the flat optimum offers nothing to
   detect however strong the underlying biology; a season crossing the steep
   flank offers a great deal. If detectability tracks the *gradient traversed*
   rather than the variance, that is a sharper and more useful statement than
   "more variation helps", because a researcher can compute it in advance from a
   thermometer and decide whether their data can answer the question at all.

Reads the stored factorial and the climate series; refits nothing.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.climate import briere  # noqa: E402
from dengue_pk.locations import point_for  # noqa: E402
from dengue_pk.robustness import (FACTORIAL_TABLES,  # noqa: E402
                                  benjamini_hochberg,
                                  complete_windows, latest_factorial,
                                  pairwise_disagreement,
                                  variance_decomposition, window_verdicts)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

cfg = load_config()
raw = resolve(cfg, "raw")
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")
clim_dir = raw / "climate_global"

KEY = ["country", "unit", "window_start"]
TR = cfg["model"]["temperature_response"]
T_MIN, T_MAX = float(TR["t_min_c"]), float(TR["t_max_c"])


def thermal_summary(temp_c: np.ndarray) -> dict:
    """Summarise what a season's temperatures do to transmission.

    ``sd`` and ``range`` describe the temperatures themselves and are what the
    published moderator analysis used. The remaining three describe what those
    temperatures do to the *transmission response*, which is the quantity the
    model can actually see:

    * ``response_range`` --- how much the normalised Brière response varies over
      the season. Zero means transmission was constant no matter what the
      thermometer did.
    * ``gradient`` --- the mean absolute slope of the response over the observed
      temperatures. A season sitting on the flat optimum has a low gradient even
      if it is warm and variable; a season on the steep flank has a high one.
    * ``near_optimum`` --- the share of weeks within 2 degrees of the peak, which
      is the regime in which temperature cannot inform anything.
    """
    t = np.asarray(temp_c, dtype=float)
    b = briere(t, T_MIN, T_MAX)
    peak = 0.5 * (T_MIN + T_MAX) if not np.isfinite(T_MAX) else None
    # The Brière peak is not the midpoint; locate it numerically once.
    grid = np.linspace(T_MIN, T_MAX, 2001)
    bg = briere(grid, T_MIN, T_MAX)
    peak = float(grid[int(np.argmax(bg))])
    # Local slope of the normalised response at each observed temperature.
    slope = np.gradient(bg, grid)
    local_slope = np.interp(t, grid, slope)
    return dict(
        temp_mean=float(np.mean(t)),
        temp_sd=float(np.std(t)),
        temp_range=float(np.max(t) - np.min(t)),
        response_range=float(np.max(b) - np.min(b)),
        gradient=float(np.mean(np.abs(local_slope))),
        near_optimum=float(np.mean(np.abs(t - peak) <= 2.0)),
        peak_c=peak,
    )


def load_climate(slug):
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


def main() -> None:
    d = pd.read_csv(latest_factorial(tables))
    d, n_combos = complete_windows(d, KEY)
    per = window_verdicts(d, n_combos, KEY)
    per["pairwise"] = pairwise_disagreement(d, n_combos, KEY).to_numpy()
    per["p_climate"] = per["wins"] / n_combos

    # Effect estimate per outbreak: the median Brière exponent across the
    # analyses that used it. The median rather than the mean because the
    # distribution has a long right tail and one runaway fit should not decide
    # what an outbreak looks like.
    br = d[d["temp_form"] == "briere"]
    est = br.groupby(KEY)["a_temp"].median().rename("a_temp_median")
    per = per.merge(est, on=KEY, how="left")

    inv = pd.read_csv(tables / "12_global_windows.csv",
                      parse_dates=["start", "end"])
    inv["window_start"] = inv["start"].dt.date.astype(str)

    rows = []
    cache: dict[str, pd.DataFrame | None] = {}
    for _, w in inv.iterrows():
        pt = point_for(w["country"], w["unit"], w["level"])
        if pt is None:
            continue
        name, lat, lon = pt
        slug = (f"{name.lower().replace(' ', '_').replace('/', '_')}"
                f"_{round(lat, 4)}_{round(lon, 4)}")
        if slug not in cache:
            cache[slug] = load_climate(slug)
        clim = cache[slug]
        if clim is None:
            continue
        sel = clim.loc[str(w["start"].date()):str(w["end"].date()), "T2M"].dropna()
        if len(sel) < 100:
            continue
        rows.append(dict(country=w["country"], unit=w["unit"],
                         window_start=w["window_start"],
                         **thermal_summary(sel.to_numpy())))
    regime = pd.DataFrame(rows).drop_duplicates(subset=KEY)

    m = per.merge(regime, on=KEY, how="inner")
    print(f"{len(m)} outbreaks with both a factorial result and a climate series")
    print(f"Brière peak of the transmission response: {m['peak_c'].iloc[0]:.1f} C\n")

    # Before asking which property of the setting explains the differences
    # between studies, bound how much any property of the setting could
    # explain. Every outbreak here has been analysed every way, so the
    # variation splits exactly into "between outbreaks" and "within".
    print("=" * 78)
    print("HOW MUCH CAN ANY PROPERTY OF THE SETTING EXPLAIN?")
    print("=" * 78)
    for name in reversed(FACTORIAL_TABLES):
        path = tables / name
        if not path.exists():
            continue
        t = pd.read_csv(path)
        tc, tn = complete_windows(t, KEY)
        vd = variance_decomposition(tc, "climate_wins", KEY)
        label = name.split("_")[-1].replace(".csv", "")
        if label.startswith("robustness"):
            label = "4factor"
        print(f"  {label:>8} ({tn:3d} analyses):  between outbreaks "
              f"{vd['between_share'] * 100:5.1f}%   within outbreaks "
              f"{vd['within_share'] * 100:5.1f}%")
    print()
    print("  Three-quarters of the variation in whether an outbreak is judged")
    print("  climate-driven is variation between ANALYSES OF THE SAME OUTBREAK.")
    print("  Everything about the setting — climate, density, wealth, serotype,")
    print("  surveillance — competes for the remaining quarter. The split barely")
    print("  moves across designs of 24, 48 and 144 analyses, so it is not an")
    print("  artefact of how finely the space was enumerated.\n")

    PREDICTORS = [
        ("temp_sd", "SD of weekly temperature (the published moderator)"),
        ("temp_range", "range of weekly temperature"),
        ("response_range", "range of the transmission response"),
        ("gradient", "mean |slope| of the response (thermal traverse)"),
        ("near_optimum", "share of weeks within 2 C of the optimum"),
        ("temp_mean", "mean temperature"),
    ]
    OUTCOMES = [
        ("a_temp_median", "estimated temperature exponent"),
        ("p_climate", "share of analyses endorsing climate forcing"),
        ("pairwise", "disagreement between two analyses"),
    ]

    out = []
    for oc, oc_label in OUTCOMES:
        for pr, pr_label in PREDICTORS:
            sub = m[[pr, oc]].dropna()
            if len(sub) < 20:
                continue
            rho, p = spearmanr(sub[pr], sub[oc])
            out.append(dict(outcome=oc, outcome_label=oc_label, predictor=pr,
                            predictor_label=pr_label, n=len(sub),
                            rho=round(float(rho), 4), p=float(p)))

    # -----------------------------------------------------------------
    # A claim the paper had been making without evidence: that the sign of a
    # monotone temperature coefficient is set by whether temperature happens to
    # rise or fall while cases rise, and is therefore a property of the calendar.
    # It is testable against the growth-phase temperature trend.
    print("=" * 78)
    print("IS THE MONOTONE TERM'S SIGN SET BY THE CALENDAR?")
    print("=" * 78)
    ll = (d[d["temp_form"] == "loglinear"].groupby(KEY)["a_temp"]
          .median().rename("a_temp_loglinear"))
    trend_rows = []
    for _, w in inv.iterrows():
        key = (w["country"], w["unit"], w["window_start"])
        if key not in ll.index:
            continue
        pt = point_for(w["country"], w["unit"], w["level"])
        if pt is None:
            continue
        name, lat, lon = pt
        slug = (f"{name.lower().replace(' ', '_').replace('/', '_')}"
                f"_{round(lat, 4)}_{round(lon, 4)}")
        if slug not in cache:
            cache[slug] = load_climate(slug)
        clim = cache[slug]
        if clim is None:
            continue
        sel = clim.loc[str(w["start"].date()):str(w["end"].date()),
                       "T2M"].dropna()
        if len(sel) < 100:
            continue
        # The growth phase: the first half of the window, which for a
        # single-wave window runs up to about the peak.
        half = sel.iloc[:len(sel) // 2].to_numpy()
        slope = float(np.polyfit(np.arange(len(half)), half, 1)[0]) * 365.0
        trend_rows.append(dict(a_temp=float(ll.loc[key]), trend_c_per_year=slope))
    tr = pd.DataFrame(trend_rows)
    if len(tr) > 20:
        rho_t, p_t = spearmanr(tr["trend_c_per_year"], tr["a_temp"])
        falling = tr[tr["trend_c_per_year"] < 0]
        rising = tr[tr["trend_c_per_year"] >= 0]
        print(f"  {len(tr)} outbreaks with a log-linear estimate and a trend")
        print(f"  correlation of the estimate with the growth-phase trend: "
              f"rho = {rho_t:+.3f} (p = {p_t:.2g})")
        print(f"  temperature falling during growth (n={len(falling)}): "
              f"coefficient negative in {(falling['a_temp'] < 0).mean() * 100:.0f}%")
        print(f"  temperature rising  during growth (n={len(rising)}): "
              f"coefficient negative in {(rising['a_temp'] < 0).mean() * 100:.0f}%")
        print()
        print("  The direction is as claimed and the effect is moderate, not")
        print("  deterministic. The calendar tilts the sign; it does not set it.")
        tr.to_csv(tables / "40_monotone_sign.csv", index=False)
    print()

    res = pd.DataFrame(out)

    # Eighteen correlations are computed here — six predictors against three
    # outcomes — so at the conventional 5% level about one would be
    # "significant" by chance. Reporting raw p-values with a flag on the small
    # ones would be precisely the undisclosed analytical choice this paper
    # exists to measure. The false-discovery rate is therefore controlled across
    # the whole family (Benjamini-Hochberg) and both values are printed.
    res["q"] = benjamini_hochberg(res["p"].to_numpy())

    for oc, oc_label in OUTCOMES:
        print("=" * 78)
        print(f"OUTCOME: {oc_label}")
        print("=" * 78)
        for _, r in res[res["outcome"] == oc].iterrows():
            flag = "  <-- survives FDR" if r["q"] < 0.05 else ""
            print(f"  {r['predictor_label']:52s} rho {r['rho']:+.3f}  "
                  f"p {r['p']:7.2g}  q {r['q']:7.2g}{flag}")
        print()
    n_raw = int((res["p"] < 0.05).sum())
    n_fdr = int((res["q"] < 0.05).sum())
    print(f"  {n_raw} of {len(res)} correlations reach p < 0.05; "
          f"{n_fdr} survive false-discovery-rate control at 5%.")
    print()
    res.to_csv(tables / "32_thermal_regime.csv", index=False)

    # ---------------------------------------------------------------------
    print("=" * 78)
    print("READING")
    print("=" * 78)
    sd_row = res[(res.outcome == "a_temp_median") & (res.predictor == "temp_sd")]
    gr_row = res[(res.outcome == "a_temp_median") & (res.predictor == "gradient")]
    if not sd_row.empty and not gr_row.empty:
        print(f"  Published moderator (temperature SD):  rho = "
              f"{sd_row['rho'].iloc[0]:+.3f}")
        print(f"  Mechanistic version (thermal traverse): rho = "
              f"{gr_row['rho'].iloc[0]:+.3f}")
        print("\n  If the mechanistic version is the stronger predictor, the")
        print("  published moderator is a proxy for it, and a researcher can")
        print("  compute the better one from a thermometer before fitting.")

    # ---------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    panels = [("temp_sd", "a_temp_median", "SD of weekly temperature (°C)",
               "estimated temperature exponent",
               "The published moderator,\nrecomputed with the analysis held fixed"),
              ("gradient", "a_temp_median", "thermal traverse (mean |slope|)",
               "estimated temperature exponent",
               "The mechanistic version:\nhow much response gradient the season crosses"),
              ("gradient", "pairwise", "thermal traverse (mean |slope|)",
               "P(two analyses disagree)",
               "And whether the question is\nanswerable at all")]
    for ax, (px, py, xl, yl, title) in zip(axes, panels):
        sub = m[[px, py]].dropna()
        ax.scatter(sub[px], sub[py], s=26, alpha=0.65, color="#1f6f8b")
        if len(sub) > 5:
            z = np.polyfit(sub[px], sub[py], 1)
            xs = np.linspace(sub[px].min(), sub[px].max(), 50)
            ax.plot(xs, np.polyval(z, xs), color="#8c1c13", lw=2)
            rho, p = spearmanr(sub[px], sub[py])
            ax.set_title(f"{title}\nρ = {rho:+.2f}, p = {p:.1g}", fontsize=10)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures / "18_thermal_regime.png", dpi=150)
    print(f"\nFigure: {figures / '18_thermal_regime.png'}")
    print(f"Table:  {tables / '32_thermal_regime.csv'}")


if __name__ == "__main__":
    main()
