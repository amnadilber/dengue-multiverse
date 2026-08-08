"""Pipeline step 32 — does the decomposition depend on which factors we chose?

The headline of Section \\ref{sec:decomposition} is that roughly three-quarters of
the variation in whether an outbreak is judged climate-driven sits between
analyses of the same outbreak rather than between outbreaks. The obvious
objection, and the one a referee will raise first, is that the number is a
property of the factor list: enumerate twenty choices and the analytical share
rises; enumerate one and it collapses.

The defence offered so far — that the split reads 78.3%, 79.7% and 76.5% on
designs of 24, 48 and 144 cells — is weak, and it is worth saying why rather than
leaving it to be found. Those three designs are *nested*. All three contain the
observation model, which is the dominant factor. Showing that the answer survives
adding cells to the same factor set is not showing that it survives choosing a
different factor set.

This step does the harder version.

**Leave one out.** Drop each factor in turn, recompute on the surviving
combinations, and report the range. If the result rests on one factor, removing
that factor will show it.

**Every subset.** Compute the decomposition for all 63 non-empty subsets of the
six factors and report the distribution. This bounds the claim honestly: the
lowest subset value is what a sceptic could obtain by choosing the factor list
adversarially.

**Interval.** Bootstrap over outbreaks, which are the independent units, so the
headline carries an uncertainty like any other estimate. It had none.

**Excluded windows.** Sixteen of 237 windows were dropped, four because no fit
converged and twelve because only part of the design completed. If those differ
systematically from the rest, the sample is not what it claims.

Reads the stored factorial; refits nothing.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (complete_windows, latest_factorial,  # noqa: E402
                                  variance_decomposition)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")
KEY = ["country", "unit", "window_start"]

FACTORS = ["observation", "temp_form", "rain_lag", "train_frac", "structure",
           "params"]


def decompose_on_subset(d: pd.DataFrame, keep: list[str]) -> dict | None:
    """Decomposition using only the combinations spanned by ``keep``.

    Factors not in ``keep`` are held at their first level rather than averaged
    over, because averaging over a dropped factor would leave its variation in
    the within-outbreak term and defeat the purpose.
    """
    sub = d
    for f in FACTORS:
        if f not in keep:
            first = sorted(d[f].unique(), key=str)[0]
            sub = sub[sub[f] == first]
    if sub.empty:
        return None
    counts = sub.groupby(KEY).size()
    n = int(counts.max())
    if n < 2:
        return None
    sub = sub.set_index(KEY).loc[counts[counts == n].index].reset_index()
    r = variance_decomposition(sub, "climate_wins", KEY)
    r["n_combos"] = n
    r["n_windows"] = sub.groupby(KEY).ngroups
    return r


def main() -> None:
    raw = pd.read_csv(latest_factorial(tables))
    d, n_combos = complete_windows(raw, KEY)
    present = [f for f in FACTORS if f in d.columns]
    base = variance_decomposition(d, "climate_wins", KEY)
    print(f"{len(d):,} fits, {d.groupby(KEY).ngroups} outbreaks, "
          f"{n_combos} combinations\n")

    print("=" * 78)
    print("1. THE HEADLINE, WITH AN INTERVAL")
    print("=" * 78)
    rng = np.random.default_rng(cfg["seed"])
    groups = [g for _, g in d.groupby(KEY)]
    boot = []
    for _ in range(500):
        pick = rng.integers(0, len(groups), len(groups))
        # Relabel the resampled outbreaks so repeats are distinct units.
        frames = []
        for j, i in enumerate(pick):
            g = groups[i].copy()
            g["country"] = f"boot{j}"
            frames.append(g)
        boot.append(variance_decomposition(pd.concat(frames, ignore_index=True),
                                           "climate_wins", KEY)["within_share"])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f"  within-outbreak share: {base['within_share'] * 100:.1f}% "
          f"[95% CI {lo * 100:.1f}–{hi * 100:.1f}]")

    print("\n" + "=" * 78)
    print("2. LEAVE ONE FACTOR OUT")
    print("=" * 78)
    rows = [dict(dropped="none (all six)", within=base["within_share"],
                 n_combos=n_combos)]
    for f in present:
        keep = [x for x in present if x != f]
        r = decompose_on_subset(d, keep)
        if r is None:
            continue
        rows.append(dict(dropped=f, within=r["within_share"],
                         n_combos=r["n_combos"]))
        print(f"  without {f:12s} {r['within_share'] * 100:5.1f}%   "
              f"({r['n_combos']:3d} combinations)")
    loo = pd.DataFrame(rows)
    span = loo[loo["dropped"] != "none (all six)"]["within"]
    print(f"\n  range across the six leave-one-out designs: "
          f"{span.min() * 100:.1f}% to {span.max() * 100:.1f}%")

    print("\n" + "=" * 78)
    print("3. EVERY SUBSET OF THE SIX FACTORS")
    print("=" * 78)
    subs = []
    for k in range(1, len(present) + 1):
        for keep in itertools.combinations(present, k):
            r = decompose_on_subset(d, list(keep))
            if r is None:
                continue
            subs.append(dict(k=k, factors="+".join(keep),
                             within=r["within_share"],
                             n_combos=r["n_combos"]))
    sd = pd.DataFrame(subs)
    print(f"  {len(sd)} non-empty subsets")
    print(f"  within-outbreak share: min {sd['within'].min() * 100:.1f}%, "
          f"median {sd['within'].median() * 100:.1f}%, "
          f"max {sd['within'].max() * 100:.1f}%")
    print("\n  by number of factors included:")
    for k, g in sd.groupby("k"):
        print(f"    {k} factor(s):  median {g['within'].median() * 100:5.1f}%   "
              f"range {g['within'].min() * 100:5.1f}–{g['within'].max() * 100:5.1f}%")
    worst = sd.loc[sd["within"].idxmin()]
    print(f"\n  The most favourable single choice a sceptic could make is to keep")
    print(f"  only [{worst['factors']}], which still leaves "
          f"{worst['within'] * 100:.1f}% within outbreaks.")
    sd.to_csv(tables / "33_decomposition_subsets.csv", index=False)

    print("\n" + "=" * 78)
    print("4. THE SAME SPLIT, APPLIED TO THE NUMBER THE FIELD REPORTS")
    print("=" * 78)
    print("  The verdict is a yes or a no. R0 is what these papers publish, and")
    print("  the decomposition applies to it unchanged. Four treatments of the")
    print("  long right tail, because the contrast and not the figure is the")
    print("  finding.\n")
    print(f"  {'treatment':34s} {'model':9s} {'between':>8} {'within':>8}")
    r0_rows = []
    for label, restrict in (("raw", None),
                            ("log scale", "log"),
                            ("restricted to R0 in [0.5, 20]", (0.5, 20.0)),
                            ("restricted to R0 <= 10", (0.0, 10.0))):
        for col in ("R0_constant", "R0_climate"):
            sub = d.assign(v=d[col].replace([np.inf, -np.inf], np.nan)).dropna(
                subset=["v"])
            if restrict == "log":
                sub = sub.assign(v=np.log(sub["v"].clip(lower=1e-6)))
            elif restrict is not None:
                lo_r, hi_r = restrict
                sub = sub[(sub["v"] >= lo_r) & (sub["v"] <= hi_r)]
            r = variance_decomposition(sub, "v", KEY)
            r0_rows.append(dict(treatment=label, model=col,
                                between=round(r["between_share"] * 100, 1),
                                within=round(r["within_share"] * 100, 1)))
            print(f"  {label:34s} {col.replace('R0_', ''):9s} "
                  f"{r['between_share'] * 100:7.1f}% {r['within_share'] * 100:7.1f}%")
    pd.DataFrame(r0_rows).to_csv(tables / "35_r0_decomposition.csv", index=False)
    spread = d.groupby(KEY)["R0_climate"].agg(lambda s: s.max() - s.min())
    print(f"\n  within-outbreak R0 range: median {spread.median():.2f}"
          f"   (median R0 across all fits {d['R0_climate'].median():.2f})")

    print("\n" + "=" * 78)
    print("5. IS THIS A FACT ABOUT ONE REGION?")
    print("=" * 78)
    per_ct = d.groupby(KEY).size().reset_index().groupby("country").size()
    print(f"  {d.groupby(KEY).ngroups} outbreaks, "
          f"{d.groupby(['country', 'unit']).ngroups} reporting units, "
          f"{d['country'].nunique()} countries")
    print(f"  three best represented supply "
          f"{per_ct.nlargest(3).sum() / per_ct.sum() * 100:.0f}% of the sample: "
          f"{', '.join(per_ct.nlargest(3).index)}\n")
    loo_rows = []
    for country in per_ct.nlargest(8).index:
        sub = d[d["country"] != country]
        r = variance_decomposition(sub, "climate_wins", KEY)
        loo_rows.append(dict(dropped=country, within=r["within_share"],
                             n_windows=sub.groupby(KEY).ngroups))
        print(f"  without {country:28s} {r['within_share'] * 100:5.1f}%   "
              f"(n={sub.groupby(KEY).ngroups})")
    big3 = list(per_ct.nlargest(3).index)
    sub = d[~d["country"].isin(big3)]
    r = variance_decomposition(sub, "climate_wins", KEY)
    print(f"\n  dropping all three at once: {r['within_share'] * 100:.1f}% "
          f"(n={sub.groupby(KEY).ngroups})")
    for lvl, g in d.groupby("level"):
        r = variance_decomposition(g, "climate_wins", KEY)
        print(f"  {lvl:12s} only: {r['within_share'] * 100:.1f}% "
              f"(n={g.groupby(KEY).ngroups})")
    pd.DataFrame(loo_rows).to_csv(tables / "36_leave_one_country_out.csv",
                                  index=False)

    print("\n" + "=" * 78)
    print("6. AND WITH THE OUTBREAKS TREATED AS CLUSTERED WITHIN COUNTRIES")
    print("=" * 78)
    cl = []
    country_groups = {k: g for k, g in d.groupby("country")}
    ckeys = list(country_groups)
    for _ in range(400):
        pick = rng.integers(0, len(ckeys), len(ckeys))
        frames = []
        for j, i in enumerate(pick):
            g = country_groups[ckeys[i]].copy()
            g["country"] = f"boot{j}"
            frames.append(g)
        cl.append(variance_decomposition(pd.concat(frames, ignore_index=True),
                                         "climate_wins", KEY)["within_share"])
    clo, chi = np.percentile(cl, [2.5, 97.5])
    print(f"  resampling outbreaks: [{lo * 100:.1f}, {hi * 100:.1f}]")
    print(f"  resampling countries: [{clo * 100:.1f}, {chi * 100:.1f}]  "
          f"<- the conservative interval")

    print("\n" + "=" * 78)
    print("7. ARE THE EXCLUDED WINDOWS DIFFERENT?")
    print("=" * 78)
    counts = raw.groupby(KEY).size()
    kept = set(counts[counts == n_combos].index)
    inv = pd.read_csv(tables / "12_global_windows.csv", parse_dates=["start"])
    inv["window_start"] = inv["start"].dt.date.astype(str)
    inv["kept"] = [tuple(x) in kept for x in
                   inv[["country", "unit", "start"]].assign(
                       start=inv["window_start"]).itertuples(index=False,
                                                             name=None)]
    for col in ("weeks", "cases", "peak"):
        if col not in inv.columns:
            continue
        a = inv.loc[inv["kept"], col]
        b = inv.loc[~inv["kept"], col]
        if len(b) == 0:
            continue
        print(f"  {col:8s} kept median {a.median():10.1f} (n={len(a)})   "
              f"excluded median {b.median():10.1f} (n={len(b)})")
    print("\n  Excluded windows are reported so a reader can judge whether the")
    print("  sample is what it claims rather than being told that it is.")

    # ---------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    axes[0].barh(loo["dropped"], loo["within"] * 100, color="#1f6f8b")
    axes[0].axvline(base["within_share"] * 100, color="#8c1c13", lw=2,
                    label=f"all six: {base['within_share'] * 100:.1f}%")
    axes[0].set_xlabel("% of variation within outbreaks")
    axes[0].set_title("Dropping any one factor")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3, axis="x")

    axes[1].hist(sd["within"] * 100, bins=20, color="#3f7d58")
    axes[1].axvline(base["within_share"] * 100, color="#8c1c13", lw=2,
                    label="all six factors")
    axes[1].set_xlabel("% of variation within outbreaks")
    axes[1].set_ylabel("subsets of the factor list")
    axes[1].set_title(f"All {len(sd)} subsets of the six factors")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(figures / "19_decomposition_robustness.png", dpi=150)
    print(f"\nFigure: {figures / '19_decomposition_robustness.png'}")
    print(f"Table:  {tables / '33_decomposition_subsets.csv'}")


if __name__ == "__main__":
    main()
