"""Pipeline step 23 — every number the paper quotes, printed in one place.

The paper had drifted from the results twice: once when the factorial grew from
four factors to five, and once when a hard-coded bound made step 18 report 11.5%
instability where the true figure was 88.0%. Both were caught by reading, which
is not a method. This script prints the quoted numbers directly from the result
tables so that a mismatch between the paper and the pipeline is visible in one
place rather than distributed across nine sections.

It computes nothing new. Anything here that disagrees with paper/paper.tex means
the paper is stale.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (FACTORIAL_TABLES,  # noqa: E402
                                  complete_windows, latest_factorial,
                                  pairwise_disagreement, window_verdicts)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")

# The richest factorial available; see dengue_pk.robustness.FACTORIAL_TABLES.
SRC = latest_factorial(tables)
raw = pd.read_csv(SRC)
KEY = ["country", "unit", "window_start"]
d, N_COMBOS = complete_windows(raw, KEY)
per = window_verdicts(d, N_COMBOS, KEY)

print("=" * 74)
print("SCOPE")
print("=" * 74)
print(f"  source table              {SRC.name}")
print(f"  factorial size            {N_COMBOS} combinations")
print(f"  fits (all windows)        {len(raw):,}")
print(f"  windows attempted         {raw.groupby(KEY).ngroups}")
print(f"  windows complete          {len(per)}")
print(f"  countries                 {d['country'].nunique()}")
print(f"  median weeks / cases      {per['weeks'].median():.0f} / "
      f"{per['cases'].median():,.0f}")

print("\n" + "=" * 74)
print("HEADLINE")
print("=" * 74)
print(f"  always favours climate    {per['always_climate'].sum()} "
      f"({per['always_climate'].mean() * 100:.1f}%)")
print(f"  never favours climate     {per['never_climate'].sum()} "
      f"({per['never_climate'].mean() * 100:.1f}%)")
print(f"  VERDICT CHANGES           {per['unstable'].sum()} "
      f"({per['unstable'].mean() * 100:.1f}%)")
split = per.loc[per["unstable"], "wins"]
print(f"  median split (unstable)   {split.median():.0f} of {N_COMBOS}")

# Per-country instability, for the "ranges from x to y" sentence. Countries with
# a single window give a share of 0 or 100 that carries no information, so the
# range is reported over countries with enough windows to mean something.
MIN_WINDOWS = 5
by_country = per.groupby("country")["unstable"].agg(["mean", "size"])
big = by_country[by_country["size"] >= MIN_WINDOWS]
print(f"  per-country range         {big['mean'].min() * 100:.0f}%–"
      f"{big['mean'].max() * 100:.0f}% "
      f"(over {len(big)} countries with >= {MIN_WINDOWS} windows)")

print("\n" + "=" * 74)
print("HOW THE HEADLINE MOVED AS THE DESIGN GREW")
print("=" * 74)
print("  A study that keeps adding factors until its number is large is doing")
print("  the thing it criticises. Every design run is shown, with the")
print("  design-invariant statistic beside the one that grows mechanically.")
print()
print(f"  {'design':>10}  {'combos':>7}  {'windows':>8}  {'unstable':>9}  "
      f"{'P(two disagree)':>16}")
for name in reversed(FACTORIAL_TABLES):
    path = tables / name
    if not path.exists():
        continue
    t = pd.read_csv(path)
    try:
        tc, tn = complete_windows(t, KEY)
    except ValueError:
        continue
    tp = window_verdicts(tc, tn, KEY)
    tpair = pairwise_disagreement(tc, tn, KEY)
    label = name.split("_")[-1].replace(".csv", "")
    if label.startswith("robustness"):
        label = "4factor"
    print(f"  {label:>10}  {tn:7d}  {len(tp):8d}  "
          f"{tp['unstable'].mean() * 100:8.1f}%  {tpair.mean() * 100:15.1f}%")
print()
print("  The first column rises with the number of cells by construction:")
print("  asking whether ANY two disagree is easier the more there are. The")
print("  last does not — it is a function of the proportion, not the count.")
print()
print("  Each factor added closed a gap the previous version named in its own")
print("  limitations, rather than being sought once the result was known.")

print("\n" + "=" * 74)
print("ENDORSEMENT RATE ON REAL DATA, BY OBSERVATION MODEL")
print("=" * 74)
for obs, g in d.groupby("observation"):
    print(f"  {obs:8s} climate wins in {g['climate_wins'].mean() * 100:5.1f}% "
          f"of {len(g):,} fits")

print("\n" + "=" * 74)
print("STORED TABLES THE PAPER QUOTES")
print("=" * 74)
for name, cols in (("15_choice_sensitivity.csv", None),
                   ("16_stability_rules.csv", None),
                   ("17_criteria_comparison.csv", None)):
    path = tables / name
    if not path.exists():
        print(f"  MISSING: {name}")
        continue
    t = pd.read_csv(path)
    print(f"\n  {name}")
    print("   " + t.to_string(index=False, max_rows=12).replace("\n", "\n   "))

# The simulation study is keyed off a different table and window set.
fp = tables / "22_operating_characteristics.csv"
if fp.exists():
    oc = pd.read_csv(fp)
    print("\n  22_operating_characteristics.csv (margins 0, 4, best)")
    for obs in ("poisson", "nb"):
        sub = oc[oc["observation"] == obs]
        best = sub.loc[sub["youden"].idxmax()]
        for _, r in sub[sub["margin"].isin([0, 4, best["margin"]])].iterrows():
            print(f"    {obs:8s} margin {r['margin']:2.0f}  "
                  f"FP {r['false_positive_pct']:5.1f}%  "
                  f"power {r['power_pct']:5.1f}%  "
                  f"separation {r['youden']:5.1f}")

print("\n" + "=" * 74)
print("Read these against paper/paper.tex before submitting.")
print("=" * 74)
