"""Pipeline step 25 — the specification curve.

A multiverse analysis has a canonical figure, and this study did not have it. The
specification curve plots every analysis in a defensible set, sorted by the result
it produces, above a panel showing which choices produced each one. It is how a
reader checks the claim for themselves: whether the analyses that favour one
answer share an identifiable set of choices, or whether the result wanders for no
reason a reader can name.

Two panels are produced, because the question has a within-outbreak form and an
across-outbreak form and they are not the same question.

**Left: one outbreak, every specification.** The window chosen is the one whose
disagreement is most evenly split, since a curve drawn from a window that agrees
with itself shows nothing. This is the figure that makes the headline concrete:
a single dataset, one epidemic, and a set of defensible analyses that do not
agree on whether climate forcing improves the fit.

**Right: all outbreaks, one point per specification.** For each specification,
the share of outbreaks in which it endorses climate forcing. If the choices were interchangeable this would be a flat line. It is not:
the spread between the most and least credulous specification is the size of the
problem, expressed as a rate rather than as a count of flips.

Reads the stored factorial; refits nothing.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (  # noqa: E402
    complete_windows, latest_factorial)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")

# The richest factorial available; see dengue_pk.robustness.FACTORIAL_TABLES.
SRC = latest_factorial(tables)
raw = pd.read_csv(SRC)
KEY = ["country", "unit", "window_start"]
d, N_COMBOS = complete_windows(raw, KEY)
n_windows = d.groupby(KEY).ngroups
print(f"{len(d):,} fits, {n_windows} windows, {N_COMBOS} specifications\n")

# The factors, in the order they are drawn. Levels are read from the data
# rather than written down: a hard-coded list silently mislabels the dot matrix
# when the design grows, and this figure is read by eye, so a wrong label would
# not announce itself the way a wrong number might.
PRETTY = {"observation": "observation model", "temp_form": "temperature",
          "rain_lag": "rainfall lag", "train_frac": "series used",
          "structure": "structure", "params": "fixed parameters"}
FACTORS = [(col, sorted(d[col].unique(), key=str), PRETTY[col])
           for col in PRETTY if col in d.columns]
SPEC = [f for f, _, _ in FACTORS]
print("factors: " + ", ".join(f"{f} ({len(lv)})" for f, lv, _ in FACTORS) + "\n")


def spec_id(frame: pd.DataFrame) -> pd.Series:
    """A stable label for one combination of the analysis choices."""
    return frame[SPEC].astype(str).agg("|".join, axis=1)


d = d.copy()
d["spec"] = spec_id(d)

# ---------------------------------------------------------------------------
# Pick the outbreak to draw: the most evenly split, since that is the case the
# headline is about. Ties are broken by case count so the choice is reproducible
# and lands on a well-observed epidemic rather than an arbitrary one.
wins = d.groupby(KEY).agg(wins=("climate_wins", "sum"),
                          cases=("cases", "first"))
wins["imbalance"] = (wins["wins"] - N_COMBOS / 2).abs()
pick = wins.sort_values(["imbalance", "cases"], ascending=[True, False]).index[0]
one = d.set_index(KEY).loc[pick].reset_index().sort_values("delta_aic")
label = f"{pick[0]} / {pick[1]}, {pick[2]}"
print(f"outbreak drawn: {label}")
print(f"  {int(wins.loc[pick, 'wins'])} of {N_COMBOS} specifications favour climate")
print(f"  delta AIC ranges {one['delta_aic'].min():.1f} to "
      f"{one['delta_aic'].max():.1f}\n")

# ---------------------------------------------------------------------------
# Across outbreaks: how credulous is each specification?
by_spec = (d.groupby("spec")["climate_wins"].mean()
           .sort_values(ascending=False))
print("Most and least credulous specifications, by share of outbreaks endorsed:")
for name, share in list(by_spec.items())[:3]:
    print(f"  {share * 100:5.1f}%   {name}")
print("  ...")
for name, share in list(by_spec.items())[-3:]:
    print(f"  {share * 100:5.1f}%   {name}")
print(f"\n  spread: {by_spec.max() * 100:.1f}% to {by_spec.min() * 100:.1f}% "
      f"({(by_spec.max() - by_spec.min()) * 100:.0f} points between the most and "
      f"least credulous analysis)")

out = pd.DataFrame({"specification": by_spec.index,
                    "share_endorsing_climate": by_spec.round(4).values})
for f in SPEC:
    out[f] = [s.split("|")[SPEC.index(f)] for s in out["specification"]]
out.to_csv(tables / "24_specification_curve.csv", index=False)

# ---------------------------------------------------------------------------
CLIMATE, CONSTANT = "#1f6f8b", "#8c1c13"
fig = plt.figure(figsize=(15, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[2.1, 1.5], hspace=0.08, wspace=0.30)

# --- left: one outbreak ----------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(one))
colours = [CLIMATE if v else CONSTANT for v in one["climate_wins"]]
ax.bar(x, one["delta_aic"], color=colours, width=0.85)
ax.axhline(0, color="k", lw=1)
# Symmetric log, linear inside +/-4. On a linear scale a single bar at -2000
# compresses every specification favouring the constant model into the axis line
# and the panel appears to show unanimity, which is the opposite of its point.
# The linear region is set to exactly the margin the paper recommends, so the
# grey band and the scale change are the same boundary.
ax.set_yscale("symlog", linthresh=4, linscale=0.9)
ax.axhspan(-4, 4, color="0.75", alpha=0.45, zorder=0)
ax.set_ylabel("ΔAIC  (climate − constant)\nsymmetric log, linear within ±4")
ax.set_title(f"One outbreak, {N_COMBOS} defensible analyses\n{label}", fontsize=11)
ax.set_xlim(-0.7, len(one) - 0.3)
ax.set_xticks([])
ax.grid(alpha=0.3, axis="y")
# A legend rather than annotations on the axes: at this aspect ratio every
# in-axes position collides with a bar. The empty lower-right quadrant is empty
# precisely because those specifications favour the constant model by small
# margins, so it stays empty whatever the data does.
handles = [Patch(facecolor=CLIMATE, label="climate model preferred"),
           Patch(facecolor=CONSTANT, label="constant model preferred"),
           Patch(facecolor="0.75", alpha=0.45,
                 label="|ΔAIC| < 4: not distinguishable")]
ax.legend(handles=handles, fontsize=8, loc="lower right", framealpha=0.92)

axl = fig.add_subplot(gs[1, 0])
rows, ylabels = [], []
for factor, settings, pretty in FACTORS:
    for setting in settings:
        rows.append((factor, setting))
        ylabels.append(f"{pretty}: {setting}")
for r, (factor, setting) in enumerate(rows):
    on = one[factor].astype(str).to_numpy() == str(setting)
    axl.scatter(x[on], np.full(on.sum(), r), s=13, marker="s",
                color=[CLIMATE if v else CONSTANT
                       for v in one["climate_wins"].to_numpy()[on]])
axl.set_yticks(range(len(rows)))
axl.set_yticklabels(ylabels, fontsize=7.5)
axl.set_xlim(-0.7, len(one) - 0.3)
axl.set_ylim(-0.8, len(rows) - 0.2)
axl.invert_yaxis()
axl.set_xticks([])
axl.set_xlabel(f"the {N_COMBOS} analyses, ordered by the result they produce",
               fontsize=9)
for r in range(len(rows)):
    axl.axhline(r, color="0.9", lw=0.6, zorder=0)

# --- right: every outbreak, one point per specification --------------------
ax2 = fig.add_subplot(gs[0, 1])
order = by_spec.index.tolist()
x2 = np.arange(len(order))
vals = by_spec.to_numpy()
is_poisson = np.array([s.split("|")[0] == "poisson" for s in order])
ax2.bar(x2[is_poisson], vals[is_poisson] * 100, color=CONSTANT, width=0.85,
        label="Poisson")
ax2.bar(x2[~is_poisson], vals[~is_poisson] * 100, color=CLIMATE, width=0.85,
        label="negative binomial")
ax2.axhline(50, color="k", lw=0.8, ls=":")
ax2.set_ylabel("% of outbreaks in which this analysis\nendorses climate forcing")
ax2.set_title(f"Each of the {N_COMBOS} analyses, applied to all "
              f"{n_windows} outbreaks", fontsize=11)
ax2.set_xlim(-0.7, len(order) - 0.3)
ax2.set_xticks([])
ax2.legend(fontsize=8.5, loc="lower left")
ax2.grid(alpha=0.3, axis="y")

axr = fig.add_subplot(gs[1, 1])
for r, (factor, setting) in enumerate(rows):
    on = np.array([s.split("|")[SPEC.index(factor)] == str(setting)
                   for s in order])
    axr.scatter(x2[on], np.full(on.sum(), r), s=13, marker="s",
                color=[CONSTANT if p else CLIMATE for p in is_poisson[on]])
axr.set_yticks(range(len(rows)))
axr.set_yticklabels(ylabels, fontsize=7.5)
axr.set_xlim(-0.7, len(order) - 0.3)
axr.set_ylim(-0.8, len(rows) - 0.2)
axr.invert_yaxis()
axr.set_xticks([])
axr.set_xlabel(f"the {N_COMBOS} analyses, ordered by how often they endorse climate forcing",
               fontsize=9)
for r in range(len(rows)):
    axr.axhline(r, color="0.9", lw=0.6, zorder=0)

fig.savefig(figures / "17_specification_curve.png", dpi=150, bbox_inches="tight")
print(f"\nFigure: {figures / '17_specification_curve.png'}")
print(f"Table:  {tables / '24_specification_curve.csv'}")
