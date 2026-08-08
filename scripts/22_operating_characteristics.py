"""
Pipeline step 22 — the operating characteristics, as a figure and a table.

Steps 21 and 29 measured false positives and power for each combination of
observation model and decision rule. This presents them the way a reader has to
see them to act on them: the two error rates together, since either alone is
misleading.

A rule with 95% power sounds excellent until its false-positive rate is 86%, at
which point it is not detecting anything — it is agreeing with everything. The
figure plots the pair so that this is impossible to miss.

Built from whichever simulation table is newest; nothing is refitted.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")

# The newest simulation table available. Step 29 supersedes step 21 with the
# six-factor design; naming only the old file here would have redrawn the
# paper's figure from superseded results while every number beside it moved.
SIM_TABLES = ("28_false_positive_6factor.csv", "21_false_positive_rate.csv")
SRC = next((tables / n for n in SIM_TABLES if (tables / n).exists()), None)
if SRC is None:
    raise FileNotFoundError(f"no simulation table; expected one of {SIM_TABLES}")
print(f"source: {SRC.name}")
d = pd.read_csv(SRC)
null = d[d["truth"] == "no_effect"]
alt = d[d["truth"] == "climate"]
print(f"{len(d):,} fits: {len(null):,} under the null, {len(alt):,} under a real effect")
print(f"{d.groupby(['country', 'unit', 'window_start']).ngroups} windows\n")

MARGINS = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24]
rows = []
for obs in ("poisson", "nb"):
    n_sub = null[null["observation"] == obs]
    a_sub = alt[alt["observation"] == obs]
    for m in MARGINS:
        fp = float((n_sub["delta_aic"] < -m).mean())
        pw = float((a_sub["delta_aic"] < -m).mean())
        rows.append(dict(observation=obs, margin=m,
                         false_positive_pct=round(fp * 100, 1),
                         power_pct=round(pw * 100, 1),
                         youden=round((pw - fp) * 100, 1)))
oc = pd.DataFrame(rows)
oc.to_csv(tables / "22_operating_characteristics.csv", index=False)

print("=" * 72)
print("OPERATING CHARACTERISTICS")
print("=" * 72)
for obs in ("poisson", "nb"):
    sub = oc[oc["observation"] == obs]
    best = sub.loc[sub["youden"].idxmax()]
    print(f"\n  {obs}:")
    print(f"    {'margin':>7} {'false pos':>10} {'power':>8} {'power-FP':>10}")
    for _, r in sub.iterrows():
        mark = "  <-- best separation" if r["margin"] == best["margin"] else ""
        print(f"    {r['margin']:7.0f} {r['false_positive_pct']:9.1f}% "
              f"{r['power_pct']:7.1f}% {r['youden']:9.1f}{mark}")

# ---------------------------------------------------------------------------
# Does a factorial manufacture instability by itself?
#
# The sharpest objection to the headline is arithmetical rather than
# epidemiological: run enough analyses and something will always disagree, so a
# high instability rate might say more about the size of the design than about
# the evidence. The simulation answers it directly, because here the truth is known.
# If instability were an artefact of the design, it would be as high on data
# containing a large, unambiguous climate effect as on data containing none.
print("\n" + "=" * 72)
print("IS THE INSTABILITY AN ARTEFACT OF RUNNING THE FACTORIAL?")
print("=" * 72)
KEY = ["country", "unit", "window_start"]
for truth, label in (("climate", "a real climate effect"),
                     ("no_effect", "no climate effect")):
    sub = d[d["truth"] == truth]
    if sub.empty:
        continue
    counts = sub.groupby(KEY).size()
    n_combos = int(counts.max())
    full = sub.set_index(KEY).loc[counts[counts == n_combos].index].reset_index()
    wins = full.groupby(KEY)["climate_wins"].sum()
    unstable = ((wins >= 1) & (wins <= n_combos - 1)).mean()
    print(f"  data generated with {label:22s} unstable in "
          f"{unstable * 100:5.1f}% of {len(wins)} windows")
print("\n  The factorial is the same in both rows. Only the truth differs, so")
print("  any gap between them is evidence about the data rather than about")
print("  the number of analyses run.")

p_sign = oc[(oc.observation == "poisson") & (oc.margin == 0)].iloc[0]
n_sign = oc[(oc.observation == "nb") & (oc.margin == 0)].iloc[0]
n_m4 = oc[(oc.observation == "nb") & (oc.margin == 4)].iloc[0]
print("\n" + "=" * 72)
print("READING")
print("=" * 72)
print(f"  Poisson with the conventional rule separates the two truths by only "
      f"{p_sign['youden']:.0f} points\n  ({p_sign['power_pct']:.0f}% power against "
      f"{p_sign['false_positive_pct']:.0f}% false positives): it is close to "
      f"endorsing\n  climate forcing regardless of whether any exists.")
print(f"\n  The negative binomial separates them by {n_sign['youden']:.0f} points "
      f"at the same rule,\n  and by {n_m4['youden']:.0f} points with a margin of 4 "
      f"({n_m4['power_pct']:.0f}% power, "
      f"{n_m4['false_positive_pct']:.1f}% false positives).")

# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

ax = axes[0]
for obs, colour in (("poisson", "#8c1c13"), ("nb", "#3f7d58")):
    sub = oc[oc["observation"] == obs]
    ax.plot(sub["false_positive_pct"], sub["power_pct"], "o-", color=colour,
            label=f"{obs} (margin 0 → 24)")
    for _, r in sub.iterrows():
        if r["margin"] in (0, 4, 24):
            ax.annotate(f"{r['margin']:.0f}",
                        (r["false_positive_pct"], r["power_pct"]),
                        fontsize=8, xytext=(4, -9), textcoords="offset points",
                        color=colour)
ax.plot([0, 100], [0, 100], "k:", lw=1, label="no discrimination")
ax.set_xlabel("false positives (%) — climate endorsed with no effect present")
ax.set_ylabel("power (%) — real effect detected")
ax.set_title("Poisson barely separates the two truths")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)

ax = axes[1]
labels, fps, pws = [], [], []
for obs, m in (("poisson", 0), ("poisson", 4), ("nb", 0), ("nb", 4)):
    r = oc[(oc.observation == obs) & (oc.margin == m)].iloc[0]
    labels.append(f"{obs}\n{'sign of ΔAIC' if m == 0 else f'margin {m}'}")
    fps.append(r["false_positive_pct"])
    pws.append(r["power_pct"])
x = np.arange(len(labels))
ax.bar(x - 0.2, fps, 0.4, color="#8c1c13", label="false positives")
ax.bar(x + 0.2, pws, 0.4, color="#3f7d58", label="power")
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("% of fits")
ax.set_title("The observation model is the fix;\nthe margin is a second safeguard")
ax.legend(fontsize=8)
ax.grid(alpha=0.3, axis="y")

fig.tight_layout()
fig.savefig(figures / "16_operating_characteristics.png", dpi=150)
print(f"\nFigure: {figures / '16_operating_characteristics.png'}")
print(f"Table:  {tables / '22_operating_characteristics.csv'}")
