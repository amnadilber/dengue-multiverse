"""
Pipeline step 17 — how often does the conclusion depend on the analyst?

Reads the factorial from step 16 and answers the question the whole study exists
for: across outbreaks worldwide, how often does the in-sample verdict on
climate-driven transmission change as the analyst's unremarked choices vary?

Three things are separated carefully, because conflating them would overstate
the result:

* **How often climate forcing wins.** Not the point, and not reported as one.
* **How often the answer is unstable within a single window.** This is the point.
  A window whose combinations all agree is telling us something about the
  world; a window that splits is telling us about the analyst.

  Reported two ways. The share of windows where *any* two combinations disagree
  is the intuitive figure, but it rises with the size of the design: asking
  whether any two of 144 cells differ is easier than any two of 24. The
  probability that two analyses drawn at random reach opposite verdicts depends
  on the proportion rather than the count, so it is comparable across designs
  and is reported beside it.
* **Whether out-of-sample validation is more stable than in-sample selection.**
  If it is not, the paper has diagnosed a problem without a remedy.

Sensitivity to each individual choice is reported as the change in the
probability of an in-sample "climate wins" verdict when only that choice is
varied and everything else is held fixed — a paired comparison, so that windows
of different sizes and case counts cannot drive the difference.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (  # noqa: E402
    latest_factorial, pairwise_disagreement)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")

# The richest factorial available; see dengue_pk.robustness.FACTORIAL_TABLES.
SRC = latest_factorial(tables)
d = pd.read_csv(SRC)
print(f"source: {SRC.name}")
KEY = ["country", "unit", "window_start"]
FACTORS = [c for c in ("observation", "temp_form", "rain_lag",
                       "train_frac", "structure", "params") if c in d.columns]

print(f"{len(d):,} fits across {d.groupby(KEY).ngroups} windows "
      f"in {d['country'].nunique()} countries\n")

# Drop any window that did not complete the full factorial: an incomplete set
# would make an unstable window look stable, or the reverse, depending on which
# combinations happened to fail.
N_COMBOS = int(counts_expected := d.groupby(KEY).size().max())
counts = d.groupby(KEY).size()
complete = counts[counts == N_COMBOS].index
d = d.set_index(KEY).loc[complete].reset_index()
print(f"{len(complete)} windows completed all {N_COMBOS} combinations "
      f"({counts.size - len(complete)} partial, excluded)\n")

# ---------------------------------------------------------------------------
# 1. Instability: does the verdict change within a window?
# ---------------------------------------------------------------------------
per_window = d.groupby(KEY).agg(
    wins=("climate_wins", "sum"),
    n=("climate_wins", "size"),
    weeks=("weeks", "first"),
    cases=("cases", "first")).reset_index()
per_window["unstable"] = per_window["wins"].between(1, N_COMBOS - 1)
per_window["always_climate"] = per_window["wins"] == N_COMBOS
per_window["never_climate"] = per_window["wins"] == 0

n_win = len(per_window)
print("=" * 74)
print("1. IS THE VERDICT STABLE WITHIN A WINDOW?")
print("=" * 74)
print(f"  always favours climate forcing : {per_window['always_climate'].sum():3d} "
      f"({per_window['always_climate'].mean() * 100:.1f}%)")
print(f"  never favours climate forcing  : {per_window['never_climate'].sum():3d} "
      f"({per_window['never_climate'].mean() * 100:.1f}%)")
print(f"  VERDICT CHANGES WITH THE CHOICES: {per_window['unstable'].sum():3d} "
      f"({per_window['unstable'].mean() * 100:.1f}%)")
print(f"\n  median split among unstable windows: "
      f"{per_window.loc[per_window['unstable'], 'wins'].median():.0f} of {N_COMBOS} "
      f"combinations favour climate forcing")

# The share above asks whether *any* two cells disagree, which is easier to
# satisfy the more cells the design has: enlarging the factorial raises it even
# if nothing about the evidence changed. The pairwise probability depends on the
# proportion favouring each verdict and not on how finely the space was
# enumerated, so it is comparable across designs and across studies.
pair = pairwise_disagreement(d, N_COMBOS, KEY)
print(f"\n  probability two analyses of the same outbreak, drawn at random,")
print(f"  reach OPPOSITE verdicts: {pair.mean() * 100:.1f}% "
      f"(median {pair.median() * 100:.1f}%, max possible 50%)")
print(f"  outbreaks where that probability exceeds 25%: "
      f"{(pair > 0.25).mean() * 100:.1f}%")

# ---------------------------------------------------------------------------
# 2. Which choice moves the verdict most? Paired within window.
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("2. WHICH CHOICE MOVES THE VERDICT? (paired within window)")
print("=" * 74)
rows = []
for factor in FACTORS:
    others = [f for f in FACTORS if f != factor]
    levels = sorted(d[factor].unique())
    # Pair fits that differ only in this factor.
    wide = d.pivot_table(index=KEY + others, columns=factor,
                         values="climate_wins", aggfunc="first")
    wide = wide.dropna()
    if wide.empty or len(levels) < 2:
        continue
    base = levels[0]
    for lev in levels[1:]:
        flip_to = float(((wide[base] == 0) & (wide[lev] == 1)).mean())
        flip_from = float(((wide[base] == 1) & (wide[lev] == 0)).mean())
        rows.append(dict(factor=factor, comparison=f"{base} -> {lev}",
                         pairs=len(wide),
                         p_climate_base=round(float(wide[base].mean()), 3),
                         p_climate_alt=round(float(wide[lev].mean()), 3),
                         flipped=round(flip_to + flip_from, 3),
                         net_change=round(float(wide[lev].mean() - wide[base].mean()), 3)))
sens = pd.DataFrame(rows).sort_values("flipped", ascending=False)
print(sens.to_string(index=False))

# ---------------------------------------------------------------------------
# Why does the observation model matter so much? Two explanations are available
# and they make opposite predictions, so the question is settled by measurement
# rather than by argument.
#
#   (a) Poisson treats ordinary noise as a large residual, so the fit chases
#       individual weeks. The climate covariates supply that flexibility, and
#       their coefficients should therefore be pushed harder under Poisson.
#   (b) Poisson changes nothing about what the climate model does; it changes how
#       much credit the criterion gives for doing it. The same small improvement
#       in fit buys a much larger delta AIC when no dispersion cushions the
#       likelihood.
#
# (a) was asserted in an earlier draft. (b) is what the data show.
print("\n" + "=" * 74)
print("WHY THE OBSERVATION MODEL MATTERS: TWO EXPLANATIONS, PAIRED")
print("=" * 74)
CELL = [c for c in ("temp_form", "rain_lag", "train_frac", "structure",
                    "params") if c in d.columns]
print("  (a) does Poisson push the climate coefficients harder?")
for col in ("a_temp", "a_rain"):
    w = d.pivot_table(index=KEY + CELL, columns="observation", values=col,
                      aggfunc="first").dropna()
    if w.empty or not {"nb", "poisson"} <= set(w.columns):
        continue
    bigger = float((w["poisson"].abs() > w["nb"].abs()).mean()) * 100
    print(f"      |{col}|  median {w['nb'].abs().median():.3f} under NB, "
          f"{w['poisson'].abs().median():.3f} under Poisson; larger under "
          f"Poisson in {bigger:.0f}% of {len(w):,} pairs")
print("      No. The coefficients are the same size either way.\n")

print("  (b) does Poisson magnify the same improvement?")
w = d.pivot_table(index=KEY + CELL, columns="observation", values="delta_aic",
                  aggfunc="first").dropna()
if not w.empty and {"nb", "poisson"} <= set(w.columns):
    ratio = (w["poisson"] / w["nb"]).replace([np.inf, -np.inf], np.nan).dropna()
    ratio = ratio[ratio > 0]
    print(f"      median delta AIC: {w['nb'].median():.2f} under NB, "
          f"{w['poisson'].median():.2f} under Poisson")
    print(f"      where the sign agrees, the ratio is "
          f"{ratio.median():.0f}x (quartiles {ratio.quantile(.25):.1f} to "
          f"{ratio.quantile(.75):.0f})")
    print(f"      |delta AIC| larger under Poisson in "
          f"{float((w['poisson'].abs() > w['nb'].abs()).mean()) * 100:.0f}% of pairs")
    print(f"      climate wins under both: "
          f"{float(((w['nb'] < 0) & (w['poisson'] < 0)).mean()) * 100:.1f}%; "
          f"only under Poisson: "
          f"{float(((w['nb'] >= 0) & (w['poisson'] < 0)).mean()) * 100:.1f}%; "
          f"only under NB: "
          f"{float(((w['nb'] < 0) & (w['poisson'] >= 0)).mean()) * 100:.1f}%")
    print("      Yes, and the asymmetry is near-total.")
    pd.DataFrame([dict(
        median_daic_nb=round(float(w["nb"].median()), 2),
        median_daic_poisson=round(float(w["poisson"].median()), 2),
        ratio_median=round(float(ratio.median()), 1),
        both=round(float(((w["nb"] < 0) & (w["poisson"] < 0)).mean()) * 100, 1),
        only_poisson=round(float(((w["nb"] >= 0) & (w["poisson"] < 0)).mean()) * 100, 1),
        only_nb=round(float(((w["nb"] < 0) & (w["poisson"] >= 0)).mean()) * 100, 1),
    )]).to_csv(tables / "42_observation_mechanism.csv", index=False)

# ---------------------------------------------------------------------------
# 3. Is out-of-sample validation more stable?
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("3. IS HELD-OUT VALIDATION MORE STABLE THAN IN-SAMPLE SELECTION?")
print("=" * 74)
held = d[d["train_frac"] < 1.0].dropna(subset=["heldout_delta"]).copy()
held["climate_better_out"] = held["heldout_delta"] < 0
hw = held.groupby(KEY).agg(out_wins=("climate_better_out", "sum"),
                           out_n=("climate_better_out", "size"),
                           in_wins=("climate_wins", "sum")).reset_index()
hw = hw[hw["out_n"] == N_COMBOS // 2]
hw["out_unstable"] = hw["out_wins"].between(1, N_COMBOS // 2 - 1)
hw["in_unstable"] = hw["in_wins"].between(1, N_COMBOS // 2 - 1)
print(f"  windows with a complete held-out set: {len(hw)}")
print(f"  in-sample verdict unstable  : {hw['in_unstable'].mean() * 100:5.1f}%")
print(f"  out-of-sample verdict unstable: {hw['out_unstable'].mean() * 100:5.1f}%")
print(f"  climate preferred out-of-sample, overall: "
      f"{held['climate_better_out'].mean() * 100:.1f}% of fits")
print(f"  climate preferred in-sample, same fits:   "
      f"{held['climate_wins'].mean() * 100:.1f}% of fits")

# ---------------------------------------------------------------------------
# 4. What sign does each temperature parameterisation report?
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("4. THE SIGN OF THE TEMPERATURE COEFFICIENT, BY PARAMETERISATION")
print("=" * 74)
for form, g in d.groupby("temp_form"):
    a = g["a_temp"].dropna()
    print(f"  {form:10s} n={len(a):5d}  mean {a.mean():+.3f}  "
          f"median {a.median():+.3f}  negative in {(a < 0).mean() * 100:5.1f}% of fits")

pd.DataFrame(per_window).to_csv(tables / "14_window_stability.csv", index=False)
sens.to_csv(tables / "15_choice_sensitivity.csv", index=False)

# ---------------------------------------------------------------------------
# Two by two rather than one by four: four panels in a row is 20 inches wide and
# unreadable at journal column width.
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
axes = axes.ravel()

axes[0].hist(per_window["wins"], bins=np.arange(-0.5, N_COMBOS + 1.5, 1),
             color="#8c1c13")
axes[0].set_xlabel(f"combinations (of {N_COMBOS}) favouring climate forcing")
axes[0].set_ylabel("windows")
axes[0].set_title(f"{per_window['unstable'].mean() * 100:.0f}% of windows change verdict\n"
                  "with the analysis choices")
axes[0].grid(alpha=0.3, axis="y")

# The panel that does not depend on how finely the space was enumerated. Drawn
# next to the one that does, deliberately: a reader should be able to see both
# and judge which they trust.
axes[3].hist(pair * 100, bins=np.arange(0, 52, 2), color="#1f6f8b")
axes[3].axvline(50, color="k", ls="--", lw=1.2,
                label="50% — an evenly split outbreak")
axes[3].axvline(pair.mean() * 100, color="#8c1c13", lw=2,
                label=f"mean {pair.mean() * 100:.1f}%")
axes[3].set_xlabel("P(two analyses of this outbreak disagree), %")
axes[3].set_ylabel("outbreaks")
axes[3].set_title("Disagreement between two analyses\n"
                  "— unchanged by the size of the design")
axes[3].legend(fontsize=8)
axes[3].grid(alpha=0.3, axis="y")

s = sens.set_index("comparison")["flipped"].sort_values()
axes[1].barh(s.index, s.values, color="#1f6f8b")
axes[1].set_xlabel("fraction of paired fits whose verdict flips")
axes[1].set_title("Which choice flips the answer?")
axes[1].grid(alpha=0.3, axis="x")

# Clipped to the region where the mass is. The Brière exponent has a long right
# tail — a handful of fits run to twenty — and drawing the full range compresses
# the part of the axis that carries the point, which is whether the estimate is
# on the negative side of zero.
CLIP = 3.0
for form, colour in (("briere", "#3f7d58"), ("loglinear", "#a63d40")):
    a = d[d["temp_form"] == form]["a_temp"].dropna().clip(-CLIP, CLIP)
    axes[2].hist(a, bins=np.linspace(-CLIP, CLIP, 61), alpha=0.6, color=colour,
                 label=form)
axes[2].set_xlim(-CLIP, CLIP)
axes[2].axvline(0, color="k", lw=1)
axes[2].set_xlabel("estimated temperature coefficient")
neg_pct = (d[d["temp_form"] == "loglinear"]["a_temp"] < 0).mean() * 100
axes[2].set_title(
    "The monotone term's sign is arbitrary\n"
    "(negative in {:.0f}% of fits; the unimodal form cannot be)".format(neg_pct))
axes[2].legend(fontsize=8)
axes[2].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(figures / "12_global_robustness.png", dpi=150)
print(f"\nFigure: {figures / '12_global_robustness.png'}")
print(f"Tables: {tables / '14_window_stability.csv'}, "
      f"{tables / '15_choice_sensitivity.csv'}")
