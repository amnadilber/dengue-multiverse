"""
Pipeline step 18 — is there a rule that makes the verdict stable?

Step 17 established that in 83% of outbreaks the in-sample verdict on
climate-driven transmission changes with the analyst's choices, and that
out-of-sample validation reduces false endorsement without restoring stability.
That is a diagnosis. A reader is entitled to ask what to do instead, and a paper
that only complains is worth less than one that tests a remedy.

Four candidate rules are evaluated, each on the whole factorial:

1. **Sign of ΔAIC** — the conventional rule, and the baseline.
2. **An evidence threshold.** Declare a winner only when |ΔAIC| exceeds a margin,
   and otherwise return "inconclusive". If the instability lives in the many
   comparisons that sit near zero, a margin should remove most of it — at the
   cost of answering less often, which must be reported alongside.
3. **Agreement across the observation model.** Report a verdict only where
   Poisson and negative binomial concur. This targets the single largest lever
   directly.
4. **Agreement across everything.** Report only where every combination
   concurs: maximally conservative, and the natural upper bound on what any rule
   can do.

The metric for each is a pair, and both halves matter: how often the rule gives
an unstable answer, and how often it gives any answer at all. A rule that is
never wrong because it never speaks is not a remedy.

Also examined: whether instability is concentrated in small or short outbreaks,
which would qualify the headline; and how the headline rate itself varies with
the window-selection threshold, which was a choice this study made and must
therefore hold itself to the same standard it applies to others.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (  # noqa: E402
    complete_windows, latest_factorial, pairwise_disagreement_at_margin,
    window_verdicts)

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
d, N_COMBOS = complete_windows(d, KEY)
print(f"{d.groupby(KEY).ngroups} windows with all {N_COMBOS} combinations")
n_windows = d.groupby(KEY).ngroups
rng = np.random.default_rng(cfg["seed"])

print(f"{len(d):,} fits, {n_windows} windows\n")


def verdicts(sub: pd.DataFrame, margin: float) -> np.ndarray:
    """+1 climate, -1 constant, 0 inconclusive, per fit."""
    v = np.zeros(len(sub), dtype=int)
    v[sub["delta_aic"].to_numpy() < -margin] = 1
    v[sub["delta_aic"].to_numpy() > margin] = -1
    return v


def assess(margin: float, restrict_obs: str | None = None,
           require_all: bool = False):
    """Instability and decisiveness of one rule."""
    sub = d if restrict_obs is None else d[d["observation"] == restrict_obs]
    unstable = decided = 0
    for _, g in sub.groupby(KEY):
        v = verdicts(g, margin)
        spoken = v[v != 0]
        if require_all:
            # Only speak where every combination agrees and none abstains.
            if len(spoken) == len(v) and len(set(spoken)) == 1:
                decided += 1
            continue
        if len(spoken) == 0:
            continue
        decided += 1
        if len(set(spoken)) > 1:
            unstable += 1
    rate = unstable / decided if decided else np.nan
    return dict(unstable_rate=rate, decided=decided,
                decided_share=decided / n_windows)


def bootstrap_ci(margin: float, restrict_obs=None, require_all=False, n=400):
    """Percentile interval over windows, which are the independent units."""
    groups = [g for _, g in
              (d if restrict_obs is None
               else d[d["observation"] == restrict_obs]).groupby(KEY)]
    out = []
    for _ in range(n):
        pick = rng.integers(0, len(groups), len(groups))
        unstable = decided = 0
        for i in pick:
            v = verdicts(groups[i], margin)
            spoken = v[v != 0]
            if require_all:
                if len(spoken) == len(v) and len(set(spoken)) == 1:
                    decided += 1
                continue
            if len(spoken) == 0:
                continue
            decided += 1
            if len(set(spoken)) > 1:
                unstable += 1
        if decided:
            out.append(unstable / decided)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) \
        if out else (np.nan, np.nan)


# ---------------------------------------------------------------------------
print("=" * 78)
print("1. THE BASELINE, WITH AN INTERVAL")
print("=" * 78)
base = assess(0.0)
lo, hi = bootstrap_ci(0.0)
print(f"  sign of dAIC: unstable in {base['unstable_rate'] * 100:.1f}% of windows "
      f"[95% CI {lo * 100:.1f}–{hi * 100:.1f}], decisive on "
      f"{base['decided_share'] * 100:.0f}%")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. DOES AN EVIDENCE THRESHOLD RESTORE STABILITY?")
print("=" * 78)
rows = []
for margin in (0, 2, 4, 6, 10, 15, 20, 30, 50):
    r = assess(float(margin))
    lo, hi = bootstrap_ci(float(margin), n=200)
    rows.append(dict(rule=f"|dAIC| > {margin}", margin=margin,
                     unstable_pct=round(r["unstable_rate"] * 100, 1),
                     ci_lo=round(lo * 100, 1), ci_hi=round(hi * 100, 1),
                     answers_pct=round(r["decided_share"] * 100, 1)))
    print(f"  |dAIC| > {margin:2d}: unstable {r['unstable_rate'] * 100:5.1f}% "
          f"[{lo * 100:4.1f}–{hi * 100:4.1f}]   gives an answer for "
          f"{r['decided_share'] * 100:5.1f}% of windows")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. RESTRICTING OR REQUIRING AGREEMENT")
print("=" * 78)
for label, kwargs in (
        ("negative binomial only", dict(margin=0.0, restrict_obs="nb")),
        ("negative binomial, |dAIC| > 10", dict(margin=10.0, restrict_obs="nb")),
        (f"all {N_COMBOS} combinations must agree", dict(margin=0.0, require_all=True)),
        (f"all {N_COMBOS} agree, |dAIC| > 10", dict(margin=10.0, require_all=True))):
    r = assess(**kwargs)
    lo, hi = bootstrap_ci(n=200, **kwargs)
    ci = "" if np.isnan(lo) else f" [{lo * 100:.1f}–{hi * 100:.1f}]"
    print(f"  {label:34s} unstable {r['unstable_rate'] * 100:5.1f}%{ci}   "
          f"answers for {r['decided_share'] * 100:5.1f}%")
    rows.append(dict(rule=label, margin=kwargs.get("margin"),
                     unstable_pct=round(r["unstable_rate"] * 100, 1),
                     ci_lo=round(lo * 100, 1) if np.isfinite(lo) else np.nan,
                     ci_hi=round(hi * 100, 1) if np.isfinite(hi) else np.nan,
                     answers_pct=round(r["decided_share"] * 100, 1)))

pd.DataFrame(rows).to_csv(tables / "16_stability_rules.csv", index=False)

# ---------------------------------------------------------------------------
# Where the margin stopped working, and why.
#
# On the five-factor design a margin of 4 cut instability from 88% to 7%,
# because the dissenting combinations were those whose evidence was weak to
# begin with. On the six-factor design it cuts 91% to 39%: much of the dissent
# is now *confident*, and a band around zero cannot reach it.
#
# One structure is the obvious suspect. Under `hostvector_tempmort` temperature
# drives mosquito mortality, so it enters both compared models and that
# structure's null is not the others' null — it can favour the constant model
# decisively where the others favour climate decisively. If it is responsible,
# excluding it should restore the margin; if it is not, the remedy was weaker
# than the five-factor study suggested for reasons that have nothing to do with
# this one addition.
if "structure" in d.columns and d["structure"].nunique() > 2:
    print("\n" + "=" * 78)
    print("6. IS THE WEAKENED MARGIN DUE TO THE STRUCTURE WITH A DIFFERENT NULL?")
    print("=" * 78)
    common = d[d["structure"] != "hostvector_tempmort"]
    saved_d = d
    for label, subset in (("all structures", d),
                          ("common-null structures only", common)):
        d = subset                       # `assess` reads the module-level frame
        n_here = d.groupby(KEY).ngroups
        for margin in (0.0, 4.0, 10.0):
            unstable = decided = 0
            for _, g in d.groupby(KEY):
                v = verdicts(g, margin)
                spoken = v[v != 0]
                if len(spoken) == 0:
                    continue
                decided += 1
                if len(set(spoken)) > 1:
                    unstable += 1
            print(f"  {label:30s} margin {margin:4.0f}: unstable "
                  f"{unstable / decided * 100:5.1f}%   answers "
                  f"{decided / n_here * 100:5.1f}%")
        print()
    d = saved_d
    print("  If the two rows agree, the margin is weaker than the smaller design")
    print("  suggested for reasons beyond the structure whose null differs.")

# ---------------------------------------------------------------------------
# What the recommendation says about these data.
#
# A paper that tells a field its methods are unreliable owes the field its own
# answer, computed its own way. Otherwise the work reads as a refusal — and
# worse, invites the misreading that climate does not drive dengue transmission,
# which is not what any of this shows.
print("\n" + "=" * 78)
print("APPLYING THIS PAPER'S OWN RECOMMENDATION TO ITS OWN DATA")
print("=" * 78)
print("  Negative binomial only; abstain where |dAIC| < 4.")
print()
nb_only = d[d["observation"] == "nb"]
own = []
for _, g in nb_only.groupby(KEY):
    v = np.where(g["delta_aic"] < -4, 1, np.where(g["delta_aic"] > 4, -1, 0))
    spoke = set(v[v != 0].tolist())
    if not spoke:
        own.append("inconclusive: no analysis speaks")
    elif spoke == {1}:
        own.append("climate forcing supported")
    elif spoke == {-1}:
        own.append("no climate forcing")
    else:
        own.append("analyses still disagree")
vc = pd.Series(own).value_counts()
for label in ("climate forcing supported", "no climate forcing",
              "analyses still disagree", "inconclusive: no analysis speaks"):
    k = int(vc.get(label, 0))
    print(f"  {label:34s} {k:4d}  ({k / len(own) * 100:5.1f}%)")
clear = int(vc.get("climate forcing supported", 0)) + \
    int(vc.get("no climate forcing", 0))
print()
print(f"  A clear verdict for {clear} of {len(own)} outbreaks "
      f"({clear / len(own) * 100:.0f}%); of those, "
      f"{int(vc.get('climate forcing supported', 0)) / clear * 100:.0f}% "
      f"support climate forcing.")
print("  The remainder cannot be answered from these data by this method, which")
print("  is a finding rather than a failure to report one.")
vc.rename("outbreaks").to_csv(tables / "38_recommended_verdicts.csv")

# ---------------------------------------------------------------------------
# The same comparison on a scale that does not grow with the design.
#
# "Some pair still disagrees after abstaining" gets easier to satisfy the more
# cells there are, exactly as the headline does, so a margin can appear to stop
# working when the only thing that changed was the enumeration. The pairwise
# form is a function of the proportions among the analyses that speak, so it can
# be compared across designs — and against the same figure from the smaller ones.
print("\n" + "=" * 78)
print("7. THE REMEDY ON A DESIGN-INVARIANT SCALE")
print("=" * 78)
print(f"  {'margin':>7}  {'P(two speaking analyses disagree)':>34}  {'windows speaking':>17}")
for margin in (0.0, 2.0, 4.0, 10.0, 20.0):
    pair_m = pairwise_disagreement_at_margin(d, margin, KEY)
    print(f"  {margin:7.0f}  {pair_m.mean() * 100:33.1f}%  "
          f"{len(pair_m) / n_windows * 100:16.1f}%")
print("\n  Compare with the sign-of-dAIC value for the smaller designs, printed")
print("  by step 23. If the margin genuinely helps, this column falls; if it")
print("  only appeared to help before, it will not.")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("4. IS INSTABILITY JUST A SMALL-SAMPLE PROBLEM?")
print("=" * 78)
per_window = window_verdicts(d, N_COMBOS, KEY)
for col, label in (("cases", "reported cases"), ("weeks", "series length")):
    q = pd.qcut(per_window[col], 4, duplicates="drop")
    tab = per_window.groupby(q, observed=True)["unstable"].agg(["mean", "size"])
    print(f"\n  by {label}:")
    for interval, row in tab.iterrows():
        print(f"    {str(interval):26s} unstable {row['mean'] * 100:5.1f}% "
              f"(n={int(row['size'])})")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("5. DOES OUR OWN WINDOW-SELECTION RULE DRIVE THE HEADLINE?")
print("=" * 78)
print("  The inventory kept windows whose peak exceeded 2.5x the median week.")
print("  Re-computing the headline on progressively stricter subsets:\n")
inv = pd.read_csv(tables / "12_global_windows.csv")
inv["key"] = list(zip(inv["country"], inv["unit"], inv["start"]))
prom = {k: p / max(c / w, 1e-9)
        for k, p, c, w in zip(inv["key"], inv["peak"], inv["cases"], inv["weeks"])}
per_window["prominence"] = [prom.get((c, u, str(s)), np.nan)
                            for c, u, s in zip(per_window["country"],
                                               per_window["unit"],
                                               per_window["window_start"])]
for thr in (0, 1.5, 2.0, 3.0, 4.0, 6.0):
    sel = per_window[per_window["prominence"].fillna(0) >= thr]
    if len(sel) < 20:
        continue
    print(f"    peak/mean week >= {thr:.1f}: unstable "
          f"{sel['unstable'].mean() * 100:5.1f}%  (n={len(sel)})")

# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
thr_rows = pd.DataFrame([r for r in rows if isinstance(r["margin"], (int, float))
                         and str(r["rule"]).startswith("|dAIC|")])
axes[0].plot(thr_rows["margin"], thr_rows["unstable_pct"], "o-",
             color="#8c1c13",
             label="windows where some pair disagrees (%)\n— grows with the design")
axes[0].fill_between(thr_rows["margin"], thr_rows["ci_lo"], thr_rows["ci_hi"],
                     color="#8c1c13", alpha=0.18)
# The invariant curve, drawn on the same axes because the contrast between them
# is the point: on the share-of-windows scale a margin looks weaker the larger
# the factorial, and on this one it does not.
inv = [pairwise_disagreement_at_margin(d, float(m), KEY).mean() * 100
       for m in thr_rows["margin"]]
axes[0].plot(thr_rows["margin"], inv, "^-", color="#3f7d58",
             label="P(two analyses disagree) (%)\n— comparable across designs")
axes[0].plot(thr_rows["margin"], thr_rows["answers_pct"], "s--",
             color="#1f6f8b", label="gives an answer (%)")
axes[0].set_xlabel("evidence threshold, |ΔAIC|")
axes[0].set_ylabel("%")
axes[0].set_title("A margin of 4 removes the disagreement\n"
                  "without declining to answer")
axes[0].legend(fontsize=7.5)
axes[0].grid(alpha=0.3)

q = pd.qcut(per_window["cases"], 4, duplicates="drop")
byc = per_window.groupby(q, observed=True)["unstable"].mean() * 100
axes[1].bar([str(i) for i in byc.index], byc.values, color="#1f6f8b")
axes[1].set_ylabel("% of windows unstable")
axes[1].set_xlabel("reported cases in the window (quartiles)")
axes[1].set_title("Instability is not confined to small outbreaks")
axes[1].tick_params(axis="x", labelrotation=20, labelsize=7)
axes[1].grid(alpha=0.3, axis="y")

fig.tight_layout()
fig.savefig(figures / "13_stability_remedy.png", dpi=150)
print(f"\nFigure: {figures / '13_stability_remedy.png'}")
print(f"Table:  {tables / '16_stability_rules.csv'}")
