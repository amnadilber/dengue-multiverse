"""Pipeline step 38 — the same three-way partition where the answer is known.

Step 35 partitions the verdict into an outbreak effect, an analysis effect and
their interaction, and finds the interaction largest. On real data that is a
description with nothing to compare it against. Two comparisons are available and
neither had been made.

**The simulation.** Both arms were generated from each window's own fitted
parameters and run through the identical factorial, one with no climate effect
and one with an effect of realistic size. The partition can be computed on each.
What it should show, if it measures anything: with no signal to detect, the
result is whatever the convention gives, the same way in every window --- a large
**analysis** main effect. With a real signal, the signal meets each epidemic
differently and the **interaction** should dominate. If the real data's
fingerprint matches one arm and not the other, that is evidence of a kind the
one-dimensional statistics could not provide, and it is why the paper's earlier
attempt at this comparison had to be withdrawn: a single number could not
discriminate between two arms that differ in three.

**The reproduction number.** The paper reports a one-way split for $R_0$ ---
mostly between outbreaks under a constant model, mostly not under a climate
model. That split has the same defect the verdict's did: the remainder is not
the analysis. Whether a convention could stabilise $R_0$ depends entirely on
which part of the remainder it is.

Reads stored tables; refits nothing. Writes ``47_decomposition_under_truth.csv``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (complete_windows, latest_factorial,  # noqa: E402
                                  two_way_decomposition)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")
KEY = ["country", "unit", "window_start"]
FACTORS = ["observation", "temp_form", "rain_lag", "train_frac", "structure",
           "params"]


def balanced(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    """Drop non-finite values, then keep only outbreaks still fully crossed.

    Without the second step the partition no longer sums to one and every share
    reported from it is meaningless — quietly, since nothing raises.
    """
    d = frame.assign(_v=pd.to_numeric(frame[column], errors="coerce")
                     .replace([np.inf, -np.inf], np.nan)).dropna(subset=["_v"])
    if d.empty:
        return d
    counts = d.groupby(KEY).size()
    return d.set_index(KEY).loc[counts[counts == counts.max()].index].reset_index()


def main() -> None:
    rows = []

    print("=" * 78)
    print("1. THE SAME PARTITION UNDER A KNOWN TRUTH")
    print("=" * 78)
    sim = pd.read_csv(tables / "28_false_positive_6factor.csv")
    print("  Both arms: same windows, same 144 analyses, same fitted parameters.")
    print("  Only the truth differs.\n")
    print(f"  {'data':28s} {'outbreak':>9} {'analysis':>9} {'interaction':>12}")
    for truth in ("no_effect", "climate"):
        g = sim[sim["truth"] == truth]
        cells = [f for f in FACTORS if f in g.columns]
        sub, n = complete_windows(g, KEY)
        r = two_way_decomposition(sub, cells, KEY)
        label = ("simulated, no climate effect" if truth == "no_effect"
                 else "simulated, real effect")
        print(f"  {label:28s} {r['outbreak'] * 100:8.1f}% {r['analysis'] * 100:8.1f}%"
              f" {r['interaction'] * 100:11.1f}%")
        rows.append(dict(source=label, n_windows=sub.groupby(KEY).ngroups,
                         **{k: round(v, 4) for k, v in r.items()}))

    d, n_combos = complete_windows(pd.read_csv(latest_factorial(tables)), KEY)
    cells = [f for f in FACTORS if f in d.columns]
    real = two_way_decomposition(d, cells, KEY)
    print(f"  {'real outbreaks':28s} {real['outbreak'] * 100:8.1f}% "
          f"{real['analysis'] * 100:8.1f}% {real['interaction'] * 100:11.1f}%")
    rows.append(dict(source="real outbreaks", n_windows=d.groupby(KEY).ngroups,
                     **{k: round(v, 4) for k, v in real.items()}))

    print("\n  The analysis main effect is the discriminating component. With")
    print("  nothing to detect it is the largest term, because the answer is")
    print("  then whatever the convention gives and the convention does not")
    print("  vary by window. The real data does not look like that.")

    print("\n" + "=" * 78)
    print("2. THE NUMBER THE FIELD REPORTS, PARTITIONED THE SAME WAY")
    print("=" * 78)
    print("  A one-way split cannot say whether a convention would help. This can:")
    print("  only the analysis main effect is reachable by agreeing on a method.\n")
    print(f"  {'quantity':34s} {'outbreak':>9} {'analysis':>9} {'interaction':>12}")
    for col, label in (("R0_constant", "R0, constant model"),
                       ("R0_climate", "R0, climate model")):
        if col not in d.columns:
            continue
        for scale, fn in (("raw", lambda s: s), ("log scale", np.log)):
            sub = balanced(d, col)
            sub = sub[sub["_v"] > 0] if scale == "log scale" else sub
            if sub.empty:
                continue
            counts = sub.groupby(KEY).size()
            sub = sub.set_index(KEY).loc[
                counts[counts == counts.max()].index].reset_index()
            sub["_v"] = fn(sub["_v"])
            r = two_way_decomposition(sub, cells, KEY, column="_v")
            print(f"  {label + ', ' + scale:34s} {r['outbreak'] * 100:8.1f}% "
                  f"{r['analysis'] * 100:8.1f}% {r['interaction'] * 100:11.1f}%")
            rows.append(dict(source=f"{label}, {scale}",
                             n_windows=sub.groupby(KEY).ngroups,
                             **{k: round(v, 4) for k, v in r.items()}))

    print("\n  Under a constant model R0 is mostly a property of the epidemic and")
    print("  the analysis contributes about a twentieth. Adding climate terms does")
    print("  not transfer that to the analysis: it transfers it to the interaction.")
    print("  There is no convention to agree on that would fix it.")

    out = pd.DataFrame(rows)
    out.to_csv(tables / "47_decomposition_under_truth.csv", index=False)

    fig, ax = plt.subplots(figsize=(9, 4.2))
    plot = out[out["source"].str.contains("simulated|real")].iloc[::-1]
    left = np.zeros(len(plot))
    for term, colour in (("outbreak", "#1f6f8b"), ("analysis", "#8c1c13"),
                         ("interaction", "#c9a227")):
        ax.barh(plot["source"], plot[term] * 100, left=left, color=colour,
                label=term)
        left = left + plot[term].to_numpy() * 100
    ax.set_xlabel("share of the variation in the verdict (%)")
    ax.set_xlim(0, 100)
    ax.legend(loc="lower right")
    ax.set_title("The fingerprint of a real signal, and of none")
    fig.tight_layout()
    fig.savefig(figures / "21_decomposition_under_truth.png", dpi=150)

    print(f"\nFigure: {figures / '21_decomposition_under_truth.png'}")
    print(f"Table:  {tables / '47_decomposition_under_truth.csv'}")


if __name__ == "__main__":
    main()
