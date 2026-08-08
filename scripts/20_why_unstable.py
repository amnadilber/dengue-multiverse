"""
Pipeline step 20 — why the instability exists, and why a small margin removes it.

Step 19 produced two results that needed explaining rather than reporting.
Stricter criteria were *less* stable — BIC by sign and a likelihood-ratio test
at p < 0.01 both worse than AIC by sign — and a margin of 4 applied on the BIC
scale left most windows unstable where the same margin on the AIC scale left
almost none. The figures are printed by the script rather than repeated here,
so that this description cannot go stale when the design changes.

**A hypothesis that turned out to be wrong.** The first explanation tried was
that a threshold flips whenever comparisons sit near it, so instability should
track the density of comparisons at the boundary. It does not. Scanning the
threshold from -14 to +6 gives instability of 95% where the local density is 4%
and 59% where the density is 16% — the relationship runs backwards. The
hypothesis is recorded here because it was tested and rejected, and because the
correct explanation only became visible once it was.

**What is actually happening.** The quantity to measure is not the global density
but the behaviour within each outbreak. Two things are true at once, and they
sound contradictory until separated:

* The choices move Delta AIC enormously — a median within-window range two
  orders of magnitude larger than the median |Delta AIC| being read as evidence.
* Yet the margin needed to eliminate disagreement within a window is small,
  around four units.

Both hold because most of that enormous variation stays on one side of zero. It
changes how strong the evidence looks without changing which model it favours.
Where the direction does flip, it flips marginally: the dissenting combinations
are those whose evidence was weak to begin with.

That is why abstaining below a margin of 4 works, why it is not a property of the
number 4 but of where the dissent lives, and why moving the boundary elsewhere —
which is what a stricter criterion or an off-centre band does — cannot help.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import latest_factorial  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2, spearmanr  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")

# The richest factorial available; see dengue_pk.robustness.FACTORIAL_TABLES.
SRC = latest_factorial(tables)
d = pd.read_csv(SRC)
print(f"source: {SRC.name}")
KEY = ["country", "unit", "window_start"]
counts = d.groupby(KEY).size()
N_COMBOS = int(counts.max())
d = d.set_index(KEY).loc[counts[counts == N_COMBOS].index].reset_index()
print(f"{d.groupby(KEY).ngroups} windows with all {N_COMBOS} combinations")
n_windows = d.groupby(KEY).ngroups
D_K = 2

print(f"{len(d):,} fits, {n_windows} windows\n")

# ---------------------------------------------------------------------------
print("=" * 80)
print("1. HOW MUCH DO THE CHOICES MOVE THE EVIDENCE?")
print("=" * 80)
spread = d.groupby(KEY)["delta_aic"].agg(lambda s: s.max() - s.min())
median_abs = d["delta_aic"].abs().median()
print(f"  within-window range of delta AIC across the {N_COMBOS} combinations:")
print(f"    quartiles {spread.quantile(.25):8.1f} / {spread.median():8.1f} / "
      f"{spread.quantile(.75):8.1f}   (90th pct {spread.quantile(.9):.0f})")
print(f"  median |delta AIC| actually interpreted as evidence: {median_abs:.1f}")
print(f"\n  The analyst's own choices move the evidence by "
      f"{spread.median() / median_abs:.0f}x the size of the difference being read.")

# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("2. BUT THE DISAGREEMENT IS MARGINAL")
print("=" * 80)
rows = []
for key, sub in d.groupby(KEY):
    x = sub["delta_aic"].to_numpy()
    m = 0.0
    while m <= 60:
        v = np.where(x < -m, 1, np.where(x > m, -1, 0))
        s = v[v != 0]
        if len(s) == 0 or len(set(s.tolist())) == 1:
            break
        m += 0.25
    # Dissenters under the sign rule, and separately under a margin of 4: the
    # second is what a reader following the recommendation would actually face,
    # and conflating them overstates how isolated the surviving dissent is.
    speaks4 = np.abs(x) > 4.0
    n4 = (int(min((x[speaks4] < 0).sum(), (x[speaks4] > 0).sum()))
          if speaks4.any() else 0)
    rows.append(dict(margin_needed=m, spread=float(x.max() - x.min()),
                     n_dissent=int(min((x < 0).sum(), (x > 0).sum())),
                     n_dissent_at_4=n4, n_speaking_at_4=int(speaks4.sum())))
need = pd.DataFrame(rows)
print(f"  margin needed to remove all disagreement, per window:")
print(f"    median {need['margin_needed'].median():.2f}, "
      f"90th pct {need['margin_needed'].quantile(.9):.2f}, "
      f"95th pct {need['margin_needed'].quantile(.95):.2f}, "
      f"max {need['margin_needed'].max():.2f}")
rho, p = spearmanr(need["margin_needed"], need["spread"])
print(f"  correlation with that window's spread: rho = {rho:+.3f} (p = {p:.2g})")

# The conclusion is drawn from the measurement rather than asserted alongside
# it. On the five-factor design this correlation was about zero and the text
# here read "the margin needed is unrelated to how far the choices move the
# evidence" — which stopped being true when the design grew, while the sentence
# would have gone on printing.
if abs(rho) < 0.2:
    print("\n  The margin needed is unrelated to how far the choices move the")
    print("  evidence: dissent is not produced by large movements but by")
    print("  comparisons that were weak in the first place.")
else:
    direction = "wider" if rho > 0 else "narrower"
    print(f"\n  The margin needed rises with the spread ({direction} windows need")
    print("  more), so dissent is NOT purely a property of comparisons that were")
    print("  weak to begin with. In a design this large some windows carry")
    print("  disagreement that is confident on both sides, and no band around")
    print("  zero can reach it.")

frac = float((need["margin_needed"] <= 4.0).mean())
print(f"\n  In {frac * 100:.1f}% of outbreaks, every combination that dissents from")
print("  the majority verdict does so on evidence weaker than delta AIC = 4.")
print(f"  In the remaining {100 - frac * 100:.1f}% at least one dissenter is confident.")

# Those two facts are compatible and the pair is the finding: a margin leaves
# *some* window with a dissenter more often than it used to, while the *share*
# of analyses dissenting within a window collapses. Reporting only the first
# reads as "the remedy failed"; only the second, as "the remedy is perfect".
surviving = need.loc[need["margin_needed"] > 4.0]
if len(surviving):
    share = surviving["n_dissent_at_4"] / surviving["n_speaking_at_4"].clip(lower=1)
    print(f"\n  Where dissent survives a margin of 4, how large is the minority?")
    print(f"  median {surviving['n_dissent_at_4'].median():.0f} of "
          f"{surviving['n_speaking_at_4'].median():.0f} analyses still speaking "
          f"({share.median() * 100:.1f}% of them);")
    print(f"  10th percentile {share.quantile(.1) * 100:.1f}%, "
          f"90th percentile {share.quantile(.9) * 100:.1f}%.")

need.to_csv(tables / "18_margin_needed.csv", index=False)

# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("3. THE REJECTED HYPOTHESIS, RECORDED")
print("=" * 80)
x = d["delta_aic"].to_numpy()


def instability_at(threshold: float) -> float:
    unstable = 0
    for _, g in d.groupby(KEY):
        v = (g["delta_aic"].to_numpy() < threshold)
        if v.any() and not v.all():
            unstable += 1
    return unstable / n_windows


grid = np.arange(-14.0, 6.01, 2.0)
prof = pd.DataFrame([dict(threshold=t, unstable=instability_at(float(t)),
                          density=float(np.mean(np.abs(x - t) <= 1.0)))
                     for t in grid])
rho2, p2 = spearmanr(prof["density"], prof["unstable"])
print(f"  instability against local density at the threshold: "
      f"rho = {rho2:+.3f} (p = {p2:.2g})")
print("  Negative, not positive: the density hypothesis is rejected.")
prof.to_csv(tables / "19_threshold_profile.csv", index=False)

named = {
    "LRT p<0.01": 2 * D_K - chi2.ppf(0.99, D_K),
    "BIC (38 wks)": 2 * D_K - D_K * np.log(38),
    "LRT p<0.05": 2 * D_K - chi2.ppf(0.95, D_K),
    "AIC sign": 0.0,
}
print("\n  Conventional criteria, ordered by strictness:")
for name, thr in sorted(named.items(), key=lambda kv: kv[1]):
    print(f"    {name:14s} at delta AIC {thr:+6.2f}: "
          f"instability {instability_at(thr) * 100:5.1f}%")
print("\n  Strictness makes it worse. A stricter threshold is still a threshold,")
print("  and it moves the boundary into a region the choices also straddle.")

# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

axes[0].hist(np.log10(spread.clip(lower=1e-3)), bins=40, color="#8c1c13")
axes[0].axvline(np.log10(median_abs), color="k", lw=2,
                label=f"median |ΔAIC| interpreted = {median_abs:.0f}")
axes[0].set_xlabel("log₁₀ range of ΔAIC within one outbreak")
axes[0].set_ylabel("outbreaks")
axes[0].set_title("The choices move the evidence\nby orders of magnitude")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].hist(need["margin_needed"], bins=np.arange(0, 12.5, 0.5),
             color="#3f7d58")
axes[1].axvline(4, color="k", lw=2, ls="--", label="margin of 4")
axes[1].set_xlabel("margin needed to remove all disagreement")
axes[1].set_ylabel("outbreaks")
axes[1].set_title(f"…but dissent is marginal\n{frac * 100:.0f}% of outbreaks "
                  "need less than 4")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

axes[2].scatter(need["spread"].clip(upper=1e4), need["margin_needed"], s=28,
                alpha=0.6, color="#1f6f8b")
axes[2].set_xscale("log")
axes[2].set_xlabel("range of ΔAIC within the outbreak")
axes[2].set_ylabel("margin needed")
axes[2].set_title(f"and unrelated to how far\nthe evidence moved (ρ = {rho:+.2f})")
axes[2].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(figures / "15_why_unstable.png", dpi=150)
print(f"\nFigure: {figures / '15_why_unstable.png'}")
print(f"Tables: {tables / '18_margin_needed.csv'}, "
      f"{tables / '19_threshold_profile.csv'}")
