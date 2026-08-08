"""
Pipeline step 19 — is the instability specific to AIC?

The obvious objection to step 17 is that it tested one criterion. If AIC is
peculiarly fragile, the finding says something about AIC rather than about the
evidence, and the remedy in step 18 would be a fix for the wrong problem.

Three criteria are compared on the same fits, and all three are computable
from what was already stored, so nothing is refitted:

* **AIC**, the baseline.
* **BIC**, which penalises the extra parameters by log n rather than 2 and so is
  stricter for any window longer than about seven weeks.
* **The likelihood-ratio test.** The two models are nested — the climate model
  reduces to the constant model when both climate coefficients vanish — so
  2(logL_climate - logL_constant) is asymptotically chi-square on 2 degrees of
  freedom, and the conventional 5% test is the significance test an
  epidemiologist would actually run.

The arithmetic that makes this free: AIC = -2 logL + 2k, so logL = k - AIC/2, and
with k known per fit the other two criteria follow. Delta k is 2 in every
comparison, since the climate model adds exactly the two climate coefficients.

The question is not which criterion is best. It is whether any of them, used
conventionally, gives a verdict that does not move with the analyst.
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
from scipy.stats import chi2  # noqa: E402

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
rng = np.random.default_rng(cfg["seed"])

# Parameter counts. The climate model estimates beta_0, a_temp, a_rain, pop_frac
# and i0_frac; the constant model drops the two climate coefficients. The
# negative binomial adds its dispersion, estimated from the same data.
d["k_climate"] = 5 + (d["observation"] == "nb").astype(int)
d["k_constant"] = 3 + (d["observation"] == "nb").astype(int)
D_K = 2

# delta AIC = AIC_c - AIC_n = -2 dlogL + 2 dk  =>  dlogL = (2 dk - dAIC) / 2
d["delta_loglik"] = (2 * D_K - d["delta_aic"]) / 2.0
d["delta_bic"] = -2 * d["delta_loglik"] + D_K * np.log(d["weeks"])
d["lrt_stat"] = 2 * d["delta_loglik"]
d["lrt_p"] = chi2.sf(np.maximum(d["lrt_stat"], 0.0), D_K)

print(f"{len(d):,} fits, {n_windows} windows\n")
print("Sanity: delta BIC - delta AIC should equal 2(log n - 2).")
chk = d["delta_bic"] - d["delta_aic"]
expect = D_K * (np.log(d["weeks"]) - 2)
print(f"  max discrepancy {np.max(np.abs(chk - expect)):.2e}\n")

CRITERIA = {
    "AIC, sign": lambda g: np.sign(-g["delta_aic"]).astype(int),
    "BIC, sign": lambda g: np.sign(-g["delta_bic"]).astype(int),
    "LRT at p < 0.05": lambda g: np.where(
        (g["lrt_stat"] > 0) & (g["lrt_p"] < 0.05), 1, -1),
    "LRT at p < 0.01": lambda g: np.where(
        (g["lrt_stat"] > 0) & (g["lrt_p"] < 0.01), 1, -1),
    "AIC, margin 4": lambda g: np.where(
        g["delta_aic"] < -4, 1, np.where(g["delta_aic"] > 4, -1, 0)),
    "BIC, margin 4": lambda g: np.where(
        g["delta_bic"] < -4, 1, np.where(g["delta_bic"] > 4, -1, 0)),
}


def assess(rule):
    unstable = decided = climate = 0
    for _, g in d.groupby(KEY):
        v = np.asarray(rule(g))
        spoken = v[v != 0]
        if len(spoken) == 0:
            continue
        decided += 1
        if len(set(spoken.tolist())) > 1:
            unstable += 1
        elif spoken[0] == 1:
            climate += 1
    return (unstable / decided if decided else np.nan,
            decided / n_windows,
            climate / decided if decided else np.nan)


def boot(rule, n=200):
    groups = [g for _, g in d.groupby(KEY)]
    out = []
    for _ in range(n):
        pick = rng.integers(0, len(groups), len(groups))
        unstable = decided = 0
        for i in pick:
            v = np.asarray(rule(groups[i]))
            spoken = v[v != 0]
            if len(spoken) == 0:
                continue
            decided += 1
            if len(set(spoken.tolist())) > 1:
                unstable += 1
        if decided:
            out.append(unstable / decided)
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


print("=" * 82)
print("INSTABILITY UNDER EACH CRITERION, USED CONVENTIONALLY")
print("=" * 82)
rows = []
for name, rule in CRITERIA.items():
    u, share, clim = assess(rule)
    lo, hi = boot(rule)
    print(f"  {name:18s} unstable {u * 100:5.1f}% [{lo * 100:4.1f}-{hi * 100:4.1f}]   "
          f"answers {share * 100:5.1f}%   "
          f"climate endorsed in {clim * 100:5.1f}% of decided windows")
    rows.append(dict(criterion=name, unstable_pct=round(u * 100, 1),
                     ci_lo=round(lo * 100, 1), ci_hi=round(hi * 100, 1),
                     answers_pct=round(share * 100, 1),
                     climate_pct=round(clim * 100, 1)))

pd.DataFrame(rows).to_csv(tables / "17_criteria_comparison.csv", index=False)

print("\n" + "=" * 82)
print("WHY: WHERE THE THREE CRITERIA PUT THEIR DECISION BOUNDARY")
print("=" * 82)
med_n = float(np.median(d["weeks"]))
lrt5 = chi2.ppf(0.95, D_K)
lrt1 = chi2.ppf(0.99, D_K)
print(f"  median window length: {med_n:.0f} weeks")
print(f"  AIC sign          <=> delta AIC < {0.0:+.2f}")
print(f"  LRT p<0.05        <=> delta AIC < {2 * D_K - lrt5:+.2f}")
print(f"  LRT p<0.01        <=> delta AIC < {2 * D_K - lrt1:+.2f}")
print(f"  BIC sign          <=> delta AIC < {2 * D_K - D_K * np.log(med_n):+.2f} "
      f"(at the median length)")
print("\n  The conventional significance test sits at a margin of about 2, which")
print("  step 18 showed still leaves most windows unstable. BIC lands near the")
print("  margin that works — not by design, but because log n happens to be the")
print("  right order of penalty for series of this length.")

# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8))

names = [r["criterion"] for r in rows]
vals = [r["unstable_pct"] for r in rows]
err = np.array([[r["unstable_pct"] - r["ci_lo"] for r in rows],
                [r["ci_hi"] - r["unstable_pct"] for r in rows]])
colours = ["#8c1c13" if v > 40 else "#3f7d58" for v in vals]
axes[0].barh(names, vals, xerr=err, color=colours, capsize=3)
axes[0].set_xlabel("% of outbreaks whose verdict changes with the analysis choices")
axes[0].set_title("No conventional criterion is stable;\na margin is what fixes it")
axes[0].grid(alpha=0.3, axis="x")

axes[1].hist(d["delta_aic"].clip(-40, 40), bins=80, color="#c9ccd1")
for x, lab, col in ((0, "AIC sign", "#8c1c13"),
                    (-(lrt5 - 2 * D_K), "LRT p<0.05", "#e08a1e"),
                    (-4, "margin 4", "#3f7d58")):
    axes[1].axvline(x, color=col, lw=1.8, label=lab)
axes[1].set_xlabel("ΔAIC (climate − constant), clipped to ±40")
axes[1].set_ylabel("fits")
axes[1].set_title("Most comparisons sit where no criterion\nshould be declaring a winner")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(figures / "14_criteria_comparison.png", dpi=150)
print(f"\nFigure: {figures / '14_criteria_comparison.png'}")
print(f"Table:  {tables / '17_criteria_comparison.csv'}")
