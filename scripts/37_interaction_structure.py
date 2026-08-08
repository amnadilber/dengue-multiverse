"""Pipeline step 37 — is the outbreak-by-analysis interaction predictable?

Step 35 established that three-fifths of the variation in the verdict is an
interaction: which analysis endorses depends on which outbreak. That is a
negative result as it stands, and the obvious referee question is whether it is
merely unexplained or genuinely irreducible. The two have opposite consequences.
If the interaction can be predicted from observable properties of an outbreak,
then a practitioner can anticipate which choices will matter for their data,
which is a tool. If it cannot, the case for reporting the spread rather than a
verdict is much stronger, because there is nothing else to report.

One hypothesis has a direction fixed in advance rather than fitted afterwards.
The observation model is the dominant factor, and the reason a Poisson likelihood
misbehaves is that it denies overdispersion. Its effect should therefore be
largest in the outbreaks whose counts are most overdispersed --- and dispersion
is measured in step 34 by a model-free index that no part of the factorial uses.
If the prediction holds, the interaction is not arbitrary: it is largest exactly
where the assumption it rests on is most violated.

The other five factors have no comparable prior, so they are tested against the
same outbreak descriptors and reported with false-discovery-rate control, as in
step 31. A study whose subject is undisclosed analytical freedom does not get to
run thirty correlations and quote the significant ones.

Per outbreak, for each factor, the effect is the difference in the share of
analyses endorsing climate forcing between two levels of that factor, averaged
over every combination of the others. That is exactly the outbreak's own slice of
the interaction.

Writes ``46_interaction_structure.csv``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (benjamini_hochberg,  # noqa: E402
                                  complete_windows, latest_factorial)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")
KEY = ["country", "unit", "window_start"]

#: Each factor reduced to one contrast, so that every outbreak yields one number
#: per factor. Where a factor has three levels the extremes are used.
CONTRASTS = {
    "observation": ("nb", "poisson"),
    "temp_form": ("briere", "loglinear"),
    "rain_lag": (3, 7),
    "train_frac": (0.75, 1.0),
    "structure": ("hostvector", "seir"),
    "params": ("central", "alt"),
}


def effect_per_outbreak(d: pd.DataFrame, factor: str, lo, hi) -> pd.Series:
    """How much moving this factor moves the endorsement rate, per outbreak."""
    a = d[d[factor] == lo].groupby(KEY)["climate_wins"].mean()
    b = d[d[factor] == hi].groupby(KEY)["climate_wins"].mean()
    return (b - a).dropna()


def main() -> None:
    d, n_combos = complete_windows(pd.read_csv(latest_factorial(tables)), KEY)
    print(f"{len(d):,} fits, {d.groupby(KEY).ngroups} outbreaks, "
          f"{n_combos} combinations\n")

    eff = pd.DataFrame({f: effect_per_outbreak(d, f, lo, hi)
                        for f, (lo, hi) in CONTRASTS.items()
                        if f in d.columns})

    # Outbreak descriptors. The dispersion index is model-free and computed by a
    # step that never touches the factorial, so it cannot be circular.
    desc = (d.groupby(KEY)
            .agg(weeks=("weeks", "first"), cases=("cases", "first"))
            .assign(log_cases=lambda x: np.log10(x["cases"])))
    disp = pd.read_csv(tables / "41_dispersion.csv").set_index(KEY)
    desc = desc.join(disp[["dispersion_index", "mean_weekly"]], how="left")
    desc["log_dispersion"] = np.log10(desc["dispersion_index"].clip(lower=1e-3))
    desc["log_mean_weekly"] = np.log10(desc["mean_weekly"].clip(lower=1e-3))

    joined = eff.join(desc, how="inner").dropna(
        subset=["log_dispersion", "log_cases", "weeks"])
    print(f"{len(joined)} outbreaks have both an effect and a dispersion index\n")

    print("=" * 78)
    print("THE PRE-REGISTERED ONE: DOES POISSON HURT MOST WHERE DISPERSION IS WORST?")
    print("=" * 78)
    r, p = stats.spearmanr(joined["log_dispersion"], joined["observation"])
    print(f"  observation-model effect against log dispersion: "
          f"rho = {r:+.3f}  (p = {p:.2g})")
    lo_half = joined[joined["log_dispersion"] <= joined["log_dispersion"].median()]
    hi_half = joined[joined["log_dispersion"] > joined["log_dispersion"].median()]
    print(f"  below median dispersion: switching to Poisson moves the endorsement "
          f"rate by {lo_half['observation'].mean():+.3f}")
    print(f"  above median dispersion: {hi_half['observation'].mean():+.3f}")
    print(f"  ratio {hi_half['observation'].mean() / max(lo_half['observation'].mean(), 1e-9):.2f}x")

    print("\n" + "=" * 78)
    print("EVERY FACTOR AGAINST EVERY DESCRIPTOR, WITH FDR CONTROL")
    print("=" * 78)
    descriptors = ["log_dispersion", "log_cases", "weeks", "log_mean_weekly"]
    rows = []
    for f in eff.columns:
        for c in descriptors:
            sub = joined[[f, c]].dropna()
            if len(sub) < 20:
                continue
            rho, pv = stats.spearmanr(sub[c], sub[f])
            rows.append(dict(factor=f, descriptor=c, n=len(sub),
                             rho=round(float(rho), 3), p=float(pv)))
    res = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    m = len(res)
    res["q"] = np.round(benjamini_hochberg(res["p"].to_numpy()), 4)
    res["p"] = res["p"].round(5)
    print(f"  {m} correlations; about {m * 0.05:.1f} would reach p < 0.05 by chance.")
    print(f"  {(res['p'] < 0.05).sum()} do; {(res['q'] < 0.05).sum()} survive "
          f"false-discovery-rate control.\n")
    print(res.head(10).to_string(index=False))
    res.to_csv(tables / "46_interaction_structure.csv", index=False)

    print("\n" + "=" * 78)
    print("HOW MUCH OF THE INTERACTION DOES ANY OF THIS EXPLAIN?")
    print("=" * 78)
    print("  The honest summary is not a p-value but a share. For each factor,")
    print("  the fraction of the between-outbreak variance in its effect that")
    print("  the four descriptors account for, by ordinary least squares.\n")
    X = np.column_stack([np.ones(len(joined))]
                        + [joined[c].to_numpy(float) for c in descriptors])
    shares = []
    for f in eff.columns:
        y = joined[f].to_numpy(float)
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        r2 = 1.0 - resid.var(ddof=0) / max(y.var(ddof=0), 1e-12)
        shares.append(dict(factor=f, r2=round(float(r2), 3),
                           sd_effect=round(float(y.std(ddof=0)), 3)))
        print(f"  {f:12s} R^2 = {r2:5.3f}   sd of the effect across outbreaks "
              f"{y.std(ddof=0):.3f}")
    sh = pd.DataFrame(shares)
    print(f"\n  Median R^2 across the six factors: {sh['r2'].median():.3f}.")
    print("  Even the factor with a mechanism behind it leaves most of its")
    print("  outbreak-to-outbreak variation unexplained by anything an analyst")
    print("  can see before fitting.")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    axes[0].scatter(joined["log_dispersion"], joined["observation"], s=16,
                    alpha=0.6, color="#1f6f8b")
    axes[0].set_xlabel("log$_{10}$ index of dispersion (model-free)")
    axes[0].set_ylabel("effect of switching to Poisson\n(change in endorsement rate)")
    axes[0].set_title(f"the one prediction with a direction (rho = {r:+.2f})")
    axes[0].grid(alpha=0.3)
    order = sh.sort_values("r2")
    axes[1].barh(order["factor"], order["r2"], color="#8c1c13")
    axes[1].set_xlabel("share of the effect's variation explained by outbreak size,\n"
                       "length, mean count and dispersion")
    axes[1].set_xlim(0, 1)
    axes[1].set_title("how predictable each factor's effect is")
    axes[1].grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(figures / "20_interaction_structure.png", dpi=150)

    print(f"\nFigure: {figures / '20_interaction_structure.png'}")
    print(f"Table:  {tables / '46_interaction_structure.csv'}")


if __name__ == "__main__":
    main()
