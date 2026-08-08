"""
Pipeline step 12 — pooling every PINN run to ask whether any diagnostic works.

Twenty-one PINN fits now exist across steps 07 and 10, all on the same synthetic
data with the same known truth. That is enough to answer the question a
practitioner actually faces: with no ground truth available, is there anything
computable from a finished run that indicates whether to trust it?

Three candidates are tested — the data likelihood, the physics residual, and the
estimated population at risk relative to the census total. The first two are what
a PINN paper conventionally reports. The third is not a diagnostic so much as a
sanity check, and it is included because it turned out to be the only informative
one.

The step-10 runs added a pattern worth isolating: the configurations that drove
the physics residual lowest returned population estimates covering most of the
country, while those that left the residual high returned populations closer to
the truth. If that holds across all runs it is a strong statement, because it
means the diagnostic most closely associated with "the PINN is working" is
actively anti-correlated with the parameters being right.
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

POP_FRAC_TRUTH = 0.0025
R0_TRUTH = 1.7
CLASSICAL_ERR = 3.1

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")

frames = []
for fname, source in (("06_pinn_tuning.csv", "step 07"),
                      ("08_pinn_timemarching.csv", "step 10")):
    path = tables / fname
    if path.exists():
        d = pd.read_csv(path)
        d["source"] = source
        frames.append(d)

if not frames:
    raise SystemExit("no PINN result tables found; run steps 07 and 10 first")

tab = pd.concat(frames, ignore_index=True)
tab["pop_frac_ratio"] = tab["pop_frac"] / POP_FRAC_TRUTH
tab["pop_frac_log_err"] = np.abs(np.log10(tab["pop_frac_ratio"]))
tab = tab[np.isfinite(tab["R0_rel_err"]) & np.isfinite(tab["physics_loss"])]

print(f"{len(tab)} PINN runs pooled "
      f"({', '.join(f'{n}: {c}' for n, c in tab['source'].value_counts().items())})\n")

print(f"R0 error:      min {tab['R0_rel_err'].min():5.1f}%  "
      f"median {tab['R0_rel_err'].median():5.1f}%  "
      f"max {tab['R0_rel_err'].max():5.1f}%   (classical {CLASSICAL_ERR}%)")
print(f"runs beating the classical estimator: "
      f"{int((tab['R0_rel_err'] < CLASSICAL_ERR).sum())} of {len(tab)}")
print(f"population at risk: estimated between "
      f"{tab['pop_frac_ratio'].min():.2f}x and {tab['pop_frac_ratio'].max():.0f}x "
      f"the truth\n")


def report(x_name, y_name, label):
    x, y = tab[x_name].to_numpy(float), tab[y_name].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y) & (y > 0)
    r = np.corrcoef(x[ok], np.log(y[ok]))[0, 1]
    direction = ("points the WRONG way" if r > 0.15 else
                 "points the right way" if r < -0.15 else "is uninformative")
    print(f"  {label:38s} corr = {r:+.2f}   {direction}")
    return r


print("Does any diagnostic identify the accurate runs?")
print("  (correlation with log R0 error; negative would mean a better "
      "diagnostic value goes with a better estimate)")
r_ll = report("loglik", "R0_rel_err", "log-likelihood on the data")
r_ph = report("physics_loss", "R0_rel_err", "physics residual (lower = better)")
r_pop = report("pop_frac_log_err", "R0_rel_err", "implausibility of the population")

print("\nAnd the pattern step 10 suggested:")
x = np.log10(tab["physics_loss"].to_numpy(float))
y = tab["pop_frac_log_err"].to_numpy(float)
ok = np.isfinite(x) & np.isfinite(y)
r_cross = np.corrcoef(x[ok], y[ok])[0, 1]
print(f"  corr(log physics residual, error in population) = {r_cross:+.2f}")
if r_cross < -0.3:
    print("  Negative: the runs that satisfy the equations best are the ones")
    print("  that place most of the country in the transmission catchment.")

low = tab[tab["physics_loss"] < 10]
high = tab[tab["physics_loss"] >= 10]
if len(low) and len(high):
    print(f"\n  physics residual < 10  ({len(low):2d} runs): population "
          f"{low['pop_frac_ratio'].median():8.1f}x truth, "
          f"R0 error {low['R0_rel_err'].median():5.1f}%")
    print(f"  physics residual >= 10 ({len(high):2d} runs): population "
          f"{high['pop_frac_ratio'].median():8.1f}x truth, "
          f"R0 error {high['R0_rel_err'].median():5.1f}%")

tab.to_csv(tables / "11_pinn_pooled.csv", index=False)

# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
colours = {"step 07": "#8c1c13", "step 10": "#1f6f8b"}

for ax, (xcol, xlabel, logx, corr) in zip(axes, [
        ("loglik", "log-likelihood on the data", False, r_ll),
        ("physics_loss", "physics residual", True, r_ph),
        ("pop_frac_ratio", "estimated population / truth", True, r_pop)]):
    for source, grp in tab.groupby("source"):
        ax.scatter(grp[xcol], grp["R0_rel_err"], s=60, alpha=0.8,
                   color=colours.get(source, "grey"), label=source)
    ax.axhline(CLASSICAL_ERR, color="k", ls="--", lw=1.2,
               label=f"classical ({CLASSICAL_ERR}%)")
    if logx:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("relative error in $R_0$ (%)")
    ax.set_title(f"corr with log error: {corr:+.2f}")
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=8)
fig.suptitle("Twenty-one PINN runs: is any finished-run diagnostic informative?",
             y=1.02)
fig.tight_layout()
fig.savefig(figures / "10_pinn_diagnostics.png", dpi=150)
print(f"\nFigure: {figures / '10_pinn_diagnostics.png'}")
print(f"Table:  {tables / '11_pinn_pooled.csv'}")
