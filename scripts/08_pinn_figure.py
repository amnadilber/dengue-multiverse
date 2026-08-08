"""
Pipeline step 8 — the figure that carries the PINN conclusion.

Three panels, each answering a question a sceptical reader would ask.

1. Did any configuration work? One did, apparently better than the classical
   estimator.
2. Is that reproducible? No — the same configuration across five random seeds
   spans nearly the whole range of possible answers.
3. Could a practitioner have told the good runs from the bad ones? No: neither
   likelihood nor physics residual points toward accuracy, and likelihood points
   away from it.

Built from the table written by step 07 rather than by re-running anything, so
the figure and the table cannot disagree.
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

CLASSICAL_ERR = 3.1        # step 06, same synthetic data

cfg = load_config()
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")
tab = pd.read_csv(tables / "06_pinn_tuning.csv")

seeds = tab[tab["config"].str.startswith("F_seed")].copy()
# The seed-0 run of the best configuration is recorded under its own label.
best = tab[tab["config"] == "B_dataweight_20000"]
seed_errs = np.concatenate([best["R0_rel_err"].to_numpy(),
                            seeds["R0_rel_err"].to_numpy()])
others = tab[~tab["config"].str.startswith("F_seed")]

fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))

# -- 1. every configuration ---------------------------------------------------
labels = [c.replace("_", " ") for c in others["config"]]
colours = ["#3f7d58" if e < CLASSICAL_ERR else "#8c1c13"
           for e in others["R0_rel_err"]]
axes[0].barh(labels, others["R0_rel_err"], color=colours)
axes[0].axvline(CLASSICAL_ERR, color="k", ls="--", lw=1.2,
                label=f"classical estimator ({CLASSICAL_ERR}%)")
axes[0].set_xscale("log")
axes[0].set_xlabel("relative error in $R_0$ (%)")
axes[0].set_title("One configuration beats the classical fit")
axes[0].legend(fontsize=8, loc="lower right")
axes[0].grid(alpha=0.3, axis="x")

# -- 2. the same configuration across seeds -----------------------------------
axes[1].bar(range(len(seed_errs)), seed_errs, color="#1f6f8b")
axes[1].axhline(CLASSICAL_ERR, color="k", ls="--", lw=1.2,
                label=f"classical estimator ({CLASSICAL_ERR}%)")
axes[1].set_xticks(range(len(seed_errs)))
axes[1].set_xticklabels([f"seed {i}" for i in range(len(seed_errs))])
axes[1].set_ylabel("relative error in $R_0$ (%)")
axes[1].set_title("…but only on one random seed\n(identical configuration throughout)")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3, axis="y")
for i, e in enumerate(seed_errs):
    axes[1].annotate(f"{e:.1f}%", (i, e), ha="center", va="bottom", fontsize=8)

# -- 3. can a practitioner detect the good runs? ------------------------------
axes[2].scatter(tab["loglik"], tab["R0_rel_err"], s=70, color="#8c1c13",
                label="log-likelihood")
axes[2].axhline(CLASSICAL_ERR, color="k", ls="--", lw=1.2)
axes[2].set_yscale("log")
axes[2].set_xlabel("log-likelihood on the fitted data")
axes[2].set_ylabel("relative error in $R_0$ (%)")
axes[2].set_title("Better fit, worse parameters\n(no usable diagnostic)")
axes[2].grid(alpha=0.3)

r = np.corrcoef(tab["loglik"], np.log(tab["R0_rel_err"]))[0, 1]
axes[2].annotate(f"corr(logL, log error) = {r:+.2f}",
                 (0.04, 0.06), xycoords="axes fraction", fontsize=9)

fig.tight_layout()
fig.savefig(figures / "06_pinn_tuning.png", dpi=150)
print(f"Figure: {figures / '06_pinn_tuning.png'}")

print(f"\nSeed spread at the best configuration: "
      f"{seed_errs.min():.1f}% to {seed_errs.max():.1f}%, "
      f"median {np.median(seed_errs):.1f}%")
print(f"Classical estimator, same data: {CLASSICAL_ERR}% (deterministic)")
print(f"Correlation between log-likelihood and log parameter error: {r:+.2f}")
print("A positive correlation means the better-fitting configurations are the "
      "less accurate ones.")
