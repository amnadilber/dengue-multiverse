"""
Pipeline step 10 — the remaining PINN remedy: progressive time-window training.

Step 07's docstring listed sequential time-marching among the configurations to
try and then did not implement it. This closes that gap, and the omission is worth
being precise about, because the standard remedy does not transfer unchanged.

Krishnapriyan et al. (2021) recommend sequence-to-sequence training for PINNs on
problems with sharp transients: partition time, solve each window with its own
network, and hard-constrain each to begin at the previous window's endpoint. That
prescription is for *forward* problems, where the coefficients are known. Here the
coefficients are what is being estimated, and they are global — a separate network
per window with its own free parameters would not be estimating one epidemic's
transmission rate, it would be estimating a different one per window.

The transferable part is the curriculum: begin on a short interval where the
dynamics are gentle, and extend it, carrying the network and the shared parameters
forward. Collocation points and observations are both restricted to the current
window, so the network is never asked to represent behaviour it has no data for.

Two variants are tested against the best configuration from step 07:

* **expanding** — the window grows from 20% to 100% of the horizon in stages.
* **expanding + restart** — the same, but the optimiser state is rebuilt at each
  stage, so momentum accumulated on a shorter window cannot carry a stale
  direction into a longer one.

Both are run across the same five seeds used in step 07, because a single seed
tells us nothing: that was the lesson of step 07.
"""

from __future__ import annotations

import copy
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.inference import (Dataset, estimate_dispersion_k,  # noqa: E402
                                 nb_loglik, predict)
from dengue_pk.models import FixedParams, basic_reproduction_number  # noqa: E402
from dengue_pk.pinn import InversePINN  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tensorflow as tf  # noqa: E402

BASE = load_config()
fixed = FixedParams.from_config(BASE)
processed = resolve(BASE, "processed")
figures = resolve(BASE, "figures")
tables = resolve(BASE, "tables")

CLASSICAL_ERR = 3.1
SEEDS = (0, 1, 2, 3, 4)

TRUTH = dict(beta_0=0.17, a_temp=0.2, a_rain=0.1, pop_frac=0.0025, i0_frac=5e-6)
R0_TRUTH = basic_reproduction_number(TRUTH["beta_0"], fixed)

ref = pd.read_csv(processed / "national_2013.csv", parse_dates=["week_start"])
template = Dataset.from_frame(ref, label="synthetic")
mu_true = predict(TRUTH, template, fixed, "climate", BASE)
rng = np.random.default_rng(BASE["seed"])
lam = rng.gamma(shape=5.0, scale=np.maximum(mu_true, 1e-12) / 5.0)
DATA = Dataset(template.days, rng.poisson(lam).astype(float), template.temp_c,
               template.z_rain, template.population, "synthetic")

print(f"synthetic truth R0 = {R0_TRUTH:.3f}; classical estimator "
      f"{CLASSICAL_ERR}% error on this data\n")


class CurriculumPINN(InversePINN):
    """Train on a growing time window, carrying network and parameters forward."""

    def train_curriculum(self, n_stages: int = 8, restart_optimiser: bool = False):
        pc = self.cfg["inference"]["pinn"]
        epochs_per_stage = max(self.epochs // n_stages, 1)
        n_weeks = len(self.data.days)

        full_col = self.tau_col_tf
        full_b, full_zr = self.b_col, self.zr_col
        full_edges, full_cases = self.tau_edges, self.cases

        t0 = time.time()
        for stage in range(1, n_stages + 1):
            frac = 0.2 + 0.8 * stage / n_stages
            # Restrict collocation points to the current window.
            keep = full_col.numpy().ravel() <= frac
            keep[0] = True
            self.tau_col_tf = tf.constant(full_col.numpy()[keep], dtype=tf.float32)
            self.b_col = tf.constant(full_b.numpy()[keep], dtype=tf.float32)
            self.zr_col = tf.constant(full_zr.numpy()[keep], dtype=tf.float32)

            # Restrict the observations likewise: the network must not be scored
            # on weeks whose dynamics it has not yet been asked to represent.
            n_obs = max(int(round(n_weeks * frac)), 3)
            self.tau_edges = tf.constant(full_edges.numpy()[:n_obs + 1],
                                         dtype=tf.float32)
            self.cases = tf.constant(full_cases.numpy()[:n_obs], dtype=tf.float32)

            if restart_optimiser:
                self.opt = tf.keras.optimizers.Adam(
                    tf.keras.optimizers.schedules.ExponentialDecay(
                        float(pc["learning_rate"]), int(pc["lr_decay_steps"]),
                        float(pc["lr_decay_rate"])))

            step = tf.function(self._step.python_function) \
                if hasattr(self._step, "python_function") else self._step
            for _ in range(epochs_per_stage):
                step()

        # Restore the full window so that evaluation uses all the data.
        self.tau_col_tf, self.b_col, self.zr_col = full_col, full_b, full_zr
        self.tau_edges, self.cases = full_edges, full_cases
        self.seconds = time.time() - t0
        self.history = np.zeros((1, 1))
        return self


def evaluate(label, pinn, notes=""):
    p = pinn.params()
    mu = pinn.weekly_mu()
    k = estimate_dispersion_k(DATA.cases, mu)
    ll = nb_loglik(DATA.cases, mu, k)
    _, physics, _ = pinn._loss()
    r0 = basic_reproduction_number(p["beta_0"], fixed)
    err = abs(r0 - R0_TRUTH) / R0_TRUTH * 100
    print(f"  {label:26s} R0 {r0:7.4f}  err {err:6.1f}%  "
          f"pop_frac {p['pop_frac']:.5f}  physics {float(physics):.2e}  "
          f"logL {ll:8.1f}  {pinn.seconds:5.0f} s")
    return dict(config=label, R0=round(r0, 4), R0_rel_err=round(err, 1),
                pop_frac=round(p["pop_frac"], 6), a_rain=round(p["a_rain"], 3),
                physics_loss=float(f"{float(physics):.4g}"),
                loglik=round(ll, 1), seconds=round(pinn.seconds, 1), notes=notes)


def cfg_with(**overrides):
    c = copy.deepcopy(BASE)
    c["inference"]["pinn"].update(overrides)
    return c


rows = []
CFG = cfg_with(data_weight=20_000.0)          # best weighting from step 07

print("expanding window, optimiser carried forward")
for seed in SEEDS:
    rows.append(evaluate(f"curriculum_seed{seed}",
                         CurriculumPINN(DATA, CFG, fixed, seed=seed)
                         .train_curriculum(restart_optimiser=False),
                         "expanding window"))

print("\nexpanding window, optimiser restarted at each stage")
for seed in SEEDS:
    rows.append(evaluate(f"curriculum_restart_seed{seed}",
                         CurriculumPINN(DATA, CFG, fixed, seed=seed)
                         .train_curriculum(restart_optimiser=True),
                         "expanding window, optimiser restarted"))

tab = pd.DataFrame(rows)
tab.to_csv(tables / "08_pinn_timemarching.csv", index=False)

print("\n" + "=" * 78)
for group, label in (("curriculum_seed", "expanding"),
                     ("curriculum_restart_seed", "expanding + restart")):
    sub = tab[tab["config"].str.startswith(group)]
    e = sub["R0_rel_err"]
    print(f"{label:22s} R0 error {e.min():5.1f}% to {e.max():5.1f}%, "
          f"median {e.median():5.1f}%   (classical {CLASSICAL_ERR}%)")

# Step 07's best configuration, for reference, across the same seeds.
PRIOR = [0.1, 8.0, 50.1, 41.1, 14.9]
print(f"{'step 07 best (global)':22s} R0 error {min(PRIOR):5.1f}% to "
      f"{max(PRIOR):5.1f}%, median {np.median(PRIOR):5.1f}%")

fig, ax = plt.subplots(figsize=(9, 4.6))
groups = {"global (step 07)": PRIOR,
          "expanding window":
              tab[tab["config"].str.startswith("curriculum_seed")]["R0_rel_err"].tolist(),
          "expanding + restart":
              tab[tab["config"].str.startswith("curriculum_restart")]["R0_rel_err"].tolist()}
positions = range(len(groups))
for i, (name, vals) in zip(positions, groups.items()):
    ax.scatter([i] * len(vals), vals, s=70, alpha=0.75, color="#8c1c13")
    ax.plot([i - 0.18, i + 0.18], [np.median(vals)] * 2, color="k", lw=2)
ax.axhline(CLASSICAL_ERR, color="#1f6f8b", ls="--", lw=1.5,
           label=f"classical estimator ({CLASSICAL_ERR}%)")
ax.set_xticks(list(positions))
ax.set_xticklabels(list(groups))
ax.set_yscale("log")
ax.set_ylabel("relative error in $R_0$ (%)")
ax.set_title("Time-window curriculum against global training, five seeds each\n"
             "(bars are medians)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(figures / "08_pinn_timemarching.png", dpi=150)
print(f"\nFigure: {figures / '08_pinn_timemarching.png'}")
print(f"Table:  {tables / '08_pinn_timemarching.csv'}")
