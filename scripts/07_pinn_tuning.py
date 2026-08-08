"""
Pipeline step 7 — a systematic attempt to make the PINN work before concluding
that it cannot.

The first PINN configuration failed badly on synthetic data: it recovered
R0 = 0.02 against a truth of 1.70, driving transmission to zero while inflating
the population at risk fortyfold, and its physics residual never approached zero.
Reporting that as "PINNs do not work here" after a single attempt would be as
unsound as reporting a spurious positive — the result would say more about the
configuration than about the method.

This script therefore varies the choices most likely to be responsible and
records every outcome, successful or not. Each configuration is judged on
synthetic data where the truth is known:

* **Loss weighting.** If the data term is too weak relative to the physics
  residual, the network can satisfy neither and settle between them. The failure
  signature — transmission collapsing while the population inflates — is what a
  weakly constrained scale parameter looks like.
* **Two-stage training.** Fitting the network to satisfy the equations at fixed
  parameters first, then releasing the parameters, avoids asking the optimiser to
  learn a solution and its coefficients simultaneously from a random start.
* **Sequential time-marching.** Short windows solved in order, each starting from
  the previous window's endpoint. This is the documented remedy for PINN failure
  on problems with sharp transients (Krishnapriyan et al., 2021).
* **Collocation density and training length.** The cheapest explanations, and so
  worth eliminating rather than assuming.

Whatever the outcome, the table this produces is the evidence for the claim the
paper eventually makes.
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

BASE = load_config()
fixed = FixedParams.from_config(BASE)
processed = resolve(BASE, "processed")
figures = resolve(BASE, "figures")
tables = resolve(BASE, "tables")

# Same calibrated synthetic truth as step 06, so results are comparable.
TRUTH = dict(beta_0=0.17, a_temp=0.2, a_rain=0.1, pop_frac=0.0025,
             i0_frac=5e-6)
R0_TRUTH = basic_reproduction_number(TRUTH["beta_0"], fixed)

ref = pd.read_csv(processed / "national_2013.csv", parse_dates=["week_start"])
template = Dataset.from_frame(ref, label="synthetic")
mu_true = predict(TRUTH, template, fixed, "climate", BASE)

rng = np.random.default_rng(BASE["seed"])
lam = rng.gamma(shape=5.0, scale=np.maximum(mu_true, 1e-12) / 5.0)
DATA = Dataset(template.days, rng.poisson(lam).astype(float), template.temp_c,
               template.z_rain, template.population, "synthetic")

print(f"synthetic truth: R0 = {R0_TRUTH:.3f}, pop_frac = {TRUTH['pop_frac']}, "
      f"a_rain = {TRUTH['a_rain']}")
print(f"synthetic data: {DATA.cases.sum():.0f} cases, peak "
      f"{DATA.cases.max():.0f}\n")


def cfg_with(**overrides) -> dict:
    cfg = copy.deepcopy(BASE)
    cfg["inference"]["pinn"].update(overrides)
    return cfg


PARTIAL = tables / "06_pinn_tuning.csv"


def evaluate(label: str, pinn: InversePINN, notes: str = "") -> dict:
    p = pinn.params()
    mu = pinn.weekly_mu()
    k = estimate_dispersion_k(DATA.cases, mu)
    ll = nb_loglik(DATA.cases, mu, k)
    _, physics, data_loss = pinn._loss()
    r0 = basic_reproduction_number(p["beta_0"], fixed)
    row = dict(config=label, R0=round(r0, 4), R0_truth=round(R0_TRUTH, 3),
               R0_rel_err=round(abs(r0 - R0_TRUTH) / R0_TRUTH * 100, 1),
               pop_frac=round(p["pop_frac"], 6),
               pop_frac_truth=TRUTH["pop_frac"],
               a_rain=round(p["a_rain"], 3),
               physics_loss=float(f"{float(physics):.4g}"),
               loglik=round(ll, 1), seconds=round(pinn.seconds, 1), notes=notes)
    print(f"  -> R0 {r0:8.4f} (err {row['R0_rel_err']:6.1f}%)   "
          f"pop_frac {p['pop_frac']:.5f}   physics {float(physics):.3e}   "
          f"logL {ll:9.1f}   {pinn.seconds:.0f} s")
    # Written after every configuration rather than at the end: a crash in a
    # later configuration should not discard the ones already run, and each takes
    # minutes.
    rows.append(row)
    pd.DataFrame(rows).to_csv(PARTIAL, index=False)
    return row


rows = []

# ---------------------------------------------------------------------------
# A. Baseline, as run in step 06
# ---------------------------------------------------------------------------
print("A. baseline (data_weight 200, 40k epochs)")
evaluate("A_baseline",
         InversePINN(DATA, cfg_with(), fixed).train(verbose_every=0),
         "as in step 06")

# ---------------------------------------------------------------------------
# B. Loss weighting sweep
# ---------------------------------------------------------------------------
for w in (2_000.0, 20_000.0, 200_000.0):
    print(f"B. data_weight = {w:,.0f}")
    cfg = cfg_with(data_weight=w)
    evaluate(f"B_dataweight_{int(w)}",
             InversePINN(DATA, cfg, fixed).train(verbose_every=0),
             "stronger data term")

# ---------------------------------------------------------------------------
# C. Denser collocation and longer training
# ---------------------------------------------------------------------------
print("C. 1200 collocation points, 80k epochs")
cfg = cfg_with(collocation_points=1200, epochs=80_000)
evaluate("C_dense_long",
         InversePINN(DATA, cfg, fixed).train(verbose_every=0),
         "denser grid, longer run")

# ---------------------------------------------------------------------------
# D. Two-stage: satisfy the physics first, then release the parameters
# ---------------------------------------------------------------------------
print("D. two-stage (physics first at fixed parameters, then joint)")


class TwoStagePINN(InversePINN):
    """Freeze the transmission parameters for an initial phase.

    Asking the optimiser to learn both the solution and its coefficients from a
    random initialisation means the physics residual is being computed against a
    solution that is meaningless early on, so the parameter gradients early in
    training are noise. Holding the parameters still until the network represents
    *some* valid trajectory gives them a meaningful gradient when released.
    """

    def train_two_stage(self, freeze_epochs: int):
        import tensorflow as tf

        t0 = time.time()
        n_net = len(self.net.trainable_variables)

        # The optimiser must see the same variable list throughout: a Keras 3
        # optimiser is built against the variables of its first call and rejects
        # any others afterwards. Freezing is therefore done by zeroing the
        # parameter gradients rather than by omitting the parameters, which keeps
        # the variable list constant and leaves the optimiser state valid when
        # they are released.
        @tf.function
        def step(freeze: bool):
            variables = (self.net.trainable_variables
                         + [self.log_beta0, self.a_rain, self.logit_popfrac,
                            self.log_i0]
                         + ([self.log_atemp] if self.model == "climate" else []))
            with tf.GradientTape() as tape:
                total, _, _ = self._loss()
            grads = tape.gradient(total, variables)
            if freeze:
                grads = grads[:n_net] + [tf.zeros_like(g) for g in grads[n_net:]]
            self.opt.apply_gradients(zip(grads, variables))
            return total

        for _ in range(freeze_epochs):
            step(True)
        for _ in range(self.epochs - freeze_epochs):
            step(False)
        self.seconds = time.time() - t0
        self.history = np.zeros((1, 1))
        return self


cfg = cfg_with(data_weight=20_000.0)
evaluate("D_two_stage",
         TwoStagePINN(DATA, cfg, fixed).train_two_stage(15_000),
         "15k epochs with parameters frozen")

# ---------------------------------------------------------------------------
# E. Wider network — is capacity the constraint?
# ---------------------------------------------------------------------------
print("E. wider network (128 x 4), data_weight 20,000")
cfg = cfg_with(layers=[128, 128, 128, 128], data_weight=20_000.0)
evaluate("E_wide",
         InversePINN(DATA, cfg, fixed).train(verbose_every=0),
         "more capacity")

# ---------------------------------------------------------------------------
# F. Is the best weighting reproducible, or was one run lucky?
#
# This is the decisive question. If a configuration recovers the truth on one
# seed and not on others, it has not solved the problem — it has sampled a
# favourable initialisation, and a practitioner with real data and no truth to
# check against would have no way to tell which run they had got.
# ---------------------------------------------------------------------------
print("F. best weighting (20,000) repeated across seeds")
for seed in (1, 2, 3, 4):
    cfg = cfg_with(data_weight=20_000.0)
    evaluate(f"F_seed{seed}",
             InversePINN(DATA, cfg, fixed, seed=seed).train(verbose_every=0),
             "reproducibility of the best weighting")

# ---------------------------------------------------------------------------
tab = pd.DataFrame(rows)
tab.to_csv(PARTIAL, index=False)
print("\n" + "=" * 92)
print(tab.to_string(index=False))
print(f"\nTable: {tables / '06_pinn_tuning.csv'}")

best = tab.loc[tab["R0_rel_err"].idxmin()]
print(f"\nBest configuration: {best['config']} — R0 error {best['R0_rel_err']}%")
print("Classical estimator on the same data, for reference: R0 error 3.1% "
      "(R0 = 1.752 against a truth of 1.700)")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
axes[0].barh(tab["config"], tab["R0_rel_err"], color="#8c1c13")
axes[0].axvline(3.1, color="k", ls="--", lw=1, label="classical estimator")
axes[0].set_xlabel("relative error in $R_0$ (%)")
axes[0].set_xscale("log")
axes[0].set_title("Parameter recovery on synthetic data")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3, axis="x")

axes[1].scatter(tab["physics_loss"], tab["R0_rel_err"], s=60, color="#1f6f8b")
for _, r in tab.iterrows():
    axes[1].annotate(r["config"].split("_", 1)[0],
                     (r["physics_loss"], r["R0_rel_err"]), fontsize=8,
                     xytext=(4, 3), textcoords="offset points")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlabel("final physics residual")
axes[1].set_ylabel("relative error in $R_0$ (%)")
axes[1].set_title("Does a smaller residual mean a better estimate?")
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(figures / "06_pinn_tuning.png", dpi=150)
print(f"Figure: {figures / '06_pinn_tuning.png'}")
