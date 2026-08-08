"""
Pipeline step 3 — fit the climate-forced and constant-transmission models to
every window by classical inversion.

Both models are fitted to every window, and both are compared on two criteria
that answer different questions:

* **AIC** asks whether the climate terms earn their parameters on the data used
  to fit them.
* **Held-out forecast error** asks whether they help predict weeks the fit never
  saw. A model can win on AIC and still forecast worse, and only the second
  question matters for whether the climate terms describe a mechanism rather
  than absorbing noise.

The held-out split is chronological. Randomly withholding weeks would let the
model interpolate between observations on either side, which is not forecasting
and would flatter both models.

Results are written as CSV tables and a figure; nothing is printed that is not
also saved.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.inference import (Dataset, fit, nb_deviance_residuals,  # noqa: E402
                                 nb_loglik, poisson_deviance_residuals,
                                 poisson_loglik, predict)
from dengue_pk.models import FixedParams, basic_reproduction_number  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

TRAIN_FRAC = 0.75

cfg = load_config()
OBSERVATION = cfg["inference"]["observation"]
fixed = FixedParams.from_config(cfg)
processed = resolve(cfg, "processed")
figures = resolve(cfg, "figures")
tables = resolve(cfg, "tables")

param_rows, compare_rows = [], []
windows = list(cfg["windows"].items())
fig, axes = plt.subplots(len(windows), 1, figsize=(11, 3.8 * len(windows)))

for ax, (name, w) in zip(np.atleast_1d(axes), windows):
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])
    full = Dataset.from_frame(df, label=w["label"])
    k_train = int(round(len(full.days) * TRAIN_FRAC))

    print(f"\n=== {w['label']} ({len(full.days)} weeks, "
          f"{k_train} for fitting, {len(full.days) - k_train} held out) ===")

    train = full.head(k_train)

    fits, curves = {}, {}
    for model in ("climate", "constant"):
        # Two fits per model, for two different purposes, and conflating them
        # would be a reporting error rather than a modelling one:
        #
        #   * the FULL-window fit is the estimate. It uses all the data, and it
        #     is what the bootstrap and the profile likelihoods are computed
        #     from, so it is what the paper must quote.
        #   * the TRAINING fit exists only to score held-out weeks. Its
        #     parameters are not reported.
        #
        # They differ materially: nationally the climate model's thermal
        # exponent is 0.215 on the full window, with a profile interval
        # excluding zero, but 0.023 +/- 0.033 on the first 75%, which includes
        # it. Quoting whichever supports the preferred conclusion would be
        # indefensible, so both are printed and only one is reported.
        res_full = fit(full, cfg, fixed, model=model, observation=OBSERVATION)
        res = fit(train, cfg, fixed, model=model, observation=OBSERVATION)
        fits[model] = res_full
        fits[model + "_train"] = res
        print(f"  [full window — reported estimate]")
        print("  " + res_full.summary().replace("\n", "\n  "))
        print(f"  [first {k_train} weeks — used only to score the held-out tail]")
        print(f"    R0 = {basic_reproduction_number(res.theta['beta_0'], fixed):.3f}, "
              f"AIC {res.aic:.1f}")
        print(f"    {res.n_forward_solves + res_full.n_forward_solves:,} forward "
              f"solves in {res.seconds + res_full.seconds:.1f} s")

        # Predict across the whole window, including the held-out tail. The
        # forcing is rebuilt on the full climate record, so the model is
        # genuinely extrapolating in time rather than being refitted.
        mu_train_fit = predict(res.theta, full, fixed, model, cfg)
        # Curves shown are from the full-window fit, matching the reported
        # estimate; the training fit's curve is used only for scoring.
        curves[model] = predict(res_full.theta, full, fixed, model, cfg)
        mu_full = mu_train_fit

        # Held-out scores are computed under the same observation model the fit
        # used, with the dispersion held at its fitted value: re-estimating it on
        # the held-out weeks would let the model tune itself to the very data
        # being used to judge it.
        held_out = slice(k_train, None)
        if OBSERVATION == "nb":
            dev_out = float(np.sum(nb_deviance_residuals(
                full.cases[held_out], mu_full[held_out], res.nb_k) ** 2))
            ll_out = nb_loglik(full.cases[held_out], mu_full[held_out], res.nb_k)
        else:
            dev_out = float(np.sum(poisson_deviance_residuals(
                full.cases[held_out], mu_full[held_out]) ** 2))
            ll_out = poisson_loglik(full.cases[held_out], mu_full[held_out])
        mae_out = float(np.mean(np.abs(full.cases[held_out] - mu_full[held_out])))

        r0 = basic_reproduction_number(res_full.theta["beta_0"], fixed)
        at_risk = res_full.theta["pop_frac"] * full.population
        compare_rows.append(dict(
            window=name, model=model, aic=round(res_full.aic, 1),
            aic_trainfit=round(res.aic, 1),
            loglik=round(res_full.loglik, 1),
            dispersion=round(res_full.dispersion, 1),
            heldout_deviance=round(dev_out, 1), heldout_loglik=round(ll_out, 1),
            heldout_mae=round(mae_out, 1),
            R0_optimum=round(r0, 2),
            pop_at_risk_millions=round(at_risk / 1e6, 2),
            pop_frac=round(res_full.theta["pop_frac"], 5),
            a_temp=round(res_full.theta.get("a_temp", float("nan")), 3),
            a_rain=round(res_full.theta.get("a_rain", float("nan")), 3),
            nb_k=round(res_full.nb_k, 2) if res_full.nb_k else float("nan"),
            starts_converged=f"{res_full.n_converged}/{res_full.n_starts}",
            forward_solves=res_full.n_forward_solves,
            seconds=round(res_full.seconds, 1)))

        for p in res_full.names:
            param_rows.append(dict(window=name, model=model, parameter=p,
                                   estimate=res_full.theta[p],
                                   stderr=res_full.stderr.get(p, float("nan")),
                                   across_start_sd=res_full.start_spread.get(p),
                                   estimate_trainfit=res.theta[p]))

    # --- Correlation between transmission and the population at risk ---------
    # These two remain partially confounded even with the reporting fraction
    # fixed: a larger epidemic in a smaller catchment resembles a smaller one in
    # a larger catchment. How partially is worth reporting.
    res_c = fits["climate"]
    if res_c.corr is not None:
        i, j = res_c.names.index("beta_0"), res_c.names.index("pop_frac")
        print(f"    corr(beta_0, pop_frac) = {res_c.corr[i, j]:+.3f}"
              "   (near +/-1 would mean they are not separately identified)")

    ax.bar(df["week_start"], full.cases, width=6, color="#c9ccd1",
           label="reported cases")
    ax.plot(df["week_start"], curves["climate"], color="#8c1c13", lw=2,
            label="climate-forced")
    ax.plot(df["week_start"], curves["constant"], color="#1f6f8b", lw=1.6, ls="--",
            label="constant transmission")
    ax.axvline(df["week_start"].iloc[k_train], color="k", lw=1, ls=":")
    ax.annotate("held out →", (df["week_start"].iloc[k_train], ax.get_ylim()[1]),
                fontsize=8, va="top", ha="left")
    ax.set_title(f"{w['label']}")
    ax.set_ylabel("cases / week")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, axis="y")

fig.tight_layout()
fig.savefig(figures / "02_classical_fits.png", dpi=150)

pd.DataFrame(compare_rows).to_csv(tables / "02_model_comparison.csv", index=False)
pd.DataFrame(param_rows).to_csv(tables / "02_parameters.csv", index=False)

print("\n" + "=" * 72)
comp = pd.DataFrame(compare_rows)
print(comp[["window", "model", "aic", "aic_trainfit", "heldout_deviance",
            "heldout_mae", "R0_optimum", "pop_at_risk_millions",
            "a_temp", "a_rain", "nb_k"]].to_string(index=False))

print("\nModel selection:")
for name in comp["window"].unique():
    sub = comp[comp["window"] == name].set_index("model")
    d_aic = sub.loc["climate", "aic"] - sub.loc["constant", "aic"]
    d_aic_train = (sub.loc["climate", "aic_trainfit"]
                   - sub.loc["constant", "aic_trainfit"])
    d_out = sub.loc["climate", "heldout_deviance"] - sub.loc["constant", "heldout_deviance"]
    verdict = ("climate forcing is supported both in-sample and out-of-sample"
               if d_aic < 0 and d_out < 0 else
               "climate forcing wins in-sample but not out-of-sample"
               if d_aic < 0 else
               "climate forcing is not supported")
    print(f"  {name}: dAIC {d_aic:+.1f} (full window), "
          f"{d_aic_train:+.1f} (training fit), "
          f"held-out deviance {d_out:+.1f} — {verdict}")

print(f"\nFigure: {figures / '02_classical_fits.png'}")
print(f"Tables: {tables / '02_model_comparison.csv'}, {tables / '02_parameters.csv'}")
