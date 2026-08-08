"""
Pipeline step 5 — confidence intervals by parametric bootstrap.

The standard errors reported so far come from the Jacobian at the optimum, which
assumes the likelihood is locally quadratic and the parameters are well
identified. Neither assumption is safe here: the thermal exponent was shown to
sit on a plateau, and transmission and the population at risk remain partially
confounded. Where those assumptions fail, asymptotic intervals are too narrow
and symmetric about a point that may not be central.

A parametric bootstrap makes no such assumption. Counts are simulated from the
fitted model under its own negative binomial observation process, the whole
estimation is repeated on each synthetic dataset, and the spread of the
resulting estimates is the sampling distribution. Percentile intervals then need
no normality.

Each replicate is warm-started from the point estimate and the dispersion is
held at its fitted value: replicates differ only by resampling noise, so the
full multi-start search would cost hours and change nothing. That choice is a
limitation worth stating — it assumes the optimiser would not find a distant
better optimum on a resampled dataset, which the multi-start agreement on the
real data makes plausible but does not prove.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.inference import Dataset, fit, predict  # noqa: E402
from dengue_pk.models import (FixedParams, IntegrationFailure,  # noqa: E402
                              basic_reproduction_number)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
OBS = cfg["inference"]["observation"]
N_BOOT = cfg["inference"]["uncertainty"]["n_bootstrap"]
fixed = FixedParams.from_config(cfg)
processed = resolve(cfg, "processed")
figures = resolve(cfg, "figures")
tables = resolve(cfg, "tables")

rng = np.random.default_rng(cfg["seed"])
rows, draws = [], []

for name, w in cfg["windows"].items():
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])
    data = Dataset.from_frame(df, label=w["label"])

    for model in ("climate", "constant"):
        print(f"\n=== {name} / {model} — {N_BOOT} bootstrap replicates ===")
        point = fit(data, cfg, fixed, model=model, observation=OBS)
        mu_hat = predict(point.theta, data, fixed, model, cfg)
        k_hat = point.nb_k if point.nb_k else 1e6
        print(f"  point estimate: R0 = "
              f"{basic_reproduction_number(point.theta['beta_0'], fixed):.3f}, "
              f"k = {k_hat:.2f}")

        t0 = time.time()
        ok = 0
        for b in range(N_BOOT):
            # Simulate under the fitted negative binomial: a gamma-Poisson
            # mixture reproduces NB2 with variance mu + mu^2/k.
            lam = rng.gamma(shape=k_hat, scale=np.maximum(mu_hat, 1e-12) / k_hat)
            synth = Dataset(data.days, rng.poisson(lam).astype(float),
                            data.temp_c, data.z_rain, data.population, name)
            try:
                res = fit(synth, cfg, fixed, model=model, observation=OBS,
                          start_from=point.theta, n_starts_override=2,
                          fixed_nb_k=k_hat, seed=cfg["seed"] + b)
            except (RuntimeError, IntegrationFailure):
                continue
            ok += 1
            draws.append(dict(
                window=name, model=model, replicate=b,
                R0=basic_reproduction_number(res.theta["beta_0"], fixed),
                pop_at_risk_M=res.theta["pop_frac"] * data.population / 1e6,
                **{p: res.theta[p] for p in res.names}))
            if (b + 1) % 25 == 0:
                el = time.time() - t0
                print(f"  {b + 1}/{N_BOOT} replicates, {el / (b + 1):.1f} s each, "
                      f"{el / 60:.1f} min elapsed")

        sub = pd.DataFrame([d for d in draws
                            if d["window"] == name and d["model"] == model])
        print(f"  {ok}/{N_BOOT} replicates converged in "
              f"{(time.time() - t0) / 60:.1f} min")

        for param in ["R0", "pop_at_risk_M"] + list(point.names):
            if param not in sub.columns:
                continue
            v = sub[param].to_numpy(float)
            v = v[np.isfinite(v)]
            if len(v) < 10:
                continue
            pt = (basic_reproduction_number(point.theta["beta_0"], fixed)
                  if param == "R0" else
                  point.theta["pop_frac"] * data.population / 1e6
                  if param == "pop_at_risk_M" else point.theta[param])
            lo, hi = np.percentile(v, [2.5, 97.5])
            rows.append(dict(window=name, model=model, parameter=param,
                             point=pt, boot_median=float(np.median(v)),
                             ci_lo=float(lo), ci_hi=float(hi),
                             boot_sd=float(np.std(v)),
                             asympt_se=point.stderr.get(param, np.nan),
                             n_replicates=len(v)))

tab = pd.DataFrame(rows)
tab.to_csv(tables / "04_bootstrap_ci.csv", index=False)
pd.DataFrame(draws).to_csv(tables / "04_bootstrap_draws.csv", index=False)

print("\n" + "=" * 78)
show = tab[tab["parameter"].isin(["R0", "pop_at_risk_M", "a_temp", "a_rain"])]
print(show.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

print("\nAsymptotic standard errors versus the bootstrap:")
for _, r in tab.iterrows():
    if np.isfinite(r["asympt_se"]) and r["asympt_se"] > 0:
        ratio = r["boot_sd"] / r["asympt_se"]
        if ratio > 1.5 or ratio < 0.67:
            print(f"  {r['window']}/{r['model']}/{r['parameter']}: bootstrap SD "
                  f"is {ratio:.1f}x the asymptotic one — the asymptotic interval "
                  f"is {'too narrow' if ratio > 1 else 'too wide'}")

# --- Figure: bootstrap distribution of R0 per window ------------------------
d = pd.DataFrame(draws)
if not d.empty:
    windows = list(cfg["windows"])
    fig, axes = plt.subplots(1, len(windows), figsize=(4.6 * len(windows), 4))
    for ax, name in zip(np.atleast_1d(axes), windows):
        for model, colour in (("climate", "#8c1c13"), ("constant", "#1f6f8b")):
            v = d[(d["window"] == name) & (d["model"] == model)]["R0"]
            v = v[np.isfinite(v)]
            if len(v):
                ax.hist(v, bins=30, alpha=0.55, color=colour, label=model)
        ax.axvline(1.0, color="k", ls=":", lw=1)
        ax.set_title(f"{name}\nbootstrap distribution of $R_0$")
        ax.set_xlabel("$R_0$")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures / "04_bootstrap_R0.png", dpi=150)
    print(f"\nFigure: {figures / '04_bootstrap_R0.png'}")

print(f"Tables: {tables / '04_bootstrap_ci.csv'}, "
      f"{tables / '04_bootstrap_draws.csv'}")
