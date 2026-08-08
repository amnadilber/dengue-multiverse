"""
Pipeline step 11 — practical identifiability, and an independent check on R0.

The bootstrap showed which parameters have wide intervals. It does not show why,
and the two possible reasons call for different responses. A parameter can be
poorly determined because the likelihood is genuinely flat along it — no amount
of data of this kind will fix that — or because the surface is multimodal and the
optimiser lands in different basins on different resamples, which more careful
optimisation would fix.

Profiling separates them. Each parameter is fixed across a grid, everything else
refitted at each point, and the resulting curve inspected. A flat profile means
practical non-identifiability. A profile with a clear minimum but a long tail
means the parameter is identified but skewed. Two minima mean multimodality.

The conventional threshold is a rise of 1.92 log-likelihood units above the
optimum, which bounds an approximate 95% interval for one parameter.

The script also computes a model-free check on R0. From the exponential growth
rate r of the early epidemic and the generation interval implied by the fixed
parameters, R0 can be estimated without fitting the compartmental model at all.
If the two disagree badly, the model — not the data — is the problem.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.inference import (CLIMATE_PARAMS, NULL_PARAMS, Dataset,  # noqa: E402
                                 fit, nb_loglik, pack, predict, unpack)
from dengue_pk.models import (FixedParams, IntegrationFailure,  # noqa: E402
                              basic_reproduction_number)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402
from scipy.stats import linregress  # noqa: E402

from dengue_pk.inference import nb_deviance_residuals  # noqa: E402

cfg = load_config()
OBS = cfg["inference"]["observation"]
fixed = FixedParams.from_config(cfg)
processed = resolve(cfg, "processed")
figures = resolve(cfg, "figures")
tables = resolve(cfg, "tables")

CHI2_95 = 1.92          # half of the 95% chi-square quantile with one d.f.
N_GRID = 15

rows, curves = [], {}


def profile(data, model, base_theta, k_hat, param, grid):
    """Refit everything except ``param`` at each grid point."""
    names = CLIMATE_PARAMS if model == "climate" else NULL_PARAMS
    free = [n for n in names if n != param]
    out = []
    for value in grid:
        def residual(x):
            theta = unpack(x, free)
            theta[param] = value
            try:
                mu = predict(theta, data, fixed, model, cfg)
            except IntegrationFailure:
                return np.full(len(data.cases), 1e3)
            r = nb_deviance_residuals(data.cases, mu, k_hat)
            return np.nan_to_num(r, nan=1e3, posinf=1e3, neginf=-1e3)

        best_ll = -np.inf
        for jitter in (1.0, 0.6, 1.6):
            start = {n: base_theta[n] * jitter for n in free}
            try:
                sol = least_squares(residual, pack(start, free), method="lm",
                                    max_nfev=2500)
            except Exception:
                continue
            theta = unpack(sol.x, free)
            theta[param] = value
            try:
                ll = nb_loglik(data.cases, predict(theta, data, fixed, model, cfg),
                               k_hat)
            except IntegrationFailure:
                continue
            best_ll = max(best_ll, ll)
        out.append(best_ll)
    return np.array(out)


for name, w in cfg["windows"].items():
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])
    data = Dataset.from_frame(df, label=w["label"])
    print(f"\n{'=' * 74}\n{w['label']}\n{'=' * 74}")

    for model in ("climate", "constant"):
        res = fit(data, cfg, fixed, model=model, observation=OBS)
        k_hat = res.nb_k if res.nb_k else 1e6
        ll_max = nb_loglik(data.cases,
                           predict(res.theta, data, fixed, model, cfg), k_hat)
        print(f"\n  [{model}] optimum logL {ll_max:.2f}, "
              f"R0 = {basic_reproduction_number(res.theta['beta_0'], fixed):.3f}")

        for param in res.names:
            centre = res.theta[param]
            if param in ("beta_0", "pop_frac", "i0_frac"):
                grid = np.exp(np.linspace(np.log(centre) - 1.6,
                                          np.log(centre) + 1.6, N_GRID))
            elif param == "a_temp":
                grid = np.linspace(0.0, max(centre * 3, 1.5), N_GRID)
            else:
                grid = np.linspace(centre - 1.0, centre + 1.0, N_GRID)

            lls = profile(data, model, res.theta, k_hat, param, grid)
            drop = ll_max - lls
            within = np.isfinite(drop) & (drop <= CHI2_95)
            if within.any():
                lo, hi = grid[within].min(), grid[within].max()
                # A parameter whose whole grid stays inside the threshold is not
                # bounded by the data at all over the range examined.
                unbounded = within.all()
            else:
                lo = hi = np.nan
                unbounded = False

            span = (hi / lo if (np.isfinite(lo) and lo > 0 and hi > 0)
                    else np.nan)
            flag = ("UNBOUNDED over the grid" if unbounded
                    else "wide" if np.isfinite(span) and span > 5
                    else "identified")
            print(f"    {param:10s} profile 95% [{lo:.4g}, {hi:.4g}]  "
                  f"ratio {span:6.2f}  {flag}")

            rows.append(dict(window=name, model=model, parameter=param,
                             estimate=centre, prof_lo=lo, prof_hi=hi,
                             ratio=span, unbounded=unbounded,
                             max_drop=float(np.nanmax(drop)) if np.isfinite(drop).any() else np.nan))
            curves[(name, model, param)] = (grid, drop)

# ---------------------------------------------------------------------------
# Model-free check: exponential growth rate of the early epidemic
# ---------------------------------------------------------------------------
print(f"\n{'=' * 74}\nModel-free R0 from the early exponential growth rate\n{'=' * 74}")
print("R0 is bounded by the growth rate and the generation interval, but the")
print("bound depends on how that interval is distributed, and the two extremes")
print("are far apart:")
print("    fixed interval        R0 = exp(r T)      -- the upper bound")
print("    exponentially spread   R0 = 1 + r T       -- the lower bound")
print("Reporting only the first, as is common, overstates R0 substantially when")
print("r T is large. Both are given below, and the fitted value should fall")
print("between them if the model is describing the same epidemic.\n")

T_GEN = 1.0 / fixed.gamma_h + 1.0 / fixed.mu_v + 1.0 / fixed.sigma_h
print(f"generation interval from the fixed parameters: {T_GEN:.1f} days\n")

growth_rows = []
for name, w in cfg["windows"].items():
    df = pd.read_csv(processed / f"{name}.csv", parse_dates=["week_start"])
    cases = df["cases"].to_numpy(float)
    days = df["days_from_start"].to_numpy(float)
    peak = int(np.argmax(cases))
    # Growth phase: from the first week above 5% of the peak up to the peak.
    start = int(np.argmax(cases > 0.05 * cases.max()))
    sel = slice(start, max(peak, start + 4))
    x, yv = days[sel], cases[sel]
    keep = yv > 0
    if keep.sum() < 4:
        print(f"  {name}: too few growth-phase weeks")
        continue
    fitres = linregress(x[keep], np.log(yv[keep]))
    r = fitres.slope
    r0_upper = float(np.exp(r * T_GEN))         # fixed generation interval
    r0_lower = float(1.0 + r * T_GEN)           # exponentially distributed
    model_r0 = [row for row in rows
                if row["window"] == name and row["model"] == "constant"
                and row["parameter"] == "beta_0"]
    r0_model = (basic_reproduction_number(model_r0[0]["estimate"], fixed)
                if model_r0 else np.nan)
    inside = r0_lower <= r0_model <= r0_upper
    # The attack rate the fitted R0 implies, against what was observed. A
    # homogeneous model ties growth rate and final size together; if the
    # observed epidemic grows quickly and infects few, no single R0 satisfies
    # both, and that is a statement about the model rather than the data.
    print(f"  {name:15s} r = {r:.4f}/day over {int(keep.sum())} weeks "
          f"(R^2 {fitres.rvalue ** 2:.2f})")
    print(f"    {'':13s} model-free R0 between {r0_lower:.2f} and "
          f"{r0_upper:.2f};  fitted {r0_model:.2f}  "
          f"{'consistent' if inside else 'BELOW the lower bound'}")
    growth_rows.append(dict(window=name, growth_rate_per_day=round(r, 5),
                            r_squared=round(fitres.rvalue ** 2, 3),
                            weeks_used=int(keep.sum()),
                            R0_lower_bound=round(r0_lower, 3),
                            R0_upper_bound=round(r0_upper, 3),
                            R0_fitted=round(r0_model, 3),
                            consistent=bool(inside)))

pd.DataFrame(rows).to_csv(tables / "09_profile_identifiability.csv", index=False)
pd.DataFrame(growth_rows).to_csv(tables / "10_growth_rate_check.csv", index=False)

# ---------------------------------------------------------------------------
key_params = ["beta_0", "pop_frac", "a_temp", "a_rain"]
windows = list(cfg["windows"])
fig, axes = plt.subplots(len(windows), len(key_params),
                         figsize=(3.6 * len(key_params), 3.2 * len(windows)),
                         squeeze=False)
for i, wname in enumerate(windows):
    for j, param in enumerate(key_params):
        ax = axes[i][j]
        drawn = False
        for model, colour in (("climate", "#8c1c13"), ("constant", "#1f6f8b")):
            item = curves.get((wname, model, param))
            if item is None:
                continue
            grid, drop = item
            ax.plot(grid, drop, "o-", ms=3, color=colour, label=model)
            drawn = True
        if not drawn:
            ax.set_axis_off()
            continue
        ax.axhline(CHI2_95, color="k", ls="--", lw=1)
        ax.set_ylim(-0.2, 12)
        if param in ("beta_0", "pop_frac"):
            ax.set_xscale("log")
        if i == 0:
            ax.set_title(param)
        if j == 0:
            ax.set_ylabel(f"{wname}\ndrop in logL")
        ax.grid(alpha=0.25)
        if i == 0 and j == 0:
            ax.legend(fontsize=7)
fig.suptitle("Profile likelihoods — a curve staying below the dashed line is "
             "not bounded by the data", y=1.01)
fig.tight_layout()
fig.savefig(figures / "09_profile_identifiability.png", dpi=150)
print(f"\nFigure: {figures / '09_profile_identifiability.png'}")
print(f"Tables: {tables / '09_profile_identifiability.csv'}, "
      f"{tables / '10_growth_rate_check.csv'}")
