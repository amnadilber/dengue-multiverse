"""
Pipeline step 13 — testing the explanation, not just the observation.

Two results point at spatial heterogeneity: aggregation depresses the fitted
transmission rate, and the observed growth rates exceed what the fitted rates
permit. Heterogeneity has been offered as the reason for both. This script asks
whether that reason survives being tested.

A two-patch model is fitted to the Khyber Pakhtunkhwa aggregate — the case with
the largest aggregation bias, −29% — using only the aggregate series. It is told
nothing about which districts contributed or when they peaked. Three predictions
follow if heterogeneity is the explanation:

  1. it should fit materially better than one patch;
  2. its recovered transmission rates should bracket the separately fitted
     district values (1.25 and 1.96) rather than falling below both, as the
     single-patch fit does at 1.14;
  3. it should reconcile growth with final size, so its faster patch should be
     consistent with the growth-rate bound of 1.79 that the single-patch fit
     violates.

A failure on any of these would mean the explanation is wrong, and the paper
would have to describe the aggregation bias without claiming to know its cause.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.climate import lagged_smoothed_rain, standardise  # noqa: E402
from dengue_pk.inference import (Dataset, estimate_dispersion_k,  # noqa: E402
                                 nb_deviance_residuals, nb_loglik, fit)
from dengue_pk.models import (FixedParams, IntegrationFailure,  # noqa: E402
                              basic_reproduction_number)
from dengue_pk.multipatch import patch_names, predict_multipatch  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

cfg = load_config()
fixed = FixedParams.from_config(cfg)
RHO = cfg["model"]["fixed"]["rho_fixed"]
raw = resolve(cfg, "raw")
figures = resolve(cfg, "figures")
tables = resolve(cfg, "tables")

# Established in step 09.
DISTRICT_R0 = {"HARIPUR": 1.251, "LAKKI MARAWAT": 1.957}
SINGLE_PATCH_R0 = 1.142
GROWTH_LOWER_BOUND = 1.79       # step 11, KP

NAME_VARIANTS = {"LAKKIMARWAT": "LAKKI MARAWAT", "LAKIMARWAT": "LAKKI MARAWAT"}
POP_M = {"HARIPUR": 1.00, "LAKKI MARAWAT": 0.88}


def canonical(n):
    n = " ".join(str(n).upper().split())
    return NAME_VARIANTS.get(n, n)


# --- rebuild the KP aggregate exactly as step 09 did ------------------------
pk = pd.read_csv(raw / cfg["data"]["opendengue"]["csv_name"],
                 usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                          "calendar_start_date", "calendar_end_date",
                          "dengue_total"], low_memory=False)
pk = pk[pk["adm_0_name"].astype(str).str.upper() == "PAKISTAN"].copy()
pk["start"] = pd.to_datetime(pk["calendar_start_date"], errors="coerce")
pk["end"] = pd.to_datetime(pk["calendar_end_date"], errors="coerce")
pk = pk[(pk["end"] - pk["start"]).dt.days + 1 == 7]
prov = pk[pk["adm_1_name"].astype(str).str.upper() == "KHYBER PAKHTUNKHWA"].copy()
prov["district"] = prov["adm_2_name"].map(canonical)

series = {}
for nm in DISTRICT_R0:
    s = (prov[prov["district"] == nm].groupby("start")["dengue_total"].sum()
         .sort_index())
    series[nm] = s
weeks = pd.Series(sorted(set.intersection(*(set(s.index) for s in series.values()))))
agg = sum(s.reindex(weeks).to_numpy() for s in series.values())
days = (weeks - weeks.min()).dt.days.to_numpy(float)
population = sum(POP_M.values()) * 1e6

print(f"KP aggregate: {len(weeks)} weeks, {agg.sum():.0f} cases")
print(f"single-patch fit gave R0 = {SINGLE_PATCH_R0}, below both districts "
      f"({DISTRICT_R0['HARIPUR']}, {DISTRICT_R0['LAKKI MARAWAT']})")
print(f"growth-rate lower bound for this window: {GROWTH_LOWER_BOUND}\n")


# --- fitting machinery ------------------------------------------------------
def unpack_patches(x, names):
    """Transform to physically admissible values.

    The population fraction is bounded to (0, 1) by a logit rather than merely
    kept positive by a log. An unbounded first attempt placed 279 million people
    in a province of 1.9 million, producing a sub-critical patch that acted as a
    slowly decaying background rather than an epidemic — a better fit with no
    physical meaning. A model comparison won on such a fit would be worthless,
    so the constraint is part of the hypothesis being tested, not a convenience.
    """
    theta = {}
    for n, v in zip(names, x):
        if n.startswith("offset_days"):
            theta[n] = float(np.clip(v, -60.0, 120.0))
        elif n.startswith("pop_frac"):
            theta[n] = float(1.0 / (1.0 + np.exp(-np.clip(v, -30.0, 30.0))))
        else:
            theta[n] = float(np.exp(np.clip(v, -25.0, 5.0)))
    return theta


def pack_patches(theta, names):
    out = []
    for n in names:
        v = theta[n]
        if n.startswith("offset_days"):
            out.append(v)
        elif n.startswith("pop_frac"):
            p = min(max(v, 1e-9), 1 - 1e-9)
            out.append(np.log(p / (1 - p)))
        else:
            out.append(np.log(v))
    return np.array(out)


def fit_patches(n_patches, k_hat, n_starts=24, seed=0):
    names = patch_names(n_patches)
    rng = np.random.default_rng(seed)
    best = None

    def residual(x):
        theta = unpack_patches(x, names)
        try:
            mu = predict_multipatch(theta, days, population, fixed, RHO,
                                    n_patches)
        except IntegrationFailure:
            return np.full(len(agg), 1e3)
        r = nb_deviance_residuals(agg, mu, k_hat)
        return np.nan_to_num(r, nan=1e3, posinf=1e3, neginf=-1e3)

    for i in range(n_starts):
        theta0 = {}
        for k in range(n_patches):
            theta0[f"beta_{k}"] = float(np.exp(rng.uniform(np.log(0.05),
                                                           np.log(0.4))))
            theta0[f"pop_frac_{k}"] = float(np.exp(rng.uniform(np.log(1e-3),
                                                              np.log(0.3))))
            theta0[f"i0_frac_{k}"] = float(np.exp(rng.uniform(np.log(1e-7),
                                                             np.log(1e-2))))
            if k > 0:
                theta0[f"offset_days_{k}"] = float(rng.uniform(-20, 80))
        try:
            sol = least_squares(residual, pack_patches(theta0, names),
                                method="lm", max_nfev=6000)
        except Exception:
            continue
        cost = float(np.sum(sol.fun ** 2))
        if best is None or cost < best[0]:
            best = (cost, sol)
    if best is None:
        raise RuntimeError(f"no start converged for {n_patches} patches")
    theta = unpack_patches(best[1].x, names)
    mu = predict_multipatch(theta, days, population, fixed, RHO, n_patches)
    return theta, mu, names


# Dispersion from the single-patch fit, held fixed so the comparison is on the
# mean structure alone rather than partly on a refitted variance.
data1 = Dataset(days, agg, np.full(len(days), 25.0), np.zeros(len(days)),
                population, "kp_aggregate")
base = fit(data1, cfg, fixed, model="constant", observation="nb")
K_HAT = base.nb_k if base.nb_k else 5.0
print(f"dispersion held at k = {K_HAT:.2f} for both models\n")

rows = []
results = {}
for n_patches in (1, 2):
    t0 = time.time()
    theta, mu, names = fit_patches(n_patches, K_HAT)
    ll = nb_loglik(agg, mu, K_HAT)
    aic = -2 * ll + 2 * len(names)
    r0s = [basic_reproduction_number(theta[f"beta_{k}"], fixed)
           for k in range(n_patches)]
    print(f"{n_patches} patch{'es' if n_patches > 1 else ''}: logL {ll:8.2f}  "
          f"AIC {aic:7.2f}  R0 {[f'{v:.3f}' for v in r0s]}  "
          f"({time.time() - t0:.0f} s)")
    total_frac = sum(theta[f"pop_frac_{k}"] for k in range(n_patches))
    for k in range(n_patches):
        print(f"    patch {k}: R0 {r0s[k]:.3f}, "
              f"pop {theta[f'pop_frac_{k}'] * population / 1e6:.3f} M "
              f"({theta[f'pop_frac_{k}'] * 100:.1f}% of the province), "
              f"offset {theta.get(f'offset_days_{k}', 0.0):+.0f} d")
    if total_frac > 1.0:
        print(f"    WARNING: patches together cover {total_frac * 100:.0f}% "
              f"of the province; they are meant to partition it")
    results[n_patches] = (theta, mu, r0s, ll, aic)
    rows.append(dict(n_patches=n_patches, loglik=round(ll, 2),
                     aic=round(aic, 2), n_params=len(names),
                     R0=";".join(f"{v:.4f}" for v in r0s)))

# --- the three predictions --------------------------------------------------
theta2, mu2, r0_2, ll2, aic2 = results[2]
theta1, mu1, r0_1, ll1, aic1 = results[1]
lo_d, hi_d = min(DISTRICT_R0.values()), max(DISTRICT_R0.values())

print(f"\n{'=' * 74}\nDoes heterogeneity explain the bias?\n{'=' * 74}")

d_aic = aic2 - aic1
p1 = d_aic < -4
print(f"1. Fit: dAIC = {d_aic:+.1f} for two patches over one   "
      f"-> {'YES' if p1 else 'NO'}")

fast, slow = max(r0_2), min(r0_2)
p2 = fast > SINGLE_PATCH_R0 and fast >= lo_d * 0.8
print(f"2. Rates: two-patch R0 {sorted(f'{v:.2f}' for v in r0_2)} against "
      f"districts [{lo_d:.2f}, {hi_d:.2f}]")
print(f"   single-patch was {SINGLE_PATCH_R0:.2f}, below both   "
      f"-> {'YES' if p2 else 'NO'}")

p3 = fast >= GROWTH_LOWER_BOUND * 0.9
print(f"3. Growth: faster patch R0 = {fast:.2f} against the growth-rate lower "
      f"bound {GROWTH_LOWER_BOUND}   -> {'YES' if p3 else 'NO'}")

verdict = ("Heterogeneity accounts for the aggregation bias."
           if p1 and p2 and p3 else
           "Heterogeneity accounts for it partially; see which tests failed."
           if p1 else
           "Heterogeneity does NOT account for it; the explanation must be "
           "withdrawn or revised.")
print(f"\nVerdict: {verdict}")

rows.append(dict(n_patches="test", loglik=None, aic=None, n_params=None,
                 R0=None, delta_aic=round(d_aic, 2),
                 better_fit=bool(p1), brackets_districts=bool(p2),
                 meets_growth_bound=bool(p3), verdict=verdict))
pd.DataFrame(rows).to_csv(tables / "12_heterogeneity_test.csv", index=False)

# --- figure -----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axes[0]
ax.bar(weeks, agg, width=6, color="#c9ccd1", label="KP aggregate")
ax.plot(weeks, mu1, color="#1f6f8b", lw=2, ls="--",
        label=f"one patch ($R_0$ {r0_1[0]:.2f})")
ax.plot(weeks, mu2, color="#8c1c13", lw=2,
        label=f"two patches ($R_0$ {max(r0_2):.2f}, {min(r0_2):.2f})")
ax.set_ylabel("cases / week")
ax.set_title("Fitting the aggregate with one and two patches")
ax.legend(fontsize=8)
ax.grid(alpha=0.25, axis="y")

ax = axes[1]
labels = ["Haripur\n(fitted alone)", "Lakki Marawat\n(fitted alone)",
          "one patch\non the sum", "two patches\non the sum"]
vals = [DISTRICT_R0["HARIPUR"], DISTRICT_R0["LAKKI MARAWAT"],
        r0_1[0], max(r0_2)]
colours = ["#c9ccd1", "#c9ccd1", "#1f6f8b", "#8c1c13"]
ax.bar(labels, vals, color=colours)
ax.axhline(GROWTH_LOWER_BOUND, color="k", ls=":", lw=1.5,
           label=f"growth-rate lower bound ({GROWTH_LOWER_BOUND})")
ax.set_ylabel("$R_0$")
ax.set_title("Does the two-patch fit recover what aggregation hid?")
ax.legend(fontsize=8)
ax.grid(alpha=0.25, axis="y")
fig.tight_layout()
fig.savefig(figures / "11_heterogeneity_test.png", dpi=150)
print(f"\nFigure: {figures / '11_heterogeneity_test.png'}")
print(f"Table:  {tables / '12_heterogeneity_test.csv'}")
