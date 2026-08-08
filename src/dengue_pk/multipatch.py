"""
Two-patch extension: does heterogeneity explain the aggregation bias?

Two findings so far both point at spatial heterogeneity as the cause. Fitting a
province's districts separately gives a higher transmission rate than fitting
their sum, and the observed early growth rate demands a higher R0 than any single
fitted value allows. In both cases the proposed explanation is the same — these
series are sums of outbreaks that peak at different times, and a homogeneous model
reads a broader curve as slower transmission.

That explanation has been asserted, not tested. This module tests it.

A two-patch model treats the aggregate as the sum of two independent epidemics,
each homogeneous within itself, with its own transmission rate, catchment size and
seeding time. Fitted to the aggregate alone — with no information about which
districts contributed or when — it should, if the explanation is right:

  1. fit the aggregate better than one patch, by enough to justify three extra
     parameters;
  2. recover transmission rates that bracket the separately fitted district
     values rather than sitting below both;
  3. reconcile the growth rate with the final size, since a fast early patch and
     a slower later one together produce a curve that rises quickly and still
     ends modestly.

If it fails these, the explanation is wrong and the paper must say so. That is
the point of running it.

The patches do not interact. Coupling them would add parameters that a single
aggregate series certainly cannot identify, and independence is the conservative
choice: it is the least favourable version of the hypothesis.
"""

from __future__ import annotations

import numpy as np

from .climate import ConstantForcing
from .models import FixedParams, initial_state, rk4_integrate, weekly_incidence

# Per-patch parameters. The offset shifts a patch's epidemic in time, which is
# what makes the two asynchronous — the feature the whole hypothesis rests on.
PATCH_PARAMS = ("beta", "pop_frac", "i0_frac", "offset_days")


def patch_names(n_patches: int) -> list[str]:
    names = []
    for k in range(n_patches):
        for p in PATCH_PARAMS:
            if k == 0 and p == "offset_days":
                continue          # the first patch defines the time origin
            names.append(f"{p}_{k}")
    return names


def predict_multipatch(theta: dict, days: np.ndarray, population: float,
                       fixed: FixedParams, rho: float, n_patches: int,
                       dt: float = 1.0) -> np.ndarray:
    """Expected weekly reported cases summed over independent patches.

    Each patch is integrated on its own shifted clock. A patch seeded later
    contributes its rise to the aggregate later, which is precisely the
    asynchrony that a single-patch fit cannot represent.
    """
    total = np.zeros(len(days))
    for k in range(n_patches):
        beta = theta[f"beta_{k}"]
        pop = theta[f"pop_frac_{k}"] * population
        i0 = theta[f"i0_frac_{k}"]
        offset = theta.get(f"offset_days_{k}", 0.0)

        # Shift the observation times rather than the model: a positive offset
        # means the patch started later, so at a given calendar week it is
        # earlier in its own epidemic.
        local_days = days - offset
        t_end = float(max(local_days[-1] + 7.0, 7.0))
        if t_end <= 0:
            continue
        t, y = rk4_integrate(initial_state(i0), t_end, dt,
                             ConstantForcing(beta), fixed)
        # Weeks before this patch was seeded contribute nothing.
        weekly = weekly_incidence(t, y, np.maximum(local_days, 0.0), pop)
        weekly = np.where(local_days < 0.0, 0.0, weekly)
        total += rho * weekly
    return np.maximum(total, 1e-9)
