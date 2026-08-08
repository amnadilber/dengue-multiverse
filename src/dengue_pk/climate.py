"""
Climate forcing of the transmission coefficient.

Dengue transmission in Pakistan is strongly seasonal, peaking after the monsoon.
Two mechanisms drive this and both are represented here:

* **Temperature** governs mosquito biting rate, development and the extrinsic
  incubation period. Transmission is negligible below roughly 18 C and falls
  again above roughly 34 C, so its effect is not monotonic.
* **Rainfall** creates breeding habitat, but with a delay: eggs laid in standing
  water need one to three weeks to become biting adults. Rainfall therefore
  enters lagged and smoothed, not contemporaneously.

The composite transmission coefficient is modelled in log space,

    beta(t) = beta_0 * exp( a_temp * z_T(t) + a_rain * z_P(t) )

with ``z_T`` and ``z_P`` standardised climate covariates. The log-linear form
keeps beta positive for any coefficient values, and standardising the covariates
puts the two coefficients on a comparable scale so that their relative magnitudes
are interpretable.

A constant-beta null model is provided as well: it is the honest baseline against
which any claimed benefit of climate forcing must be judged.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd


def lagged_smoothed_rain(daily_rain: pd.Series, lag_weeks: int,
                         smooth_weeks: int) -> pd.Series:
    """Apply the mosquito development delay to a daily rainfall series.

    The rainfall relevant to this week's biting mosquitoes fell some weeks ago,
    and it is the accumulation over a period rather than a single day's total
    that matters for habitat. Both operations are applied before any
    standardisation so that the lag is expressed in real time units.
    """
    smoothed = daily_rain.rolling(window=smooth_weeks * 7, min_periods=1).mean()
    return smoothed.shift(lag_weeks * 7).bfill()


def standardise(x: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Centre and scale, returning the constants so the transform is invertible.

    Standardisation constants are computed once on the fitting window and then
    reused; recomputing them on a forecast window would leak information from
    the period being predicted.
    """
    mu, sd = float(np.mean(x)), float(np.std(x))
    sd = sd if sd > 0 else 1.0
    return (x - mu) / sd, mu, sd


def briere(temp_c, t_min: float, t_max: float):
    """Brière temperature response, normalised to peak at one.

        B(T) = T (T - T_min) sqrt(T_max - T),   T_min < T < T_max,  else 0

    Mosquito traits — biting rate, development, survival — rise from a lower
    thermal threshold, peak, and collapse before an upper one. The Brière form is
    the standard empirical description of that shape and is what the entomological
    literature reports thresholds for.

    Using it in place of a free log-linear coefficient matters here. Over a range
    that spans both sides of the thermal optimum, a monotone term must choose a
    single direction, and fitting one to a season whose temperature falls while
    the epidemic grows returns a negative coefficient — an artefact of covariate
    shape rather than a biological effect. A unimodal response cannot produce that
    artefact, because it does not have a direction to choose.

    Normalising to a peak of one keeps the interpretation of the baseline
    transmission coefficient it multiplies: beta_0 is transmission at the thermal
    optimum.
    """
    t = np.asarray(temp_c, dtype=float)
    inside = (t > t_min) & (t < t_max)
    out = np.zeros_like(t)
    out[inside] = (t[inside] * (t[inside] - t_min)
                   * np.sqrt(np.maximum(t_max - t[inside], 0.0)))
    peak = out.max() if out.max() > 0 else 1.0
    # Normalise by the analytic maximum over the admissible range, not by the
    # maximum observed in this particular window: normalising per window would
    # make beta_0 mean something different in each one and destroy comparability.
    #
    # That maximum depends only on the thermal limits, so it is memoised. It was
    # previously recomputed on a 2001-point grid inside every call, and this
    # function is called once per residual evaluation — tens of millions of times
    # across the global study, for a constant.
    analytic_peak = _briere_peak(t_min, t_max)
    return out / analytic_peak if analytic_peak > 0 else out / peak


@lru_cache(maxsize=32)
def _briere_peak(t_min: float, t_max: float) -> float:
    grid = np.linspace(t_min + 1e-6, t_max - 1e-6, 4001)
    return float(np.max(grid * (grid - t_min) * np.sqrt(t_max - grid)))


class MechanisticForcing:
    """Transmission driven by a unimodal temperature response and lagged rainfall.

        beta(t) = beta_0 * B_norm(T(t)) ** a_temp * exp( a_rain * z_rain(t) )

    The temperature exponent ``a_temp`` makes the temperature term a nested
    hypothesis rather than an assumption: at ``a_temp = 0`` temperature drops out
    entirely, at ``a_temp = 1`` the literature response applies in full, and
    values between interpolate. Whether temperature belongs in the model is then
    something the fit can answer, instead of something the model structure
    imposes.

    Rainfall keeps its log-linear form. Unlike temperature it has a defensible
    monotone effect over the observed range — more standing water means more
    breeding habitat — so a unimodal response would be unwarranted there.
    """

    def __init__(self, days, temp_c, z_rain, beta_0, a_temp, a_rain,
                 t_min: float, t_max: float):
        self._days = np.asarray(days, dtype=float)
        self._b_temp = briere(temp_c, t_min, t_max)
        self._z_rain = np.asarray(z_rain, dtype=float)
        self.beta_0 = float(beta_0)
        self.a_temp = float(a_temp)
        self.a_rain = float(a_rain)

    def _beta(self, b, zr):
        # A floor keeps the power finite where the Brière response is zero, i.e.
        # outside the thermal limits, without letting transmission reach zero and
        # stall the epidemic irrecoverably.
        b = np.maximum(b, 1e-6)
        return self.beta_0 * b ** self.a_temp * np.exp(self.a_rain * zr)

    def __call__(self, t: float) -> float:
        b = np.interp(t, self._days, self._b_temp)
        zr = np.interp(t, self._days, self._z_rain)
        return float(self._beta(b, zr))

    def on_grid(self, days=None) -> np.ndarray:
        d = self._days if days is None else np.asarray(days, dtype=float)
        b = np.interp(d, self._days, self._b_temp)
        zr = np.interp(d, self._days, self._z_rain)
        return self._beta(b, zr)


class ClimateForcing:
    """Callable ``beta(t)`` built from standardised climate covariates.

    Parameters
    ----------
    days : ndarray
        Days from the start of the window at which the covariates are given.
    z_temp, z_rain : ndarray
        Standardised temperature and lagged rainfall on the same grid.
    beta_0, a_temp, a_rain : float
        Baseline transmission and the two forcing coefficients.

    Values between grid points are linearly interpolated. Outside the grid the
    endpoint values are held, which matters only for the final partial week.
    """

    def __init__(self, days, z_temp, z_rain, beta_0, a_temp, a_rain):
        self._days = np.asarray(days, dtype=float)
        self._z_temp = np.asarray(z_temp, dtype=float)
        self._z_rain = np.asarray(z_rain, dtype=float)
        self.beta_0 = float(beta_0)
        self.a_temp = float(a_temp)
        self.a_rain = float(a_rain)

    def __call__(self, t: float) -> float:
        zt = np.interp(t, self._days, self._z_temp)
        zr = np.interp(t, self._days, self._z_rain)
        return self.beta_0 * np.exp(self.a_temp * zt + self.a_rain * zr)

    def on_grid(self, days=None) -> np.ndarray:
        """Vectorised evaluation, for plotting the recovered beta(t)."""
        d = self._days if days is None else np.asarray(days, dtype=float)
        zt = np.interp(d, self._days, self._z_temp)
        zr = np.interp(d, self._days, self._z_rain)
        return self.beta_0 * np.exp(self.a_temp * zt + self.a_rain * zr)


class ConstantForcing:
    """Null model: transmission does not vary with climate.

    Retained deliberately. If the climate-forced model cannot beat this on
    held-out data, the climate terms are not earning their parameters and the
    paper must say so.
    """

    def __init__(self, beta_0: float):
        self.beta_0 = float(beta_0)

    def __call__(self, t: float) -> float:
        return self.beta_0

    def on_grid(self, days) -> np.ndarray:
        return np.full(len(np.asarray(days)), self.beta_0)


# Aedes aegypti adult lifespan thermal limits, from the trait-based literature
# (Mordecai et al. 2017 and the mechanistic dengue-transmission models that
# follow it). Lifespan is unimodal in temperature with a broad optimum near
# 24 C, falling to zero at both ends of this range.
LIFESPAN_T_MIN = 11.7
LIFESPAN_T_MAX = 37.2

#: Lower bound on the normalised lifespan response. Outside the thermal limits
#: the quadratic is zero, which would make mortality infinite and the epidemic
#: stop dead. The floor caps mortality at 1/FLOOR times its optimum value.
_LIFESPAN_FLOOR = 0.15


class ThermalMortality:
    """Mosquito mortality that rises as temperature departs from the optimum.

        mu_v(t) = mu_v_base / max( L_norm(T(t)), floor )

    where ``L_norm`` is a quadratic lifespan response normalised to 1 at its
    peak, so that at the thermal optimum mortality equals the constant value the
    simpler structure uses and the two models coincide there.

    This exists because the criticism a reader will make of the host--vector
    model in this study is that it is too simple to represent what applied
    dengue modelling actually does: holding mosquito lifespan constant across a
    season whose temperature swings by fifteen degrees is not what the
    mechanistic literature does. Making mortality temperature-dependent is the
    standard realism upgrade, and including it as a third structure lets the
    robustness result be tested against a model closer to current practice
    rather than asserted to hold there.

    Note what this does *not* do: it adds no estimated parameter. The thermal
    limits are fixed from vector biology exactly as the constant mortality was
    fixed from vector biology, so the comparison is between structures rather
    than between one model with more freedom than another.
    """

    def __init__(self, days, temp_c, mu_v_base: float,
                 t_min: float = LIFESPAN_T_MIN, t_max: float = LIFESPAN_T_MAX):
        self._days = np.asarray(days, dtype=float)
        temp = np.asarray(temp_c, dtype=float)
        self.mu_v_base = float(mu_v_base)
        peak = 0.5 * (t_min + t_max)
        q_peak = (peak - t_min) * (t_max - peak)
        q = np.where((temp > t_min) & (temp < t_max),
                     (temp - t_min) * (t_max - temp), 0.0) / q_peak
        self._mu = self.mu_v_base / np.maximum(q, _LIFESPAN_FLOOR)

    def __call__(self, t: float) -> float:
        return float(np.interp(t, self._days, self._mu))

    def on_grid(self, days=None) -> np.ndarray:
        d = self._days if days is None else np.asarray(days, dtype=float)
        return np.interp(d, self._days, self._mu)

    def mean(self) -> float:
        """Time-averaged mortality, used where a single R0 must be reported."""
        return float(np.mean(self._mu))
