"""
Classical inversion: recovering transmission parameters by nonlinear least squares
over repeated forward solves.

Three choices here are not incidental and are argued rather than assumed.

**The residual scale.** Fitting raw count differences would let the epidemic peak
dominate the objective: a week with 2,000 cases contributes a hundred times the
squared error of a week with 200, so the fit would be driven almost entirely by a
handful of peak weeks and would ignore the growth phase, which is precisely where
the transmission rate is most identifiable. Poisson deviance residuals are used
instead. Their sum of squares is the Poisson deviance, so minimising it is
equivalent to Poisson maximum likelihood, and each week contributes in proportion
to how surprising it is rather than how large it is.

**Unconstrained parameterisation.** Rates and fractions must stay positive, and a
reporting fraction must additionally stay below one. Rather than imposing box
constraints on the optimiser, parameters are transformed — log for positive
quantities, logit for the reporting fraction — so the search is unconstrained and
cannot stall against a boundary.

**Multi-start.** The objective is not convex, and a single descent from one guess
can settle in a local minimum that looks converged. Several starts are drawn from
the configured ranges and the best result is kept; the spread across starts is
reported, because agreement between independent starts is evidence the optimum is
global while disagreement is evidence it is not.

Overdispersion is expected in surveillance counts and is not assumed away: the
dispersion is estimated after fitting and used to widen the confidence intervals,
which would otherwise be too narrow to believe.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammaln

from .climate import (ClimateForcing, ConstantForcing, MechanisticForcing,
                      ThermalMortality, standardise)
from .models import FixedParams, IntegrationFailure, simulate

# Parameters estimated for each model. The null model omits the climate terms,
# which is what makes the two models nested.
#
# The reporting fraction is absent from both: it is perfectly confounded with the
# population at risk (see docs/MODEL.md), so it is fixed from the literature and
# the identifiable product is attributed to `pop_frac`.
CLIMATE_PARAMS = ["beta_0", "a_temp", "a_rain", "pop_frac", "i0_frac"]
NULL_PARAMS = ["beta_0", "pop_frac", "i0_frac"]

# Transform used per parameter to make the search unconstrained.
#
# `a_temp` is log-transformed, which keeps the Brière exponent non-negative —
# the property that prevents the wrong-signed artefact. Under the log-linear
# form the same parameter is a free slope that must be allowed to go negative,
# so the transform is switched by `set_temperature_form`.
_TRANSFORM = {"beta_0": "log", "pop_frac": "logit", "i0_frac": "log",
              "a_temp": "log", "a_rain": "identity"}


def set_temperature_form(form: str) -> None:
    """Match the transform on ``a_temp`` to the temperature parameterisation.

    Under `loglinear` the coefficient is a slope on standardised temperature and
    must be free to take either sign; forcing it positive would prevent the very
    artefact this comparison exists to demonstrate.
    """
    _TRANSFORM["a_temp"] = "identity" if form == "loglinear" else "log"


def _to_unconstrained(name: str, value: float) -> float:
    kind = _TRANSFORM[name]
    if kind == "log":
        return float(np.log(value))
    if kind == "logit":
        v = min(max(value, 1e-9), 1 - 1e-9)
        return float(np.log(v / (1 - v)))
    return float(value)


def _from_unconstrained(name: str, x: float) -> float:
    kind = _TRANSFORM[name]
    if kind == "log":
        return float(np.exp(np.clip(x, -30, 30)))
    if kind == "logit":
        return float(1.0 / (1.0 + np.exp(-np.clip(x, -30, 30))))
    return float(x)


def pack(theta: dict, names) -> np.ndarray:
    return np.array([_to_unconstrained(n, theta[n]) for n in names])


def unpack(x: np.ndarray, names) -> dict:
    return {n: _from_unconstrained(n, xi) for n, xi in zip(names, x)}


@dataclass
class Dataset:
    """One study window, prepared for fitting."""

    days: np.ndarray            # days from window start, one per week
    cases: np.ndarray           # reported cases per week
    temp_c: np.ndarray          # weekly mean temperature, degrees C
    z_rain: np.ndarray          # standardised lagged rainfall
    population: float           # census population of the region
    label: str = ""

    @classmethod
    def from_frame(cls, df, label: str = "") -> "Dataset":
        z_r, *_ = standardise(df["rain_lagged"].to_numpy())
        return cls(days=df["days_from_start"].to_numpy(float),
                   cases=df["cases"].to_numpy(float),
                   temp_c=df["T2M"].to_numpy(float), z_rain=z_r,
                   population=float(df["population"].iloc[0]),
                   label=label)

    def head(self, n_weeks: int) -> "Dataset":
        """The first ``n_weeks`` of the window, for chronological train/test splits.

        Splits must be by time, never at random: predicting a withheld week from
        weeks on either side of it is interpolation rather than forecasting, and
        would flatter every model tested.

        Note that the standardisation constants in ``z_temp`` and ``z_rain`` were
        computed on the full window. That is a deliberate simplification — the
        covariates are climate, which is not itself being predicted — but it does
        mean the held-out comparison is not a strict out-of-sample test of the
        standardisation, only of the fitted parameters.
        """
        return Dataset(self.days[:n_weeks], self.cases[:n_weeks],
                       self.temp_c[:n_weeks], self.z_rain[:n_weeks],
                       self.population, f"{self.label} (first {n_weeks} weeks)")


def make_forcing(theta: dict, data: Dataset, model: str, cfg: dict):
    """Build beta(t) for the requested model and temperature parameterisation.

    Two forms of the temperature term are supported, because which one is chosen
    turned out to determine the sign of the estimated effect:

    * ``briere`` — a unimodal response with thermal limits fixed from vector
      biology, raised to an estimated exponent. Cannot represent "transmission
      falls as conditions approach the optimum".
    * ``loglinear`` — a free coefficient on standardised temperature. Monotone,
      so fitted to a season whose temperature falls while cases rise it must
      return a negative coefficient whatever the biology.

    The second is the conventional choice in the applied literature and is
    retained so that the consequence of making it can be measured rather than
    asserted.
    """
    if model != "climate":
        return ConstantForcing(theta["beta_0"])

    form = cfg["model"].get("temperature_form", "briere")
    if form == "loglinear":
        z_temp, *_ = standardise(data.temp_c)
        return ClimateForcing(data.days, z_temp, data.z_rain,
                              theta["beta_0"], theta["a_temp"], theta["a_rain"])
    tr = cfg["model"]["temperature_response"]
    return MechanisticForcing(data.days, data.temp_c, data.z_rain,
                              theta["beta_0"], theta["a_temp"],
                              theta["a_rain"], tr["t_min_c"], tr["t_max_c"])


def predict(theta: dict, data: Dataset, fixed: FixedParams, model: str,
            cfg: dict, days=None) -> np.ndarray:
    """Expected reported cases per week.

    The population at risk is ``pop_frac`` times the census population, and the
    reporting fraction is fixed from configuration — the two are confounded, so
    only one of them can be estimated.

    Propagates :class:`IntegrationFailure` unchanged; callers inside the
    optimiser convert it into a finite penalty, while callers reporting results
    should let it surface rather than silently substitute a fitted curve that
    was never computed.
    """
    forcing = make_forcing(theta, data, model, cfg)
    week_days = data.days if days is None else days
    at_risk = theta["pop_frac"] * data.population
    sim_theta = {"i0_frac": theta["i0_frac"],
                 "rho": cfg["model"]["fixed"]["rho_fixed"]}
    structure = cfg["model"].get("structure", "hostvector")
    # The third structure differs from the first only in letting mosquito
    # mortality follow temperature; the thermal limits are fixed from vector
    # biology and it estimates no additional parameter.
    mu_v_t = (ThermalMortality(data.days, data.temp_c, fixed.mu_v)
              if structure == "hostvector_tempmort" else None)
    mu, _ = simulate(sim_theta, forcing, fixed, week_days, at_risk,
                     structure=structure, mu_v_t=mu_v_t)
    return np.maximum(mu, 1e-9)


def poisson_deviance_residuals(y: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """Signed square roots of the per-observation Poisson deviance.

    Their sum of squares is the deviance, so least_squares on these residuals
    performs Poisson maximum likelihood.

    The floor against division by zero is applied *only inside the logarithm*.
    Flooring ``mu`` everywhere would corrupt the quantity it was meant to
    protect: an observation of zero predicted as zero would acquire a spurious
    deviance of ``2 x floor``, and with a floor of 1e-9 that surfaces as a
    residual of about -4e-5 where the true value is exactly zero. Small, but it
    is a systematic bias applied to every zero week, and the windows here begin
    and end with many of them.
    """
    y = np.asarray(y, dtype=float)
    mu = np.asarray(mu, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(y > 0, y * np.log(y / np.maximum(mu, 1e-300)), 0.0)
    d = 2.0 * (term - (y - mu))
    d = np.maximum(d, 0.0)                 # guard against tiny negatives
    return np.sign(y - mu) * np.sqrt(d)


def poisson_loglik(y: np.ndarray, mu: np.ndarray) -> float:
    mu = np.maximum(mu, 1e-9)
    return float(np.sum(y * np.log(mu) - mu - gammaln(y + 1.0)))


# ---------------------------------------------------------------------------
# Negative binomial observation model
#
# Surveillance counts here are dispersed far beyond Poisson — the national fit
# gives a Pearson dispersion of 164, meaning the variance is over a hundred
# times the mean. Under Poisson the epidemic peak is treated as almost
# noise-free, so the fit is dragged toward reproducing a handful of peak weeks
# at the expense of everything else, and the reported precision is fiction.
#
# The NB2 parameterisation has variance mu + mu^2/k, so k controls how much
# extra variability is permitted; k -> infinity recovers Poisson exactly.
# ---------------------------------------------------------------------------
def nb_loglik(y: np.ndarray, mu: np.ndarray, k: float) -> float:
    mu = np.maximum(mu, 1e-9)
    return float(np.sum(
        gammaln(y + k) - gammaln(k) - gammaln(y + 1.0)
        + k * np.log(k / (k + mu)) + y * np.log(mu / (k + mu))))


def nb_deviance_residuals(y: np.ndarray, mu: np.ndarray, k: float) -> np.ndarray:
    """Signed square roots of the per-observation negative binomial deviance."""
    y = np.asarray(y, dtype=float)
    mu = np.asarray(mu, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        term1 = np.where(y > 0, y * np.log(y / np.maximum(mu, 1e-300)), 0.0)
    term2 = (y + k) * np.log((y + k) / (mu + k))
    d = np.maximum(2.0 * (term1 - term2), 0.0)
    return np.sign(y - mu) * np.sqrt(d)


def estimate_dispersion_k(y: np.ndarray, mu: np.ndarray,
                          bounds=(1e-3, 1e6)) -> float:
    """Maximum likelihood estimate of the NB dispersion, given fitted means.

    Estimated by profiling rather than jointly with the transmission parameters.
    The deviance that the optimiser minimises is itself a function of k, so
    minimising both at once would not maximise the likelihood; alternating
    between the two is the standard treatment and is what `glm.nb` does.
    """
    from scipy.optimize import minimize_scalar

    res = minimize_scalar(lambda lk: -nb_loglik(y, mu, float(np.exp(lk))),
                          bounds=(np.log(bounds[0]), np.log(bounds[1])),
                          method="bounded")
    return float(np.exp(res.x))


@dataclass
class FitResult:
    theta: dict
    model: str
    names: list
    loglik: float
    aic: float
    dispersion: float
    observation: str = "poisson"
    nb_k: float | None = None
    stderr: dict = field(default_factory=dict)
    corr: np.ndarray | None = None
    n_starts: int = 0
    n_converged: int = 0
    start_spread: dict = field(default_factory=dict)
    n_forward_solves: int = 0
    seconds: float = 0.0

    def summary(self) -> str:
        obs = (f"NB(k={self.nb_k:.2f})" if self.observation == "nb"
               else "quasi-Poisson")
        lines = [f"{self.model} model [{obs}] — logL {self.loglik:,.1f}, "
                 f"AIC {self.aic:,.1f}, dispersion {self.dispersion:.1f}"]
        for n in self.names:
            se = self.stderr.get(n)
            se_txt = f" ± {se:.4g}" if se is not None and np.isfinite(se) else ""
            lines.append(f"    {n:<10} {self.theta[n]:.6g}{se_txt}")
        lines.append(f"    starts converged to the best optimum: "
                     f"{self.n_converged}/{self.n_starts}")
        return "\n".join(lines)


def fit(data: Dataset, cfg: dict, fixed: FixedParams, model: str = "climate",
        seed: int | None = None, observation: str = "poisson",
        max_nb_rounds: int = 8, nb_tol: float = 1e-3,
        start_from: dict | None = None, n_starts_override: int | None = None,
        fixed_nb_k: float | None = None) -> FitResult:
    """Multi-start Levenberg--Marquardt fit of one model to one window.

    With ``observation="nb"`` the dispersion is profiled: the parameters are
    fitted given a value of k, k is then re-estimated by maximum likelihood from
    the fitted means, and the two alternate until k stops moving. Fitting both
    at once would not maximise the likelihood, because the deviance being
    minimised is itself a function of k.

    ``start_from`` seeds the first start at a known parameter set instead of the
    configured initial values, and ``n_starts_override`` reduces the number of
    starts. Both exist for the bootstrap, where each replicate differs only by
    resampling noise and the point estimate is already a good starting guess;
    running the full multi-start search on every one of several hundred
    replicates would cost hours for no gain. ``fixed_nb_k`` holds the dispersion
    at its point-estimate value, so replicates are not each re-profiling it.
    """
    import time

    names = CLIMATE_PARAMS if model == "climate" else NULL_PARAMS
    est = cfg["model"]["estimated"]
    n_starts = (n_starts_override if n_starts_override is not None
                else cfg["inference"]["classical"]["n_restarts"])
    rng = np.random.default_rng(cfg["seed"] if seed is None else seed)

    counter = {"n": 0}
    # Starting value for the dispersion. Very large is effectively Poisson, so
    # the first round reproduces the Poisson fit and subsequent rounds move away
    # from it only as far as the data require.
    nb_k = {"value": 1e6 if fixed_nb_k is None else float(fixed_nb_k)}

    # Penalty returned for parameter sets that make the model diverge. It must be
    # large enough to be rejected but finite and constant, so that the optimiser
    # sees a flat barrier rather than a cliff of arbitrary values it could try to
    # descend.
    PENALTY = 1e3

    def residual(x: np.ndarray) -> np.ndarray:
        counter["n"] += 1
        theta = unpack(x, names)
        try:
            mu = predict(theta, data, fixed, model, cfg)
        except IntegrationFailure:
            return np.full(len(data.cases), PENALTY)
        r = (nb_deviance_residuals(data.cases, mu, nb_k["value"])
             if observation == "nb"
             else poisson_deviance_residuals(data.cases, mu))
        if not np.all(np.isfinite(r)):
            return np.full(len(data.cases), PENALTY)
        return r

    starts = []
    for i in range(n_starts):
        theta0 = {}
        for n in names:
            lo, hi, init = est[n]["min"], est[n]["max"], est[n]["init"]
            if i == 0 and start_from is not None:
                theta0[n] = start_from[n]              # warm start
            elif i == 0:
                theta0[n] = init                       # configured starting point
            elif _TRANSFORM[n] in ("log", "logit"):    # sample log-uniformly
                if lo <= 0:
                    raise ValueError(
                        f"config: '{n}' is sampled log-uniformly, so its 'min' "
                        f"must be strictly positive (got {lo}). Set a small "
                        f"positive bound; the optimiser is unconstrained and can "
                        f"still search below it.")
                theta0[n] = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            else:
                theta0[n] = float(rng.uniform(lo, hi))
        starts.append(theta0)

    def run_all_starts():
        sols = []
        for theta0 in starts:
            try:
                sols.append(least_squares(residual, pack(theta0, names),
                                          method="lm", max_nfev=4000))
            except Exception:                          # a start that diverges
                continue
        return sols

    t0 = time.time()
    solutions = run_all_starts()

    if observation == "nb" and fixed_nb_k is None:
        # Alternate: refit given k, re-estimate k from the fitted means, repeat.
        for _ in range(max_nb_rounds):
            live = [s for s in solutions
                    if float(np.sum(s.fun ** 2))
                    < 0.99 * len(data.cases) * PENALTY ** 2]
            if not live:
                break
            best_so_far = min(live, key=lambda s: float(np.sum(s.fun ** 2)))
            try:
                mu_hat = predict(unpack(best_so_far.x, names), data, fixed,
                                 model, cfg)
            except IntegrationFailure:
                break
            k_new = estimate_dispersion_k(data.cases, mu_hat)
            if abs(np.log(k_new) - np.log(nb_k["value"])) < nb_tol:
                nb_k["value"] = k_new
                break
            nb_k["value"] = k_new
            solutions = run_all_starts()

    seconds = time.time() - t0

    # Discard starts that ended on the divergence barrier: their cost is the
    # penalty, not a fit, and including them would corrupt both the choice of
    # best solution and the reported spread across starts.
    diverged_cost = len(data.cases) * PENALTY ** 2
    solutions = [s for s in solutions
                 if float(np.sum(s.fun ** 2)) < 0.99 * diverged_cost]
    if not solutions:
        raise RuntimeError(
            f"every start diverged for the {model} model; the sampling ranges in "
            f"config are too wide for this window")

    best = min(solutions, key=lambda s: float(np.sum(s.fun ** 2)))
    theta = unpack(best.x, names)
    mu = predict(theta, data, fixed, model, cfg)

    if observation == "nb":
        # Re-estimate the dispersion at the final parameter values, and count it
        # as a parameter in AIC: it was estimated from the same data. Held fixed
        # when the caller supplied it, as the bootstrap does.
        if fixed_nb_k is None:
            nb_k["value"] = estimate_dispersion_k(data.cases, mu)
        ll = nb_loglik(data.cases, mu, nb_k["value"])
        k = len(names) + 1
    else:
        ll = poisson_loglik(data.cases, mu)
        k = len(names)
    aic = -2.0 * ll + 2.0 * k

    # Overdispersion relative to the assumed observation model. Under Poisson a
    # value far above one means the assumption understates variability, so
    # standard errors computed under it would be too narrow and are inflated
    # below. Under the negative binomial the variance already includes the
    # mu^2/k term, so a value near one is the sign that the dispersion has been
    # accounted for rather than merely absorbed.
    var = (np.maximum(mu, 1e-9) if observation != "nb"
           else np.maximum(mu, 1e-9) + mu ** 2 / nb_k["value"])
    resid_pearson = (data.cases - mu) / np.sqrt(var)
    dof = max(len(data.cases) - k, 1)
    dispersion = float(np.sum(resid_pearson ** 2) / dof)

    # Asymptotic standard errors from the Jacobian, inflated by the dispersion.
    stderr, corr = {}, None
    try:
        J = best.jac
        cov = np.linalg.inv(J.T @ J) * dispersion
        se_unc = np.sqrt(np.diag(cov))
        sd = np.sqrt(np.diag(cov))
        corr = cov / np.outer(sd, sd)
        # Delta method back through the transform.
        for n, xi, se in zip(names, best.x, se_unc):
            kind = _TRANSFORM[n]
            if kind == "log":
                stderr[n] = float(np.exp(xi) * se)
            elif kind == "logit":
                p = 1.0 / (1.0 + np.exp(-xi))
                stderr[n] = float(p * (1 - p) * se)
            else:
                stderr[n] = float(se)
    except np.linalg.LinAlgError:
        # A singular Jacobian means the parameters are not locally identifiable.
        stderr = {n: np.inf for n in names}

    best_cost = float(np.sum(best.fun ** 2))
    n_conv = sum(1 for s in solutions
                 if abs(float(np.sum(s.fun ** 2)) - best_cost)
                 <= 1e-3 * max(best_cost, 1.0))
    spread = {n: float(np.std([unpack(s.x, names)[n] for s in solutions]))
              for n in names}

    return FitResult(theta=theta, model=model, names=list(names), loglik=ll,
                     aic=aic, dispersion=dispersion,
                     observation=observation,
                     nb_k=nb_k["value"] if observation == "nb" else None,
                     stderr=stderr, corr=corr,
                     n_starts=len(starts), n_converged=n_conv,
                     start_spread=spread, n_forward_solves=counter["n"],
                     seconds=seconds)
