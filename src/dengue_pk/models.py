"""
Host--vector compartmental model for dengue transmission.

Dengue is not directly transmissible between humans: the virus passes from an
infectious human to a mosquito and back again. A single-population SIR model
therefore cannot represent it correctly, because it omits the delay imposed by
the mosquito's own infection cycle. The model below tracks both populations.

Humans (SEIR):
    dS_h/dt = -lambda_h S_h
    dE_h/dt =  lambda_h S_h - sigma_h E_h
    dI_h/dt =  sigma_h E_h  - gamma_h I_h
    dR_h/dt =  gamma_h I_h

Mosquitoes (SI, no recovery -- an infected mosquito stays infectious for life):
    dS_v/dt =  mu_v N_v - lambda_v S_v - mu_v S_v
    dI_v/dt =  lambda_v S_v - mu_v I_v

with forces of infection

    lambda_h(t) = beta(t) * m * I_v / N_v        (mosquito -> human)
    lambda_v(t) = beta(t) *     I_h / N_h        (human -> mosquito)

where m is the vector-to-host ratio and beta(t) is the climate-forced composite
transmission coefficient supplied by :mod:`dengue_pk.climate`.

The human exposed compartment E_h represents the intrinsic incubation period.
Omitting it biases the recovered transmission rate upward, because the model
would otherwise have to reproduce the observed epidemic growth without the delay
that genuinely slows it.

Incidence, not prevalence, is what surveillance observes, so the integrator also
accumulates new human infections; see :func:`simulate`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# State vector layout. Fractions of the respective populations, except CUM_INC
# which accumulates the human infection incidence used by the observation model.
S_H, E_H, I_H, R_H, S_V, I_V, CUM_INC = range(7)
N_STATES = 7


@dataclass(frozen=True)
class FixedParams:
    """Parameters held fixed from the literature rather than estimated.

    Human case data alone cannot identify mosquito lifespan or the incubation
    period; attempting to estimate them produces arbitrary values that trade off
    against the transmission rate. See docs/MODEL.md for sources.
    """

    gamma_h: float = 0.2            # human recovery rate, /day
    sigma_h: float = 0.2            # 1 / intrinsic incubation period, /day
    mu_v: float = 0.1               # mosquito mortality rate, /day
    vector_host_ratio: float = 2.0  # baseline mosquitoes per human

    @classmethod
    def from_config(cls, cfg: dict, key: str = "fixed") -> "FixedParams":
        """Build from configuration.

        ``key`` selects which literature parameter set to use. The alternative
        set (``fixed_alt``) sits at the other end of the same published ranges
        and is a factor in the robustness study: choosing within those ranges is
        an analytical degree of freedom like any other, and the study measures
        it rather than fixing it silently.
        """
        f = cfg["model"][key]
        return cls(gamma_h=f["gamma_h"], sigma_h=f["sigma_h"], mu_v=f["mu_v"],
                   vector_host_ratio=f["vector_host_ratio"])


def rhs(t: float, y: np.ndarray, beta_t, fixed: FixedParams) -> np.ndarray:
    """Right-hand side of the host--vector system.

    Parameters
    ----------
    t : float
        Time in days from the start of the study window.
    y : ndarray, shape (7,)
        State vector, fractions; see the module-level layout constants.
    beta_t : callable
        ``beta_t(t)`` returns the climate-forced transmission coefficient at
        time ``t``. Supplied as a callable rather than a constant so that the
        same model serves both the constant-beta null model and the
        climate-forced model.
    fixed : FixedParams
        Literature parameters.
    """
    s_h, e_h, i_h, _r_h, s_v, i_v = y[:6]
    beta = beta_t(t)
    m = fixed.vector_host_ratio

    # Forces of infection. Mosquito compartments are fractions of the mosquito
    # population, so i_v enters directly; the vector-to-host ratio converts the
    # mosquito prevalence into a per-human exposure.
    lambda_h = beta * m * i_v
    lambda_v = beta * i_h

    d = np.empty(N_STATES)
    d[S_H] = -lambda_h * s_h
    d[E_H] = lambda_h * s_h - fixed.sigma_h * e_h
    d[I_H] = fixed.sigma_h * e_h - fixed.gamma_h * i_h
    d[R_H] = fixed.gamma_h * i_h
    # Mosquito births replace deaths, holding the population constant; seasonal
    # abundance is carried by beta(t) instead, which keeps the state space small
    # and avoids a second poorly identified time-varying quantity.
    d[S_V] = fixed.mu_v * (1.0 - s_v) - lambda_v * s_v
    d[I_V] = lambda_v * s_v - fixed.mu_v * i_v
    # Incidence: the rate at which humans leave E_h and become infectious. This
    # is the quantity surveillance detects, up to the reporting fraction.
    d[CUM_INC] = fixed.sigma_h * e_h
    return d


def initial_state(i0_frac: float) -> np.ndarray:
    """Initial condition seeded with a small infected human fraction.

    The mosquito population starts fully susceptible: seeding both populations
    would introduce a second free parameter that the data cannot distinguish
    from the first.
    """
    y0 = np.zeros(N_STATES)
    y0[S_H] = 1.0 - i0_frac
    y0[E_H] = 0.0
    y0[I_H] = i0_frac
    y0[R_H] = 0.0
    y0[S_V] = 1.0
    y0[I_V] = 0.0
    y0[CUM_INC] = 0.0
    return y0


class IntegrationFailure(RuntimeError):
    """Raised when the state leaves the physically meaningful region.

    An explicit method has a finite stability region, and a sufficiently large
    transmission coefficient will drive the fixed step outside it. During
    multi-start fitting the optimiser does propose such values, so the failure
    must be signalled rather than allowed to propagate as NaN: a silent NaN
    turns into a meaningless residual, and the optimiser cannot distinguish an
    impossible parameter set from a merely poor one.
    """


def rk4_integrate(y0: np.ndarray, t_end: float, dt: float, beta_t,
                  fixed: FixedParams, check_every: int = 32, mu_v_t=None):
    """Fixed-step classical Runge--Kutta integration.

    A fixed step is used deliberately rather than an adaptive solver: the
    inference routines call this function many thousands of times, and a
    deterministic cost per call keeps timing comparisons between the classical
    and PINN approaches meaningful. Step-size adequacy is verified in
    ``tests/test_models.py`` by comparison against a tenfold finer step.

    The state is checked periodically, rather than every step, because the check
    is pure overhead on the overwhelming majority of calls that never diverge.
    Divergence is monotone once it begins, so a check every ``check_every``
    steps catches it while it is still representable.
    """
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    y = np.empty((n + 1, N_STATES))
    y[0] = y0

    # The stage evaluations below are the same RK4 as `rhs` expresses, written
    # out over scalars. Allocating four seven-element arrays per step dominated
    # the runtime — inference calls this tens of thousands of times, and the
    # global study calls it millions — and the arrays were pure overhead for a
    # system this small. `test_models.py` asserts the two forms agree to machine
    # precision, so `rhs` remains the readable definition and this the executed
    # one.
    gam, sig, mu_v, m = fixed.gamma_h, fixed.sigma_h, fixed.mu_v, \
        fixed.vector_host_ratio
    h6 = dt / 6.0
    h2 = dt / 2.0

    # beta(t) is needed at t_k, t_k + dt/2 and t_k + dt, so it is evaluated once
    # on the half-step grid and indexed thereafter. Calling the forcing inside
    # the loop cost two array interpolations per stage — twelve per step, and
    # more time than the differential equation itself.
    half_grid = np.arange(2 * n + 1) * h2
    if hasattr(beta_t, "on_grid"):
        beta_grid = np.asarray(beta_t.on_grid(half_grid), dtype=float)
    else:
        beta_grid = np.array([beta_t(x) for x in half_grid], dtype=float)
    # Converted to a Python list deliberately. Indexing a NumPy array returns a
    # np.float64, and every subsequent operation in the step then goes through
    # NumPy's scalar machinery rather than plain float arithmetic — which
    # doubled the cost of the loop when this was first written with array
    # indexing.
    beta_grid = beta_grid.tolist()

    # Mosquito mortality is constant unless a thermal response is supplied. It
    # is evaluated on the same half-step grid as beta, for the same reason.
    # Keeping the constant case on a filled grid rather than a branch inside the
    # step keeps one code path: `tests/test_models.py` asserts that a constant
    # mu_v_t reproduces the scalar path to machine precision.
    if mu_v_t is None:
        mu_grid = [mu_v] * (2 * n + 1)
    elif hasattr(mu_v_t, "on_grid"):
        mu_grid = np.asarray(mu_v_t.on_grid(half_grid), dtype=float).tolist()
    else:
        mu_grid = [float(mu_v_t(x)) for x in half_grid]

    def stages(beta, mu, s_h, e_h, i_h, s_v, i_v):
        lam_h = beta * m * i_v
        lam_v = beta * i_h
        return (-lam_h * s_h,
                lam_h * s_h - sig * e_h,
                sig * e_h - gam * i_h,
                gam * i_h,
                mu * (1.0 - s_v) - lam_v * s_v,
                lam_v * s_v - mu * i_v,
                sig * e_h)

    with np.errstate(over="ignore", invalid="ignore"):
        s_h, e_h, i_h, r_h, s_v, i_v, cum = (float(v) for v in y0)
        for k in range(n):
            b0 = beta_grid[2 * k]
            bh = beta_grid[2 * k + 1]
            b1 = beta_grid[2 * k + 2]
            m0 = mu_grid[2 * k]
            mh = mu_grid[2 * k + 1]
            m1 = mu_grid[2 * k + 2]
            a = stages(b0, m0, s_h, e_h, i_h, s_v, i_v)
            b = stages(bh, mh, s_h + h2 * a[0], e_h + h2 * a[1],
                       i_h + h2 * a[2], s_v + h2 * a[4], i_v + h2 * a[5])
            c = stages(bh, mh, s_h + h2 * b[0], e_h + h2 * b[1],
                       i_h + h2 * b[2], s_v + h2 * b[4], i_v + h2 * b[5])
            d = stages(b1, m1, s_h + dt * c[0], e_h + dt * c[1],
                       i_h + dt * c[2], s_v + dt * c[4], i_v + dt * c[5])

            s_h += h6 * (a[0] + 2 * b[0] + 2 * c[0] + d[0])
            e_h += h6 * (a[1] + 2 * b[1] + 2 * c[1] + d[1])
            i_h += h6 * (a[2] + 2 * b[2] + 2 * c[2] + d[2])
            r_h += h6 * (a[3] + 2 * b[3] + 2 * c[3] + d[3])
            s_v += h6 * (a[4] + 2 * b[4] + 2 * c[4] + d[4])
            i_v += h6 * (a[5] + 2 * b[5] + 2 * c[5] + d[5])
            cum += h6 * (a[6] + 2 * b[6] + 2 * c[6] + d[6])
            y[k + 1] = (s_h, e_h, i_h, r_h, s_v, i_v, cum)

            if (k % check_every) == 0:
                lo = min(s_h, e_h, i_h, r_h, s_v, i_v)
                hi = max(s_h, e_h, i_h, r_h, s_v, i_v)
                if not (lo == lo and hi == hi) or lo < -1e-6 or hi > 1.0 + 1e-6:
                    raise IntegrationFailure(
                        f"state left [0, 1] at t = {t[k + 1]:.1f} d")
    return t, y


def rk4_integrate_seir(y0_frac: float, t_end: float, dt: float, beta_t,
                       fixed: FixedParams, check_every: int = 32):
    """Human-only SEIR, the alternative model structure.

    Many published dengue analyses fit a directly transmitted SEIR or SIR model
    to case counts, omitting the vector entirely. That is a different structural
    assumption, not a simplification of the same one: it removes the delay the
    mosquito's own infection cycle imposes, so the same epidemic curve implies a
    different transmission rate.

    It is included so that the robustness study can vary model structure as well
    as the four fitting choices. If the instability were an artefact of the
    host--vector formulation, it would not survive here.

    State layout is (S, E, I, R, cumulative incidence); R0 = beta / gamma_h.
    """
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    y = np.empty((n + 1, 5))

    gam, sig = fixed.gamma_h, fixed.sigma_h
    h6, h2 = dt / 6.0, dt / 2.0

    half_grid = np.arange(2 * n + 1) * h2
    if hasattr(beta_t, "on_grid"):
        beta_grid = np.asarray(beta_t.on_grid(half_grid), dtype=float).tolist()
    else:
        beta_grid = [float(beta_t(x)) for x in half_grid]

    def stages(beta, s, e, i):
        lam = beta * i
        return (-lam * s, lam * s - sig * e, sig * e - gam * i, gam * i, sig * e)

    s, e, i, r, cum = 1.0 - y0_frac, 0.0, y0_frac, 0.0, 0.0
    y[0] = (s, e, i, r, cum)
    with np.errstate(over="ignore", invalid="ignore"):
        for k in range(n):
            b0, bh, b1 = beta_grid[2 * k], beta_grid[2 * k + 1], beta_grid[2 * k + 2]
            a = stages(b0, s, e, i)
            b = stages(bh, s + h2 * a[0], e + h2 * a[1], i + h2 * a[2])
            c = stages(bh, s + h2 * b[0], e + h2 * b[1], i + h2 * b[2])
            dd = stages(b1, s + dt * c[0], e + dt * c[1], i + dt * c[2])
            s += h6 * (a[0] + 2 * b[0] + 2 * c[0] + dd[0])
            e += h6 * (a[1] + 2 * b[1] + 2 * c[1] + dd[1])
            i += h6 * (a[2] + 2 * b[2] + 2 * c[2] + dd[2])
            r += h6 * (a[3] + 2 * b[3] + 2 * c[3] + dd[3])
            cum += h6 * (a[4] + 2 * b[4] + 2 * c[4] + dd[4])
            y[k + 1] = (s, e, i, r, cum)
            if (k % check_every) == 0:
                lo, hi = min(s, e, i, r), max(s, e, i, r)
                if not (lo == lo and hi == hi) or lo < -1e-6 or hi > 1.0 + 1e-6:
                    raise IntegrationFailure(
                        f"SEIR state left [0, 1] at t = {t[k + 1]:.1f} d")
    return t, y


def basic_reproduction_number(beta: float, fixed: FixedParams,
                              structure: str = "hostvector",
                              mu_v: float | None = None) -> float:
    """R0 for this parameterisation, from the next-generation argument.

    One infectious human infects mosquitoes at rate ``beta`` for ``1/gamma_h``
    days; one infectious mosquito infects humans at rate ``beta * m`` for
    ``1/mu_v`` days. R0 is the geometric mean of the two, since a full
    transmission cycle requires both steps:

        R0 = beta * sqrt( m / (gamma_h * mu_v) )

    With the default fixed parameters this is ``R0 = 10 * beta``, which is what
    makes an estimated ``beta`` interpretable and lets the plausible parameter
    range be stated in epidemiological rather than numerical terms.

    For the human-only SEIR structure the cycle has one step rather than two and
    the mosquito plays no part, giving the standard ``R0 = beta / gamma_h``, or
    ``5 * beta`` here. The two structures therefore map the same estimated beta
    to different R0 values, which is why every comparison across structures must
    be made on R0 and never on beta.
    """
    if structure == "seir":
        return beta / fixed.gamma_h
    # ``mu_v`` overrides the constant value for the temperature-dependent
    # structure, where mortality varies over the window and a single reported
    # R0 has to use its time average. Without the override that structure's R0
    # would be computed from a mortality the model never actually used.
    mu = fixed.mu_v if mu_v is None else float(mu_v)
    return beta * np.sqrt(fixed.vector_host_ratio / (fixed.gamma_h * mu))


def weekly_incidence(t: np.ndarray, y: np.ndarray, week_starts_days: np.ndarray,
                     population: float) -> np.ndarray:
    """Convert the cumulative incidence trajectory into weekly counts.

    Surveillance reports cases aggregated over a week, so the model must be
    compared against differences of the cumulative curve rather than against an
    instantaneous rate.

    Each observation covers the seven days following its own start date.
    Differencing consecutive entries instead is equivalent only while the
    observed weeks are contiguous, and silently assigns several weeks of
    incidence to a single reported week wherever they are not. The main study
    windows are contiguous by construction and the two forms agree there, but
    `scripts/09_aggregation_bias.py` intersects districts whose reporting weeks
    do not align, and there the distinction is real.
    """
    starts = np.asarray(week_starts_days, dtype=float)
    cum_start = np.interp(starts, t, y[:, CUM_INC])
    cum_end = np.interp(starts + 7.0, t, y[:, CUM_INC])
    return (cum_end - cum_start) * population


def simulate(theta: dict, beta_t, fixed: FixedParams, week_starts_days: np.ndarray,
             population: float, dt: float = 1.0, structure: str = "hostvector",
             mu_v_t=None):
    """Run the model and return expected reported cases per week.

    ``theta`` supplies the estimated parameters; only the initial infected
    fraction and the reporting fraction enter here, the transmission parameters
    having already been absorbed into ``beta_t``.

    The default step of one day is used because inference calls this function
    tens of thousands of times and the step dominates the total cost. Its
    adequacy under climate forcing — where the driver itself varies weekly — is
    verified in ``tests/test_models.py`` against an eightfold finer step, not
    assumed from the constant-coefficient case.
    """
    t_end = float(week_starts_days[-1] + 7.0)
    if structure == "seir":
        # The SEIR state vector is five wide and its cumulative incidence sits
        # in the last column, so it is padded to the host--vector layout before
        # `weekly_incidence`, which indexes that column by name.
        t, y5 = rk4_integrate_seir(theta["i0_frac"], t_end, dt, beta_t, fixed)
        y = np.zeros((len(t), N_STATES))
        y[:, [S_H, E_H, I_H, R_H]] = y5[:, :4]
        y[:, CUM_INC] = y5[:, 4]
    else:
        # ``mu_v_t`` is supplied only by the temperature-dependent-mortality
        # structure; passing None reproduces the constant-mortality path
        # exactly, which `tests/test_models.py` asserts.
        t, y = rk4_integrate(initial_state(theta["i0_frac"]), t_end, dt,
                             beta_t, fixed, mu_v_t=mu_v_t)
    incidence = weekly_incidence(t, y, week_starts_days, population)
    return theta["rho"] * incidence, (t, y)
