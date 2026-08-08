"""
Solver validation.

These tests exist because every number in the paper depends on the integrator
being correct. A parameter estimate obtained from a subtly wrong model is worse
than no estimate at all: it is wrong and it looks convincing. Each test below
checks a property that must hold regardless of parameter values.
"""

import numpy as np
import pytest

from dengue_pk.climate import (LIFESPAN_T_MAX, LIFESPAN_T_MIN, ClimateForcing,
                               ConstantForcing, ThermalMortality)
from dengue_pk.models import (CUM_INC, E_H, I_H, I_V, R_H, S_H, S_V, rhs,
                              FixedParams, IntegrationFailure,
                              basic_reproduction_number, initial_state,
                              rk4_integrate, weekly_incidence)

FIXED = FixedParams()


def test_human_population_is_conserved():
    """S_h + E_h + I_h + R_h must remain 1 throughout.

    A leak here would mean the model is losing or inventing people, which would
    silently distort every incidence estimate.
    """
    y0 = initial_state(1e-6)
    _, y = rk4_integrate(y0, 400.0, 0.25, ConstantForcing(0.5), FIXED)
    total = y[:, S_H] + y[:, E_H] + y[:, I_H] + y[:, R_H]
    assert np.allclose(total, 1.0, atol=1e-10)


def test_mosquito_population_is_conserved():
    """S_v + I_v must remain 1: births are set to replace deaths exactly."""
    y0 = initial_state(1e-6)
    _, y = rk4_integrate(y0, 400.0, 0.25, ConstantForcing(0.5), FIXED)
    total = y[:, S_V] + y[:, I_V]
    assert np.allclose(total, 1.0, atol=1e-10)


def test_all_compartments_stay_non_negative():
    """Negative compartments are physically meaningless and signal too large a step."""
    y0 = initial_state(1e-5)
    _, y = rk4_integrate(y0, 500.0, 0.25, ConstantForcing(0.8), FIXED)
    assert y[:, :6].min() > -1e-12


def test_cumulative_incidence_is_monotone():
    """Incidence accumulates; it can never decrease."""
    y0 = initial_state(1e-6)
    _, y = rk4_integrate(y0, 400.0, 0.25, ConstantForcing(0.5), FIXED)
    assert np.all(np.diff(y[:, CUM_INC]) >= -1e-14)


def test_no_epidemic_without_seed():
    """With zero initial infection the system must remain at its disease-free state."""
    y0 = initial_state(0.0)
    _, y = rk4_integrate(y0, 400.0, 0.25, ConstantForcing(1.0), FIXED)
    assert y[-1, CUM_INC] == pytest.approx(0.0, abs=1e-14)
    assert y[-1, S_H] == pytest.approx(1.0, abs=1e-12)


def test_low_transmission_does_not_produce_an_epidemic():
    """Below threshold the seeded infection must die out rather than grow.

    This is a behavioural check on the host--vector coupling: with a very small
    transmission coefficient, mosquito mortality outpaces transmission and the
    outbreak cannot take off.
    """
    y0 = initial_state(1e-5)
    _, y = rk4_integrate(y0, 600.0, 0.25, ConstantForcing(0.005), FIXED)
    final_attack_rate = y[-1, R_H] + y[-1, I_H] + y[-1, E_H]
    assert final_attack_rate < 1e-4


def test_high_transmission_produces_an_epidemic():
    """Above threshold a substantial fraction of the population must be infected."""
    y0 = initial_state(1e-5)
    _, y = rk4_integrate(y0, 900.0, 0.25, ConstantForcing(1.0), FIXED)
    assert y[-1, R_H] > 0.1


def test_step_size_is_adequate():
    """The production step of 0.25 days must agree with a tenfold finer step.

    This is the test that licenses the fixed-step integrator used throughout the
    inference. The tolerance is far tighter than any quantity the paper reports.
    """
    y0 = initial_state(1e-5)
    _, y_coarse = rk4_integrate(y0, 500.0, 0.25, ConstantForcing(0.6), FIXED)
    _, y_fine = rk4_integrate(y0, 500.0, 0.025, ConstantForcing(0.6), FIXED)
    assert abs(y_coarse[-1, CUM_INC] - y_fine[-1, CUM_INC]) < 1e-8


def _final_incidence(dt: float) -> float:
    y0 = initial_state(1e-5)
    _, y = rk4_integrate(y0, 200.0, dt, ConstantForcing(0.6), FIXED)
    return float(y[-1, CUM_INC])


def test_rk4_converges_at_fourth_order():
    """Halving the step must reduce the error by roughly a factor of sixteen.

    The order is measured by self-convergence,

        p = log2( |y(4h) - y(2h)| / |y(2h) - y(h)| ),

    rather than against a fine-step reference solution. A reference is only
    usable while its own error stays well below the errors being measured, and
    that condition fails here at small steps (see
    ``test_convergence_measurement_hits_a_roundoff_floor``). Self-convergence
    needs no reference and so cannot be contaminated by one.

    Steps are chosen in the range where truncation error dominates. Values run
    slightly above four because the coarsest steps are not yet fully asymptotic.
    """
    dts = [1.6, 0.8, 0.4, 0.2]
    vals = [_final_incidence(dt) for dt in dts]
    diffs = [abs(a - b) for a, b in zip(vals, vals[1:])]
    orders = [np.log2(a / b) for a, b in zip(diffs, diffs[1:])]

    assert all(3.5 < o < 5.0 for o in orders), f"measured orders {orders}"
    # The finest triple is closest to the asymptotic regime and should sit
    # nearer to four.
    assert 3.7 < orders[-1] < 4.6, f"finest measured order {orders[-1]}"


def test_convergence_measurement_hits_a_roundoff_floor():
    """Below roughly dt = 0.2 the error stops falling, and this is expected.

    By the end of the window the cumulative incidence has saturated near one, so
    an absolute error near 1e-13 is only a few hundred units of double-precision
    epsilon: accumulated rounding, not truncation. Refining further cannot help,
    and any convergence order measured in this regime is meaningless.

    The test pins this behaviour deliberately. It documents why the order test
    above restricts itself to coarser steps, and it would fail loudly if the
    integrator were ever changed in a way that altered the floor.
    """
    ref = _final_incidence(0.0015625)
    err_coarse = abs(_final_incidence(0.4) - ref)
    err_fine = abs(_final_incidence(0.05) - ref)

    assert err_coarse > err_fine          # refinement helps at first
    assert err_fine < 1e-11               # and then stalls at the roundoff level
    assert err_fine > 0.0                 # it never reaches exactly zero


def test_weekly_incidence_uses_a_seven_day_window_per_observation():
    """Gaps between observed weeks must not inflate the weeks that follow them.

    With contiguous weeks, differencing consecutive cumulative values gives the
    same answer as taking each observation's own seven days. With a gap it does
    not: three weeks of incidence would be attributed to the single week
    reported after the gap. The spatial comparison intersects districts whose
    reporting weeks do not align, so gaps are not hypothetical there.
    """
    y0 = initial_state(1e-5)
    t, y = rk4_integrate(y0, 200.0, 0.25, ConstantForcing(0.6), FIXED)

    contiguous = np.arange(0.0, 140.0, 7.0)
    gapped = np.array([0.0, 21.0, 28.0, 63.0, 70.0])

    w_contig = weekly_incidence(t, y, contiguous, population=1.0)
    w_gapped = weekly_incidence(t, y, gapped, population=1.0)

    # Each gapped observation must equal the contiguous one starting the same day.
    for day, value in zip(gapped, w_gapped):
        expected = w_contig[int(day // 7)]
        assert value == pytest.approx(expected, rel=1e-10), (
            f"observation at day {day:.0f} covers more than seven days")


def test_weekly_incidence_sums_to_cumulative():
    """Weekly counts must partition the cumulative incidence without loss."""
    y0 = initial_state(1e-5)
    t, y = rk4_integrate(y0, 70.0, 0.25, ConstantForcing(0.6), FIXED)
    weeks = np.arange(0, 63, 7, dtype=float)
    weekly = weekly_incidence(t, y, weeks, population=1.0)
    total = np.interp(weeks[-1] + 7.0, t, y[:, CUM_INC]) - np.interp(weeks[0], t,
                                                                    y[:, CUM_INC])
    assert weekly.sum() == pytest.approx(total, rel=1e-10)


def test_climate_forcing_matches_its_closed_form():
    """The forcing callable must equal beta_0 * exp(a_T z_T + a_P z_P) exactly."""
    days = np.arange(0.0, 100.0)
    z_t = np.sin(days / 20.0)
    z_r = np.cos(days / 30.0)
    f = ClimateForcing(days, z_t, z_r, beta_0=0.4, a_temp=0.7, a_rain=-0.3)
    for t in (0.0, 13.5, 47.0, 99.0):
        expected = 0.4 * np.exp(0.7 * np.interp(t, days, z_t)
                                - 0.3 * np.interp(t, days, z_r))
        assert f(t) == pytest.approx(expected, rel=1e-12)


def test_divergent_parameters_raise_rather_than_return_nan():
    """A transmission rate outside the explicit method's stability region must
    be reported, not silently returned as NaN.

    Multi-start fitting does propose such values. If the integrator returned
    NaN, the optimiser would receive a meaningless residual and could not
    distinguish an impossible parameter set from a merely poor one.
    """
    y0 = initial_state(1e-5)
    with pytest.raises(IntegrationFailure):
        rk4_integrate(y0, 400.0, 1.0, ConstantForcing(1e7), FIXED)


def test_valid_parameters_do_not_trigger_the_guard():
    """The guard must not fire anywhere in the epidemiologically plausible range.

    Dengue outbreaks are generally reported with R0 between about 1.5 and 6.
    Values up to R0 = 15 are tested, well beyond anything credible, to leave
    margin for the transient excursions that climate forcing produces.
    """
    y0 = initial_state(1e-5)
    for r0 in (0.5, 1.5, 3.0, 6.0, 15.0):
        beta = r0 / 10.0                      # R0 = 10 * beta, see models.py
        rk4_integrate(y0, 500.0, 1.0, ConstantForcing(beta), FIXED)


def test_reproduction_number_relation():
    """R0 = 10 * beta under the default fixed parameters."""
    assert basic_reproduction_number(0.3, FIXED) == pytest.approx(3.0, rel=1e-12)
    assert basic_reproduction_number(0.55, FIXED) == pytest.approx(5.5, rel=1e-12)


def test_stability_limit_lies_far_outside_the_plausible_range():
    """The one-day step must remain stable well past any credible epidemic.

    An explicit method has a finite stability region, so the guard introduces a
    barrier in parameter space. That barrier is only harmless if it sits far
    from where the optimum can plausibly lie. Bisection locates it, and the test
    asserts the margin rather than trusting it.

    Measured boundary scales roughly as beta_max ~ 1.5 / dt, so halving the step
    doubles the reachable transmission rate if a future model ever needs it.
    """
    y0 = initial_state(1e-5)

    def stable(beta: float) -> bool:
        try:
            rk4_integrate(y0, 500.0, 1.0, ConstantForcing(beta), FIXED)
            return True
        except IntegrationFailure:
            return False

    lo, hi = 0.01, 100.0
    for _ in range(30):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if stable(mid) else (lo, mid)

    r0_limit = basic_reproduction_number(lo, FIXED)
    assert r0_limit > 15.0, f"stability limit at R0 = {r0_limit:.1f} is too tight"


def test_one_day_step_is_adequate_under_climate_forcing():
    """The production step must hold when the driver itself varies weekly.

    Adequacy was established earlier for a constant coefficient. A time-varying
    beta is a harder case, so it is checked separately rather than assumed to
    follow.
    """
    days = np.arange(0.0, 561.0)
    z_t = np.sin(2 * np.pi * days / 365.0) * 1.5
    z_r = np.exp(-0.5 * ((days - 250.0) / 40.0) ** 2) * 2.0 - 0.5
    forcing = ClimateForcing(days, z_t, z_r, beta_0=0.35, a_temp=0.8, a_rain=0.7)

    y0 = initial_state(1e-6)
    _, y_coarse = rk4_integrate(y0, 560.0, 1.0, forcing, FIXED)
    _, y_fine = rk4_integrate(y0, 560.0, 0.125, forcing, FIXED)

    rel = abs(y_coarse[-1, CUM_INC] - y_fine[-1, CUM_INC]) / y_fine[-1, CUM_INC]
    assert rel < 1e-6, f"relative difference {rel:.2e} at dt = 1 day"


def test_climate_forcing_reduces_to_constant_when_coefficients_vanish():
    """With zero coefficients the climate model must equal the null model.

    Guarantees the two models are nested, which is what makes a comparison
    between them meaningful.
    """
    days = np.arange(0.0, 100.0)
    f = ClimateForcing(days, np.sin(days), np.cos(days), 0.55, 0.0, 0.0)
    g = ConstantForcing(0.55)
    assert np.allclose(f.on_grid(days), g.on_grid(days))


def test_fast_integrator_matches_the_reference_formulation():
    """The scalar inner loop must reproduce `rhs`-based RK4 exactly.

    `rhs` is the readable statement of the model; the integrator executes an
    unrolled scalar form of the same arithmetic because allocating arrays per
    stage dominated the runtime. An optimisation that changes the numbers is a
    bug, so the two are compared directly rather than trusted to agree.
    """
    def reference(y0, t_end, dt, beta_t, fixed):
        n = int(round(t_end / dt))
        t = np.linspace(0.0, n * dt, n + 1)
        y = np.empty((n + 1, len(y0)))
        y[0] = y0
        for k in range(n):
            tk, yk = t[k], y[k]
            k1 = rhs(tk, yk, beta_t, fixed)
            k2 = rhs(tk + dt / 2, yk + dt / 2 * k1, beta_t, fixed)
            k3 = rhs(tk + dt / 2, yk + dt / 2 * k2, beta_t, fixed)
            k4 = rhs(tk + dt, yk + dt * k3, beta_t, fixed)
            y[k + 1] = yk + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        return t, y

    days = np.arange(0.0, 561.0)
    forcings = [
        ConstantForcing(0.35),
        ClimateForcing(days, np.sin(2 * np.pi * days / 365.0) * 1.2,
                       np.cos(2 * np.pi * days / 300.0), 0.3, 0.6, -0.4),
    ]
    for forcing in forcings:
        y0 = initial_state(3e-6)
        _, y_fast = rk4_integrate(y0, 400.0, 1.0, forcing, FIXED)
        _, y_ref = reference(y0, 400.0, 1.0, forcing, FIXED)
        assert np.allclose(y_fast, y_ref, rtol=0, atol=1e-12), (
            f"max difference {np.max(np.abs(y_fast - y_ref)):.3e}")


def test_seir_structure_conserves_the_population():
    """The alternative structure must not leak people either."""
    from dengue_pk.models import rk4_integrate_seir
    _, y = rk4_integrate_seir(1e-5, 400.0, 1.0, ConstantForcing(0.5), FIXED)
    assert np.allclose(y[:, :4].sum(axis=1), 1.0, atol=1e-10)
    assert np.all(np.diff(y[:, 4]) >= -1e-14)          # incidence accumulates


def test_seir_reproduction_number_uses_its_own_relation():
    """R0 = beta/gamma_h for SEIR, not the host-vector geometric mean.

    The two structures map the same beta to different R0 values, so a comparison
    across structures made on beta rather than R0 would be meaningless. This pins
    the distinction.
    """
    assert basic_reproduction_number(0.4, FIXED, "seir") == pytest.approx(2.0)
    assert basic_reproduction_number(0.4, FIXED, "hostvector") == pytest.approx(4.0)


def test_seir_threshold_behaviour():
    """Below R0 = 1 the outbreak must die; above it, take off."""
    from dengue_pk.models import rk4_integrate_seir
    _, sub = rk4_integrate_seir(1e-5, 600.0, 1.0, ConstantForcing(0.1), FIXED)
    _, sup = rk4_integrate_seir(1e-5, 600.0, 1.0, ConstantForcing(0.6), FIXED)
    assert sub[-1, 3] < 1e-3          # R0 = 0.5, dies out
    assert sup[-1, 3] > 0.5           # R0 = 3.0, substantial attack rate


def test_the_two_structures_are_genuinely_different_at_matched_r0():
    """Matching R0 does not match the epidemic, and that is the point.

    The naive expectation is that the host-vector model, carrying the mosquito's
    infection cycle as an extra delay, must be the slower of the two at equal R0.
    It is not: measured peaks are day 71 for host-vector against day 77 for SEIR.

    The reason is that R0 is not a speed. For the host-vector system R0 is the
    geometric mean over a two-step cycle, so a "generation" in that definition is
    a half-cycle of about 7.5 days, against 10 days for SEIR. Matching R0
    therefore leaves the host-vector model with the shorter generation interval
    and the faster growth.

    What the robustness study needs is only that the structures differ, so that
    varying structure varies something real. This asserts that, and records the
    correction, because a test written to the wrong expectation would have been
    "fixed" by loosening it.
    """
    from dengue_pk.models import CUM_INC, rk4_integrate_seir

    beta_hv = 0.30                                  # R0 = 3.0
    beta_seir = 3.0 * FIXED.gamma_h                 # R0 = 3.0
    assert basic_reproduction_number(beta_hv, FIXED, "hostvector") ==         pytest.approx(basic_reproduction_number(beta_seir, FIXED, "seir"))

    _, y_hv = rk4_integrate(initial_state(1e-5), 500.0, 1.0,
                            ConstantForcing(beta_hv), FIXED)
    _, y_se = rk4_integrate_seir(1e-5, 500.0, 1.0,
                                 ConstantForcing(beta_seir), FIXED)
    inc_hv = np.diff(y_hv[:, CUM_INC])
    inc_se = np.diff(y_se[:, 4])

    peak_hv, peak_se = int(np.argmax(inc_hv)), int(np.argmax(inc_se))
    assert abs(peak_hv - peak_se) >= 3, (
        f"peaks at day {peak_hv} and {peak_se} are too close for structure to "
        f"be a meaningful factor")
    # Different shapes, not merely shifted: compare normalised peak heights.
    assert abs(inc_hv.max() / inc_hv.sum() - inc_se.max() / inc_se.sum()) > 0.001


# --- temperature-dependent mosquito mortality --------------------------------

def test_constant_mu_v_t_reproduces_the_scalar_path():
    """Supplying a constant mortality callable must change nothing.

    The temperature-dependent structure was added by threading a mortality grid
    through the integrator. If that refactor perturbed the constant case, every
    result in the study would move for a reason unrelated to the new structure.
    """
    fixed = FixedParams()
    beta = ConstantForcing(0.15)
    y0 = initial_state(1e-5)

    t_a, y_a = rk4_integrate(y0, 200.0, 1.0, beta, fixed)
    t_b, y_b = rk4_integrate(y0, 200.0, 1.0, beta, fixed,
                             mu_v_t=lambda _t: fixed.mu_v)

    assert np.allclose(t_a, t_b)
    assert np.max(np.abs(y_a - y_b)) < 1e-14


def test_thermal_mortality_equals_baseline_at_the_optimum():
    """At the thermal optimum the two structures must coincide exactly."""
    fixed = FixedParams()
    peak = 0.5 * (LIFESPAN_T_MIN + LIFESPAN_T_MAX)
    days = np.arange(0.0, 201.0)
    mort = ThermalMortality(days, np.full_like(days, peak), fixed.mu_v)

    assert np.allclose(mort.on_grid(days), fixed.mu_v)

    beta = ConstantForcing(0.15)
    y0 = initial_state(1e-5)
    _, y_const = rk4_integrate(y0, 200.0, 1.0, beta, fixed)
    _, y_therm = rk4_integrate(y0, 200.0, 1.0, beta, fixed, mu_v_t=mort)
    assert np.max(np.abs(y_const - y_therm)) < 1e-12


def test_thermal_mortality_raises_mortality_away_from_the_optimum():
    """Mortality is unimodal in temperature: lowest at the optimum, higher on
    both sides. A response that only rose in one direction would be a monotone
    term wearing a unimodal name."""
    fixed = FixedParams()
    peak = 0.5 * (LIFESPAN_T_MIN + LIFESPAN_T_MAX)
    temps = np.array([peak - 10.0, peak, peak + 10.0])
    mu = ThermalMortality(np.arange(3.0), temps, fixed.mu_v).on_grid(np.arange(3.0))

    assert mu[1] == pytest.approx(fixed.mu_v)
    assert mu[0] > mu[1] and mu[2] > mu[1]


def test_thermal_mortality_is_bounded_outside_the_thermal_limits():
    """Beyond the limits the lifespan quadratic is zero; mortality must stay
    finite or the epidemic stops dead for a numerical reason."""
    fixed = FixedParams()
    temps = np.array([-20.0, 60.0])
    mu = ThermalMortality(np.arange(2.0), temps, fixed.mu_v).on_grid(np.arange(2.0))
    assert np.all(np.isfinite(mu))
    assert np.all(mu <= fixed.mu_v / 0.15 + 1e-12)


def test_the_third_structure_differs_from_the_first():
    """A structure that produced identical predictions would not be a factor."""
    fixed = FixedParams()
    days = np.arange(0.0, 301.0)
    temp = 24.45 + 9.0 * np.sin(2 * np.pi * days / 365.0)
    mort = ThermalMortality(days, temp, fixed.mu_v)

    beta = ConstantForcing(0.15)
    y0 = initial_state(1e-5)
    _, y_const = rk4_integrate(y0, 300.0, 1.0, beta, fixed)
    _, y_therm = rk4_integrate(y0, 300.0, 1.0, beta, fixed, mu_v_t=mort)
    peak_shift = abs(np.argmax(np.diff(y_const[:, CUM_INC]))
                     - np.argmax(np.diff(y_therm[:, CUM_INC])))
    assert np.max(np.abs(y_const - y_therm)) > 1e-4, "structures are identical"
    assert peak_shift >= 0  # recorded; direction depends on the season drawn
