"""
Inference validation.

The decisive test is recovery on synthetic data: if the fitting machinery cannot
recover parameters it generated itself, from data with no model misspecification,
then nothing it reports on real data can be believed. Everything else here checks
a property the estimator must satisfy regardless of the data.
"""

import numpy as np
import pytest

from dengue_pk import load_config
from dengue_pk.climate import ClimateForcing
from dengue_pk.inference import (Dataset, estimate_dispersion_k,
                                 nb_deviance_residuals, nb_loglik,
                                 poisson_deviance_residuals, poisson_loglik,
                                 pack, unpack, predict, fit)
from dengue_pk.models import FixedParams

import copy

# Tests use fewer multi-starts than the analysis. Ten restarts on every fit in
# every test pushed the suite past nine minutes, which is long enough that it
# stops being run — and a test suite that is not run protects nothing. Four
# starts still exercises the multi-start machinery and its agreement check.
# Anything reported in the paper comes from `scripts/`, which uses the full
# configured count.
CFG = copy.deepcopy(load_config())
CFG["inference"]["classical"]["n_restarts"] = 4

FIXED = FixedParams.from_config(CFG)


def synthetic_window(n_weeks=60, beta_0=0.32, a_temp=1.0, a_rain=0.4,
                     pop_frac=0.02, i0=2e-7, seed=0):
    """Generate a window from known parameters, with Poisson observation noise."""
    rng = np.random.default_rng(seed)
    days = np.arange(n_weeks, dtype=float) * 7.0
    # A seasonal climate: temperature sweeping through the mosquito's viable
    # range, rainfall peaking mid-window.
    temp_c = 27.0 + 8.0 * np.sin(2 * np.pi * days / 365.0)
    z_rain = np.exp(-0.5 * ((days - days.mean()) / 60.0) ** 2) * 2.0 - 0.5

    truth = dict(beta_0=beta_0, a_temp=a_temp, a_rain=a_rain,
                 pop_frac=pop_frac, i0_frac=i0)
    data = Dataset(days, np.zeros(n_weeks), temp_c, z_rain,
                   population=50e6, label="synthetic")
    mu = predict(truth, data, FIXED, "climate", CFG)
    data.cases = rng.poisson(mu).astype(float)
    return data, truth, mu


def test_parameter_transforms_round_trip():
    theta = dict(beta_0=0.37, a_temp=1.2, a_rain=0.8, pop_frac=0.023,
                 i0_frac=4e-7)
    names = list(theta)
    assert unpack(pack(theta, names), names) == pytest.approx(theta, rel=1e-10)


def test_population_fraction_stays_in_range():
    """Any real number must map to a population fraction strictly inside (0, 1)."""
    for x in (-50.0, -3.0, 0.0, 3.0, 50.0):
        p = unpack(np.array([x]), ["pop_frac"])["pop_frac"]
        assert 0.0 < p < 1.0


def test_reporting_fraction_and_population_are_confounded():
    """Only the product of the two is identifiable, so only one may be estimated.

    The state variables are fractions, so the population never enters the
    dynamics — only the conversion of incidence into counts. Halving the
    population at risk while doubling the reporting fraction must therefore leave
    the prediction bit-for-bit unchanged. This is why `rho` is fixed from the
    literature and `pop_frac` estimated, rather than both being free.
    """
    days = np.arange(40.0) * 7.0
    data = Dataset(days, np.zeros(40), np.full(40, 28.0), np.zeros(40),
                   population=50e6)
    theta = dict(beta_0=0.35, a_temp=1.0, a_rain=0.0, pop_frac=0.02,
                 i0_frac=1e-6)

    cfg_a = {**CFG, "model": {**CFG["model"],
                              "fixed": {**CFG["model"]["fixed"], "rho_fixed": 0.05}}}
    cfg_b = {**CFG, "model": {**CFG["model"],
                              "fixed": {**CFG["model"]["fixed"], "rho_fixed": 0.10}}}

    mu_a = predict(theta, data, FIXED, "climate", cfg_a)
    mu_b = predict({**theta, "pop_frac": 0.01}, data, FIXED, "climate", cfg_b)
    assert np.allclose(mu_a, mu_b, rtol=1e-12)


def test_deviance_residuals_vanish_at_a_perfect_fit():
    """Including at zero observations, which the safety floor once biased.

    A floor applied to the whole expression rather than only inside the
    logarithm gave a zero week predicted as zero a residual of about -4e-5
    instead of zero — a systematic bias on every quiet week, and these windows
    contain many.
    """
    y = np.array([0.0, 5.0, 120.0, 900.0])
    assert np.allclose(poisson_deviance_residuals(y, y), 0.0, atol=1e-12)


def test_zero_weeks_are_penalised_in_proportion_to_the_prediction():
    """A zero observation must cost more the larger the prediction was."""
    y = np.zeros(3)
    r = poisson_deviance_residuals(y, np.array([0.01, 1.0, 50.0]))
    assert np.all(r < 0)                       # over-prediction is signed negative
    assert abs(r[0]) < abs(r[1]) < abs(r[2])


def test_deviance_residuals_are_signed():
    y = np.array([10.0, 10.0])
    r = poisson_deviance_residuals(y, np.array([5.0, 20.0]))
    assert r[0] > 0 and r[1] < 0          # under-prediction positive, over negative


def test_deviance_residual_sum_of_squares_equals_the_deviance():
    """The identity that makes least_squares equivalent to Poisson likelihood."""
    y = np.array([0.0, 3.0, 40.0, 260.0])
    mu = np.array([1.0, 4.5, 33.0, 300.0])
    saturated = poisson_loglik(y, np.maximum(y, 1e-12))
    deviance = 2.0 * (saturated - poisson_loglik(y, mu))
    assert np.sum(poisson_deviance_residuals(y, mu) ** 2) == pytest.approx(
        deviance, rel=1e-8)


def test_loglik_is_maximised_at_the_truth():
    """Perturbing the mean away from the data must reduce the likelihood."""
    y = np.array([2.0, 18.0, 95.0, 40.0])
    assert poisson_loglik(y, y) > poisson_loglik(y, y * 1.3)
    assert poisson_loglik(y, y) > poisson_loglik(y, y * 0.7)


def test_nb_reduces_to_poisson_as_dispersion_grows():
    """As k grows the negative binomial must converge on the Poisson.

    This is what makes the two models nested and the comparison between them
    meaningful — and it is the check that catches a mis-transcribed gamma term,
    which would otherwise produce plausible but wrong likelihoods.
    """
    y = np.array([0.0, 3.0, 40.0, 260.0])
    mu = np.array([1.0, 4.5, 33.0, 300.0])
    assert nb_loglik(y, mu, 1e9) == pytest.approx(poisson_loglik(y, mu), rel=1e-6)
    assert np.allclose(nb_deviance_residuals(y, mu, 1e9),
                       poisson_deviance_residuals(y, mu), atol=1e-4)


def test_nb_is_more_tolerant_of_overdispersion_than_poisson():
    """A wildly over-dispersed observation must cost less under the NB."""
    y = np.array([500.0])
    mu = np.array([100.0])
    assert abs(nb_deviance_residuals(y, mu, 2.0)[0]) < \
        abs(poisson_deviance_residuals(y, mu)[0])


def test_nb_deviance_residuals_vanish_at_a_perfect_fit():
    y = np.array([0.0, 5.0, 120.0, 900.0])
    assert np.allclose(nb_deviance_residuals(y, y, 3.0), 0.0, atol=1e-10)


def test_dispersion_estimator_recovers_a_known_k():
    """Simulate from a known NB and check the profiled estimate finds it."""
    rng = np.random.default_rng(0)
    mu = np.full(4000, 80.0)
    k_true = 4.0
    # NB2 as a gamma-Poisson mixture
    lam = rng.gamma(shape=k_true, scale=mu / k_true)
    y = rng.poisson(lam).astype(float)
    assert estimate_dispersion_k(y, mu) == pytest.approx(k_true, rel=0.15)


def test_dispersion_estimator_returns_a_large_k_for_poisson_data():
    """Poisson data carry no excess variance, so k must be driven upward."""
    rng = np.random.default_rng(1)
    mu = np.full(4000, 80.0)
    y = rng.poisson(mu).astype(float)
    assert estimate_dispersion_k(y, mu) > 100.0


def test_predictions_are_non_negative_and_finite():
    data, truth, _ = synthetic_window()
    mu = predict(truth, data, FIXED, "climate", CFG)
    assert np.all(np.isfinite(mu)) and np.all(mu > 0)


def test_climate_model_recovers_synthetic_parameters():
    """The test the whole analysis rests on.

    Data are generated from known parameters and refitted from scratch.

    The thermal exponent is deliberately **not** asserted here. A profile
    likelihood over it (docs/ANALYSIS_LOG.md, 2026-07-27) is flat to within 0.36
    log-likelihood units across a sixteen-fold range, far inside any confidence
    interval, while beta_0 and pop_frac stay pinned throughout. Asserting a
    tolerance on a parameter the data cannot determine would test the optimiser's
    arbitrary tie-breaking, not the estimator, and would fail or pass with the
    random seed.
    """
    data, truth, _ = synthetic_window(seed=1)
    res = fit(data, CFG, FIXED, model="climate")

    assert res.theta["beta_0"] == pytest.approx(truth["beta_0"], rel=0.25)
    assert res.theta["pop_frac"] == pytest.approx(truth["pop_frac"], rel=0.50)


def test_temperature_exponent_cannot_go_negative():
    """The Brière exponent is non-negative by construction.

    A negative exponent would mean transmission falling as conditions approach
    the thermal optimum, which was exactly the artefact the log-linear
    parameterisation produced on real data. The log transform makes it
    impossible rather than merely unlikely.
    """
    for x in (-50.0, -5.0, 0.0, 5.0):
        assert unpack(np.array([x]), ["a_temp"])["a_temp"] > 0.0


def test_multistart_reaches_a_consistent_optimum():
    """Independent starts should agree; disagreement would signal local minima."""
    data, _, _ = synthetic_window(seed=2)
    res = fit(data, CFG, FIXED, model="climate")
    assert res.n_converged >= 2, (
        f"only {res.n_converged}/{res.n_starts} starts reached the best optimum")


def test_transmission_is_recovered_even_where_the_thermal_exponent_is_not():
    """beta_0 must be recovered regardless of how well a_temp is pinned down.

    The Brière exponent is only weakly identified from a single season (see
    docs/ANALYSIS_LOG.md): beta_0 is fixed by transmission near the thermal
    optimum, while the exponent governs the fall-off away from it, which one
    epidemic barely visits. The conclusions rest on beta_0 and R0, so it is
    those that must be recoverable, and this test states that requirement
    separately rather than letting a loose tolerance on a_temp hide it.
    """
    for seed in (1, 5, 7):
        data, truth, _ = synthetic_window(seed=seed)
        res = fit(data, CFG, FIXED, model="climate")
        assert res.theta["beta_0"] == pytest.approx(truth["beta_0"], rel=0.30), (
            f"seed {seed}: beta_0 {res.theta['beta_0']:.4f} vs truth "
            f"{truth['beta_0']:.4f}")


def test_climate_model_beats_the_null_on_climate_driven_data():
    """AIC must prefer the climate model when the data really are climate-driven.

    This checks the model comparison itself: if the criterion cannot detect
    forcing that is known to be present, any preference it expresses on real
    data would be meaningless.
    """
    data, _, _ = synthetic_window(a_temp=1.5, a_rain=0.9, seed=3)
    climate = fit(data, CFG, FIXED, model="climate")
    null = fit(data, CFG, FIXED, model="constant")
    assert climate.aic < null.aic


def test_nb_fit_recovers_transmission_from_overdispersed_data():
    """Under-dispersed assumptions must not be required for a usable estimate.

    Data are generated with genuine extra-Poisson noise, which is what the real
    surveillance counts show. The negative binomial fit must still recover the
    transmission parameter; if it cannot, the observation model is not fit for
    the data it was introduced to handle.
    """
    rng = np.random.default_rng(11)
    data, truth, mu = synthetic_window(seed=11)
    k_true = 5.0
    lam = rng.gamma(shape=k_true, scale=np.maximum(mu, 1e-9) / k_true)
    data.cases = rng.poisson(lam).astype(float)

    res = fit(data, CFG, FIXED, model="climate", observation="nb")
    assert res.theta["beta_0"] == pytest.approx(truth["beta_0"], rel=0.35)
    assert res.nb_k is not None and 1.0 < res.nb_k < 50.0, (
        f"estimated k = {res.nb_k}")


def test_nb_fit_reports_dispersion_near_one():
    """The point of the NB model: residual dispersion should no longer be huge.

    Under Poisson these data give a Pearson dispersion far above one. Once the
    variance function includes the mu^2/k term, a value near one is the evidence
    that the extra variability has been modelled rather than ignored.
    """
    rng = np.random.default_rng(12)
    data, _, mu = synthetic_window(seed=12)
    lam = rng.gamma(shape=4.0, scale=np.maximum(mu, 1e-9) / 4.0)
    data.cases = rng.poisson(lam).astype(float)

    pois = fit(data, CFG, FIXED, model="climate", observation="poisson")
    nb = fit(data, CFG, FIXED, model="climate", observation="nb")
    assert nb.dispersion < pois.dispersion
    assert nb.dispersion < 5.0, f"NB dispersion still {nb.dispersion:.1f}"


def test_null_model_is_not_penalised_on_flat_data():
    """With no climate forcing present, the extra parameters must not pay off.

    The converse of the previous test, and the one that guards against a
    criterion biased toward the richer model.
    """
    data, _, _ = synthetic_window(a_temp=1e-6, a_rain=0.0, seed=4)
    climate = fit(data, CFG, FIXED, model="climate")
    null = fit(data, CFG, FIXED, model="constant")
    assert null.aic < climate.aic + 4.0


def test_temperature_form_switches_the_forcing_and_the_transform():
    """Both halves of the switch must move together.

    `loglinear` needs a free coefficient that can go negative; `briere` needs a
    non-negative exponent, which is what makes the wrong-signed artefact
    impossible. Switching the forcing without switching the transform would
    silently constrain the log-linear coefficient to be positive and destroy the
    comparison the global study rests on.
    """
    from dengue_pk.climate import ClimateForcing, MechanisticForcing
    from dengue_pk.inference import make_forcing, set_temperature_form

    days = np.arange(30.0) * 7.0
    data = Dataset(days, np.zeros(30), np.full(30, 27.0), np.zeros(30), 1e6)
    theta = dict(beta_0=0.3, a_temp=0.5, a_rain=0.1, pop_frac=0.01, i0_frac=1e-6)

    cfg_b = copy.deepcopy(CFG)
    cfg_b["model"]["temperature_form"] = "briere"
    set_temperature_form("briere")
    assert isinstance(make_forcing(theta, data, "climate", cfg_b),
                      MechanisticForcing)
    assert unpack(np.array([-5.0]), ["a_temp"])["a_temp"] > 0

    cfg_l = copy.deepcopy(CFG)
    cfg_l["model"]["temperature_form"] = "loglinear"
    set_temperature_form("loglinear")
    assert isinstance(make_forcing(theta, data, "climate", cfg_l),
                      ClimateForcing)
    assert unpack(np.array([-5.0]), ["a_temp"])["a_temp"] == pytest.approx(-5.0)

    set_temperature_form("briere")          # leave the module as found


def test_loglinear_form_can_return_a_negative_coefficient():
    """The conventional parameterisation must be able to produce the artefact.

    The global study's finding is that a monotone term reports a sign determined
    by the covariate's trajectory rather than by biology. That finding is only
    meaningful if the estimator was actually free to return either sign, so this
    fits data generated with transmission falling as temperature rises and
    checks that the recovered coefficient is negative.
    """
    from dengue_pk.climate import ClimateForcing
    from dengue_pk.inference import set_temperature_form
    from dengue_pk.models import simulate

    cfg = copy.deepcopy(CFG)
    cfg["model"]["temperature_form"] = "loglinear"
    set_temperature_form("loglinear")
    try:
        rng = np.random.default_rng(3)
        days = np.arange(52.0) * 7.0
        # Temperature falls steadily while the epidemic is seeded and grows.
        temp_c = np.linspace(32.0, 20.0, 52)
        z_temp, *_ = __import__("dengue_pk.climate", fromlist=["standardise"]) \
            .standardise(temp_c)
        forcing = ClimateForcing(days, z_temp, np.zeros(52), 0.30, -0.45, 0.0)
        mu, _ = simulate({"i0_frac": 2e-6, "rho": 0.05}, forcing, FIXED,
                         days, 2e6)
        data = Dataset(days, rng.poisson(np.maximum(mu, 1e-9)).astype(float),
                       temp_c, np.zeros(52), 2e6)

        res = fit(data, cfg, FIXED, model="climate", observation="nb")
        assert res.theta["a_temp"] < 0, (
            f"a_temp came back {res.theta['a_temp']:.3f}; the log-linear form "
            f"must be able to report a negative coefficient")
    finally:
        set_temperature_form("briere")
