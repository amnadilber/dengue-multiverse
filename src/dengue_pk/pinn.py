"""
Physics-informed neural network for the host--vector inverse problem.

The classical approach solves the model forward thousands of times inside an
optimiser. A PINN never solves it at all: a network represents the solution,
automatic differentiation supplies its derivative exactly, and training drives
the residual of the governing equations toward zero while simultaneously fitting
the observed counts. The transmission parameters are trainable variables
alongside the network weights.

The argument made for this trade is that the cost does not scale with the cost of
a forward solve. On a seven-state ODE integrated in twelve milliseconds that buys
nothing; the question this module exists to answer is what it costs in accuracy
to find out.

Three features of this problem make it considerably harder than the textbook SIR
demonstrations:

* **Scale separation.** The susceptible fraction is order one while the infected
  fraction peaks near 10^-2 and the initial seed is near 10^-6. A network with a
  single output scale cannot resolve all three, so the infected and exposed
  compartments are represented in log space.
* **The observable is not a state.** Surveillance reports weekly counts, which are
  differences of cumulative incidence. The data term must therefore difference
  the network's output at week boundaries rather than read a state directly.
* **Time-varying coefficients.** The transmission rate is driven by a climate
  series defined only at weekly resolution, so it enters as an interpolant.

The same negative binomial likelihood used by the classical fit serves as the data
term, so the two approaches are optimising comparable objectives and a difference
in result cannot be attributed to a difference in loss.
"""

from __future__ import annotations

import time

import numpy as np

# The MSVC runtime shim in dengue_pk/__init__.py must already have run, which it
# has by the time this module can be imported at all. Doing it here instead would
# be too late: NumPy would already have pinned the older runtime.
import tensorflow as tf

from .climate import briere  # noqa: E402

# State layout inside the network. The exposed and infected compartments are
# carried as natural logarithms because they span several orders of magnitude
# within a single window.
LOG_STATES = ("e_h", "i_h", "i_v")
STATE_ORDER = ("s_h", "e_h", "i_h", "r_h", "s_v", "i_v", "cum")


class InversePINN:
    """Joint estimation of the solution and the transmission parameters."""

    def __init__(self, data, cfg, fixed, model: str = "climate", seed: int = 0):
        tf.random.set_seed(seed)
        np.random.seed(seed)

        self.data = data
        self.cfg = cfg
        self.fixed = fixed
        self.model = model

        pc = cfg["inference"]["pinn"]
        self.T = float(data.days[-1] + 7.0)            # horizon, days
        self.n_col = int(pc["collocation_points"])
        self.epochs = int(pc["epochs"])
        self.data_weight = float(pc["data_weight"])

        tr = cfg["model"]["temperature_response"]
        self.b_temp = briere(data.temp_c, tr["t_min_c"], tr["t_max_c"])

        # Climate covariates on the collocation grid, interpolated once.
        self.tau_col = np.linspace(0.0, 1.0, self.n_col).reshape(-1, 1)
        days_col = self.tau_col.ravel() * self.T
        self.b_col = tf.constant(
            np.interp(days_col, data.days, self.b_temp).reshape(-1, 1),
            dtype=tf.float32)
        self.zr_col = tf.constant(
            np.interp(days_col, data.days, data.z_rain).reshape(-1, 1),
            dtype=tf.float32)
        self.tau_col_tf = tf.constant(self.tau_col, dtype=tf.float32)

        # Week boundaries, where cumulative incidence must be differenced.
        edges = np.concatenate([data.days, [data.days[-1] + 7.0]]) / self.T
        self.tau_edges = tf.constant(edges.reshape(-1, 1), dtype=tf.float32)
        self.cases = tf.constant(data.cases.reshape(-1, 1), dtype=tf.float32)

        self.net = self._build(pc)

        est = cfg["model"]["estimated"]
        self.log_beta0 = tf.Variable(np.log(est["beta_0"]["init"]), dtype=tf.float32)
        self.log_atemp = tf.Variable(np.log(est["a_temp"]["init"]), dtype=tf.float32)
        self.a_rain = tf.Variable(est["a_rain"]["init"], dtype=tf.float32)
        self.logit_popfrac = tf.Variable(
            np.log(est["pop_frac"]["init"] / (1 - est["pop_frac"]["init"])),
            dtype=tf.float32)
        self.log_i0 = tf.Variable(np.log(est["i0_frac"]["init"]), dtype=tf.float32)

        self.rho = float(cfg["model"]["fixed"]["rho_fixed"])
        self.nb_k = tf.Variable(5.0, dtype=tf.float32, trainable=False)

        self.opt = tf.keras.optimizers.Adam(
            tf.keras.optimizers.schedules.ExponentialDecay(
                float(pc["learning_rate"]), int(pc["lr_decay_steps"]),
                float(pc["lr_decay_rate"])))

    @staticmethod
    def _build(pc) -> tf.keras.Model:
        layers = [tf.keras.layers.Input(shape=(1,))]
        for width in pc["layers"]:
            layers.append(tf.keras.layers.Dense(width, activation=pc["activation"]))
        layers.append(tf.keras.layers.Dense(len(STATE_ORDER)))
        return tf.keras.Sequential(layers)

    # -- parameters ---------------------------------------------------------
    def params(self) -> dict:
        return dict(beta_0=float(tf.exp(self.log_beta0)),
                    a_temp=float(tf.exp(self.log_atemp)),
                    a_rain=float(self.a_rain),
                    pop_frac=float(tf.sigmoid(self.logit_popfrac)),
                    i0_frac=float(tf.exp(self.log_i0)))

    def _beta(self, b, zr):
        return (tf.exp(self.log_beta0)
                * tf.pow(tf.maximum(b, 1e-6), tf.exp(self.log_atemp))
                * tf.exp(self.a_rain * zr)) if self.model == "climate" \
            else tf.exp(self.log_beta0) * tf.ones_like(b)

    # -- trial solution -----------------------------------------------------
    def _solution(self, tau):
        """Hard-constrained trial solution.

        Compartments that span orders of magnitude are represented in log space,
        so the network's raw output is an increment to the log of the initial
        value rather than to the value itself. The multiplication by tau makes
        every initial condition exact for any weights, which removes the trivial
        equilibrium that a soft initial-condition penalty admits.
        """
        raw = self.net(tau)
        i0 = tf.exp(self.log_i0)

        s_h = 1.0 - i0 + tau * raw[:, 0:1]
        e_h = tf.exp(tf.math.log(i0 * 0.5) + tau * raw[:, 1:2])
        i_h = tf.exp(tf.math.log(i0) + tau * raw[:, 2:3])
        r_h = tau * tf.nn.softplus(raw[:, 3:4])
        s_v = 1.0 + tau * raw[:, 4:5]
        i_v = tf.exp(tf.math.log(i0 * 1e-2) + tau * raw[:, 5:6])
        cum = tau * tf.nn.softplus(raw[:, 6:7])
        return tf.concat([s_h, e_h, i_h, r_h, s_v, i_v, cum], axis=1)

    # -- losses -------------------------------------------------------------
    @tf.function
    def _loss(self):
        with tf.GradientTape() as tape:
            tape.watch(self.tau_col_tf)
            y = self._solution(self.tau_col_tf)
        dy = tape.batch_jacobian(y, self.tau_col_tf)[:, :, 0]

        s_h, e_h, i_h = y[:, 0], y[:, 1], y[:, 2]
        s_v, i_v = y[:, 4], y[:, 5]
        beta = self._beta(self.b_col[:, 0], self.zr_col[:, 0])

        f = self.fixed
        lam_h = beta * f.vector_host_ratio * i_v
        lam_v = beta * i_h

        # Chain rule from days to tau contributes the factor T.
        rhs = tf.stack([
            -lam_h * s_h,
            lam_h * s_h - f.sigma_h * e_h,
            f.sigma_h * e_h - f.gamma_h * i_h,
            f.gamma_h * i_h,
            f.mu_v * (1.0 - s_v) - lam_v * s_v,
            lam_v * s_v - f.mu_v * i_v,
            f.sigma_h * e_h,
        ], axis=1) * self.T

        # Residuals are scaled by the magnitude of each compartment. Without
        # this the susceptible equation, whose terms are order one, would
        # dominate the infected equations, whose terms are order 1e-4, and the
        # epidemic itself would be fitted last.
        scale = tf.maximum(tf.abs(y), 1e-6)
        physics = tf.reduce_mean(tf.square((dy - rhs) / scale))

        # Data term: weekly cases are differences of cumulative incidence.
        cum_edges = self._solution(self.tau_edges)[:, 6:7]
        weekly = (cum_edges[1:] - cum_edges[:-1]) \
            * tf.sigmoid(self.logit_popfrac) * self.data.population * self.rho
        mu = tf.maximum(weekly, 1e-6)

        k = self.nb_k
        nb = (tf.math.lgamma(self.cases + k) - tf.math.lgamma(k)
              - tf.math.lgamma(self.cases + 1.0)
              + k * tf.math.log(k / (k + mu))
              + self.cases * tf.math.log(mu / (k + mu)))
        data_loss = -tf.reduce_mean(nb)

        return physics + self.data_weight * data_loss, physics, data_loss

    @tf.function
    def _step(self):
        variables = (self.net.trainable_variables
                     + [self.log_beta0, self.a_rain, self.logit_popfrac,
                        self.log_i0]
                     + ([self.log_atemp] if self.model == "climate" else []))
        with tf.GradientTape() as tape:
            total, physics, data_loss = self._loss()
        self.opt.apply_gradients(zip(tape.gradient(total, variables), variables))
        return total, physics, data_loss

    # -- driver -------------------------------------------------------------
    def train(self, verbose_every: int = 5000, k_update_every: int = 2000):
        history = []
        t0 = time.time()
        for epoch in range(self.epochs + 1):
            total, physics, data_loss = self._step()

            # Refresh the dispersion periodically from the current fit, mirroring
            # the alternating scheme the classical estimator uses.
            if epoch and epoch % k_update_every == 0:
                mu = self.weekly_mu()
                y = self.data.cases
                from .inference import estimate_dispersion_k
                self.nb_k.assign(np.float32(
                    np.clip(estimate_dispersion_k(y, mu), 0.1, 1e4)))

            if epoch % 200 == 0:
                history.append((epoch, float(total), float(physics),
                                float(data_loss), *self.params().values()))
            if verbose_every and epoch % verbose_every == 0:
                p = self.params()
                print(f"    epoch {epoch:6d}  total {float(total):10.4f}  "
                      f"physics {float(physics):9.3e}  "
                      f"beta0 {p['beta_0']:.4f}  popfrac {p['pop_frac']:.5f}")
        self.seconds = time.time() - t0
        self.history = np.array(history)
        return self

    def weekly_mu(self) -> np.ndarray:
        cum = self._solution(self.tau_edges).numpy()[:, 6]
        weekly = np.diff(cum) * self.params()["pop_frac"] \
            * self.data.population * self.rho
        return np.maximum(weekly, 1e-9)
