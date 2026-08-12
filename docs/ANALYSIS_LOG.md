# Analysis log

A dated record of what was run, what it showed, and what was changed as a result.
Kept because the sequence of decisions is part of the result: a reader is entitled to
know which choices were made in advance and which were forced by the data.

---

## 2026-07-26 — First fit to real data

`scripts/03_fit_classical.py`, both models, all three windows, 10 multi-starts each.

### Headline

The climate-forced model wins decisively on AIC in every window (ΔAIC = −1000, −21, −35)
and **loses out-of-sample in every window** (held-out deviance +1306, +90, +366). That
combination is the signature of overfitting, not of a mechanism.

### Four results that indict the current model specification

**1. The climate coefficients have the wrong sign.**
Nationally `a_temp = −0.31 ± 0.09`, and in Khyber Pakhtunkhwa `a_temp = −0.38 ± 0.22`.
The model is saying transmission *falls* as temperature rises, which is backwards for
*Aedes*-borne disease over the temperature range observed.

The explanation is straightforward and damaging: the 2013 national epidemic peaks in
September, while temperature peaks in June. Through the entire growth phase, temperature
is falling. A negative coefficient therefore lets β rise on schedule — the covariate is
serving as a proxy for elapsed time since midsummer, not as a driver. Fitting a
descending covariate to an ascending epidemic will always find a negative coefficient,
whatever the biology.

**2. The reporting fraction is implausible.**
Estimated ρ ≈ 1 × 10⁻⁴ in every window: roughly one reported case per ten thousand
infections. Published dengue under-ascertainment is typically one in ten to one in thirty.

The cause is the population denominator. With R₀ = 2.28 and 188 million susceptible
people, the model's epidemic infects on the order of 10⁸ individuals, and the only way to
reconcile that with 17,894 reported cases is a reporting fraction three orders of
magnitude too small. The 2013 outbreak was concentrated in a few districts; treating the
entire country as one well-mixed population is wrong, and the fit is saying so through ρ.

**3. The initial condition is pathological.**
For the two 2021 windows, `i0_frac ≈ 0.065` — six and a half percent of the province
already infected in the first modelled week, against a configured sampling range with an
upper bound of 10⁻⁴. Nationally the opposite: `i0_frac ≈ 2 × 10⁻¹³`, with a standard error
larger than the estimate. Both windows begin mid-wave, and the model compensates through
whatever initial condition reproduces the observed level. The parameter is absorbing the
missing early history rather than describing it.

**4. Identifiability is weak where the window is short.**
Only 3 of 10 starts reached the best optimum for the KP climate model, against 9 or 10
for every null model. In Sindh, `corr(β₀, ρ) = −0.85`: transmission and reporting are
nearly exchangeable there, exactly the confounding anticipated in `docs/MODEL.md`.
Standard errors on the climate coefficients exceed the estimates themselves in Sindh
(`a_temp = 0.20 ± 0.60`).

Dispersion of 64 nationally also confirms that Poisson is far too tight for these counts;
confidence intervals must be inflated, as they already are, but a negative binomial
observation model would be the honest choice.

### What this does not mean

It does not mean climate is irrelevant to dengue in Pakistan; the descriptive
cross-correlations in step 02 were real. It means **this parameterisation cannot
distinguish a climate effect from the shape of a single epidemic wave.** With one wave per
window, any smooth covariate that happens to rise and fall near the right time will fit,
and the held-out test is what exposes that.

### Changes this forces, in priority order

1. **Replace the population denominator with an estimated effective susceptible
   population.** The relevant quantity is the population actually at risk in the reporting
   catchment, not the census total. This is the single change most likely to fix ρ.
2. **Constrain ρ to a defensible range** from the under-ascertainment literature rather
   than letting it absorb the denominator error.
3. **Reconsider the temperature term.** A log-linear response over a range spanning both
   sides of the optimum is not mechanistic. Either use a unimodal response peaking near
   29 °C, or drop temperature and retain lagged rainfall alone, which at least has a
   defensible monotone effect on habitat.
4. **Start each window before its wave begins**, so the initial condition describes a
   quiet baseline rather than absorbing unobserved history. The national window already
   does; the 2021 provincial windows do not, and no earlier data exists for them.
5. **Move to a negative binomial observation model**, given dispersion of 6–64.

### Status

The negative result stands until items 1–3 are addressed. If the climate model still
fails out-of-sample after them, that is the paper's finding and it will be reported as
such: that routine surveillance from a single season is insufficient to identify
climate-driven transmission, however well it fits in sample.

---

## 2026-07-27 — Reformulation, and a structural identifiability result

### The reporting fraction cannot be estimated at all

Checked numerically before attempting any fix: halving the population at risk while
doubling the reporting fraction changes the prediction by **exactly zero**, bit for bit.
The state variables are fractions, so population enters only when incidence is converted
to counts, and ρ and N appear solely as the product ρN.

This means the plan's first item — "estimate an effective susceptible population" — was
wrong as stated. Adding N as a free parameter alongside ρ could not have worked; the pair
is unidentifiable by construction, and the optimiser would simply have wandered along the
ridge. The implausible ρ ≈ 10⁻⁴ from the first fit was not a failure of estimation but
the correct answer to a badly posed question.

**Resolution:** fix ρ = 0.05 from the under-ascertainment literature and estimate
`pop_frac`, the fraction of the census population in the transmission catchment. The
identifiable product is unchanged; the error is now attributed to the quantity that
plausibly carries it. Documented in `docs/MODEL.md`; asserted in `tests/test_inference.py`.

### Temperature: unimodal response instead of a log-linear coefficient

The wrong-signed coefficient was a consequence of the parameterisation. A monotone term
fitted to a season whose temperature falls while the epidemic grows must come out
negative, whatever the biology. Replaced with a Brière response using *Aedes aegypti*
thermal limits fixed from the literature (13.35 °C, 40.55 °C; peak at 34.1 °C), entering
as `B(T) ** a_temp`. The exponent makes temperature a nested hypothesis — zero means no
effect, one means the literature response in full — and being non-negative by
construction, it cannot reproduce the artefact.

### The thermal exponent is not identifiable from one season

The recovery test failed after reformulation: `a_temp` came back as 4.96 against a truth
of 1.0. Rather than loosen the tolerance, the likelihood was profiled over it.

| a_temp | 0.5 | 1.0 (truth) | 2.0 | 5.0 | 8.0 |
|---|---|---|---|---|---|
| best logL | −77.23 | −77.15 | −77.03 | **−76.87** | −77.03 |
| β₀ | 0.319 | 0.320 | 0.321 | 0.324 | 0.327 |
| pop_frac | 0.02005 | 0.02005 | 0.02005 | 0.02005 | 0.02006 |

The profile is flat: 0.36 log-likelihood units across a sixteen-fold range of the
parameter, against the ~1.92 units that bound a 95% interval. The estimator is not
failing — the parameter is practically unidentifiable, and the optimiser is picking a
point on a plateau.

The reason is interpretable. β₀ is transmission at the thermal optimum, fixed by the
epidemic's peak; the exponent governs how sharply transmission falls away from the
optimum, which a single wave concentrated near favourable conditions barely observes.

**Consequence for the paper.** β₀ and hence R₀ are recoverable and are what the
conclusions rest on; the shape of the thermal response is not, and no estimate of it will
be reported as if it were. The test suite now states these separately: one test asserts
β₀ recovery across several noise realisations, and the recovery test no longer asserts
anything about `a_temp` — asserting a tolerance on a parameter the data cannot determine
would test the optimiser's tie-breaking rather than the estimator.

### Refit with the reformulated model

| Window | Model | R₀ | Pop. at risk | a_temp | a_rain | ΔAIC | Held-out deviance |
|---|---|---|---|---|---|---|---|
| National 2013 | climate | 1.84 | 0.39 M | 0.026 ± 0.081 | −0.065 ± 0.042 | −388 | +122 |
| National 2013 | null | 1.59 | 0.41 M | — | — | — | — |
| Sindh 2021 | climate | 1.36 | 0.06 M | 0.287 ± 1.256 | +0.241 ± 0.213 | −21 | +86 |
| Sindh 2021 | null | 1.87 | 0.04 M | — | — | — | — |
| KP 2021 | climate | 1.70 | 0.04 M | ~0 | −0.069 | +3 | −16 |
| KP 2021 | null | 1.80 | 0.04 M | — | — | — | — |

**What improved.** The reporting fraction pathology is gone: populations at risk of
40,000 to 390,000 are consistent with outbreaks concentrated in a handful of urban
districts, where 188 million was not. R₀ estimates of **1.4–1.9** now fall squarely inside
the published range for dengue, where before the parameter was absorbing denominator
error.

**What the reformulation revealed.** With the monotone parameterisation removed, the
temperature effect vanishes rather than reversing:

- national: `a_temp = 0.026 ± 0.081` — indistinguishable from zero
- Sindh: `0.287 ± 1.256` — standard error four times the estimate
- KP: `9 × 10⁻¹⁴` — driven to the boundary, and the Jacobian is singular there, so no
  standard errors exist at all

The apparent temperature effect in the first fit was **entirely an artefact of the
log-linear form**. Rainfall fares no better: −0.065, +0.241, −0.069 — inconsistent in sign
across the three windows and never more than one standard error from zero.

Out-of-sample, the climate model still loses nationally (+122, down from +1306) and in
Sindh (+86), and wins only in KP (−16), where AIC simultaneously prefers the null model.
No consistent story.

Dispersion of 164 for the national climate fit confirms the Poisson assumption is badly
violated; the quasi-Poisson inflation of standard errors is doing real work and a
negative binomial observation model remains outstanding.

### Where this leaves the paper

The result is negative and it is now clean, because the confound that produced the
spurious positive has been removed and identified. Three findings stand:

1. **R₀ for Pakistani dengue waves is 1.4–1.9**, estimated from national 2013 and
   provincial 2021 surveillance. New for this setting.
2. **The reporting fraction and the population at risk are structurally
   unidentifiable**, and treating either as free while fixing the other at a census
   total produces parameter estimates that are wrong by orders of magnitude in a way
   that looks like a plausible fit.
3. **Climate covariates cannot be identified from a single season's surveillance.**
   Apparent effects are artefacts of the parameterisation: a monotone term fitted to one
   wave returns whatever sign the covariate's trajectory implies, and a mechanistically
   correct unimodal term returns nothing at all.

That is a more useful paper than the one originally planned, and an honest one. The
question it answers is what can and cannot be learned from routine single-season
surveillance, rather than asserting a climate–dengue relationship the data cannot
support.

### Outstanding

- Negative binomial observation model (dispersion 6–164)
- Sensitivity to ρ, to the rainfall lag, and to the climate location
- Bootstrap confidence intervals
- The PINN comparison, which was the third research question and is untouched

---

## 2026-07-27 (later) — The Poisson assumption was manufacturing the result

Refitted everything with a negative binomial observation model, dispersion profiled by
alternating between the parameters and k until k stops moving.

| Window | Model | R₀ | Pop. at risk | a_temp | a_rain | k | Pearson disp. | AIC | Held-out dev. |
|---|---|---|---|---|---|---|---|---|---|
| National 2013 | climate | 1.51 | 0.44 M | 0.023 ± 0.033 | 0.023 ± 0.031 | 1.39 | 1.0 | 507.7 | 912.3 |
| National 2013 | **null** | 1.55 | 0.41 M | — | — | 1.22 | 0.9 | **506.4** | **706.7** |
| Sindh 2021 | climate | 1.19 | 0.06 M | ~0 ± 1.75 | 0.268 ± 0.374 | 30.5 | 2.2 | **141.0** | 354.1 |
| Sindh 2021 | null | 1.81 | 0.04 M | — | — | 16.0 | 1.5 | 142.8 | **257.8** |
| KP 2021 | climate | 1.46 | 0.06 M | ~0 ± 0.016 | −0.200 ± 0.383 | 22.4 | 2.1 | 113.4 | 20.2 |
| KP 2021 | **null** | 1.70 | 0.04 M | — | — | 19.1 | 1.4 | **111.2** | **0.4** |

### The headline

**The climate model's entire in-sample advantage was an artefact of assuming Poisson.**

Under Poisson the national ΔAIC was −388, an overwhelming preference for climate forcing.
Under the negative binomial it is **+1.3** — a preference for the null model. Sindh moves
from −21 to −1.8, KP from +3 to +2.2. Nothing about the data or the transmission model
changed; only the assumed variance did.

The mechanism is straightforward. Poisson asserts that variance equals mean, so a week
with 2,000 cases is treated as accurate to about ±45. The observed scatter is far larger
than that, and the only way to reduce a residual the model regards as enormous is to bend
the transmission curve toward individual weeks. The climate covariates supply exactly the
flexibility to do so. Once the negative binomial admits that a 2,000-case week carries
variance of order 2,000 + 2,000²/1.4, chasing it stops paying, and the climate terms earn
nothing.

Pearson dispersion confirms the fix: 164 under Poisson, **1.0** under the negative
binomial.

### Consolidated result

With a defensible observation model, climate forcing is not supported by any criterion in
any window:

- AIC prefers the null model in two of three windows and is within 2 units in the third
- Held-out deviance prefers the null model in **all three**
- Every climate coefficient is within one standard error of zero, and both thermal
  exponents are driven to the boundary at ~10⁻⁹ and ~10⁻¹³

R₀ from the null model is stable across windows and separated by eight years:
**1.55, 1.81, 1.70**.

### A new identifiability warning

Khyber Pakhtunkhwa, null model: `corr(β₀, pop_frac) = −0.991`. Transmission and catchment
size are almost perfectly exchangeable over a fifteen-week window, so that R₀ estimate is
far less secure than its standard error of ±0.012 on β₀ suggests. Sindh is better at
−0.71, the national window much better at +0.12 — longer series identify the pair, short
ones do not.

This is precisely what asymptotic standard errors cannot express, and it is why the
bootstrap is the next step rather than an optional refinement.

### Outstanding after this round

- Bootstrap confidence intervals — now clearly necessary, not optional
- Sensitivity to ρ, rainfall lag, climate location
- The PINN comparison

---

## 2026-07-27 (later still) — The PINN fits better and is wrong

`scripts/06_pinn_compare.py`. Both estimators optimise the same negative binomial
likelihood on the same data, so any difference is attributable to the method.

### Synthetic data, where the truth is known

| | R₀ | pop_frac | logL | Time | Forward solves |
|---|---|---|---|---|---|
| Truth | 1.700 | 0.0025 | — | — | — |
| Classical | **1.752** (3% error) | 0.00186 | −160.3 | 678 s | 19,973 |
| PINN | **0.022** (99% error) | 0.0996 (40× high) | −432.7 | 127 s | 0 |

The PINN does not recover parameters it was handed. Its trajectory is instructive:
β₀ moved from 0.30 to a reasonable 0.19 by 5,000 epochs, then collapsed through
0.13, 0.034, 0.007 to 0.002 while `pop_frac` inflated from 0.005 to 0.0996. It slid
along the confounding ridge between transmission and catchment size, and nothing
stopped it.

### Real data — where it looks like a success

| Window | Classical R₀ | Classical logL | PINN R₀ | PINN pop_frac | PINN logL |
|---|---|---|---|---|---|
| National 2013 | 1.599 | −311.9 | 1.083 | **0.865** | **−309.8** |
| Sindh 2021 | 1.025 | −100.0 | **0.024** | **0.655** | **−83.8** |
| KP 2021 | 1.669 | −67.2 | **0.267** | **0.724** | **−48.1** |

**The PINN attains a better likelihood in all three windows** — by 2, 16 and 19 units
— and every estimate it produces is impossible:

- R₀ of 0.024 and 0.267 for epidemics that demonstrably occurred. Below one, no
  outbreak can grow at all.
- 65 to 87 per cent of a province's entire population in the transmission
  catchment, against 0.1–0.4 per cent from the classical fit.
- Thermal exponents of 15.2 and 23.8, which drive the temperature response to
  numerical zero for all but the hottest weeks.

### Why, and why it matters

The physics residual never approached zero: it finished at 2.6, 7.8 and 20.5. The
network therefore never represented a solution of the model, and a curve that does
not solve the equations is not constrained by them. `figures/05_pinn_comparison.png`
shows the consequence directly — the PINN traces individual weekly bars, following
Sindh's secondary January rise and KP's dip at the start of November, neither of
which any single epidemic wave can produce. The classical curve is smooth because it
has no choice: it *is* a solution.

So the better likelihood is not evidence of a better fit to the epidemic. It is
evidence that the constraint was relaxed. And once relaxed, the parameters are
free — which is exactly why they came out impossible.

This sharpens the failure mode into something more specific than "PINNs can be
inaccurate":

> Approximate enforcement of the physics dissolves the identifiability that makes
> the inverse problem solvable. The classical estimator recovers interpretable
> parameters *because* it is forced to solve the model exactly; the PINN's freedom
> to satisfy the equations only approximately is precisely what destroys the
> parameter estimates, while improving the very likelihood one would use to judge
> the fit.

A practitioner comparing the two on likelihood alone would choose the PINN and
report an R₀ below one for an epidemic that infected thousands.

### Before concluding

One configuration failing is not evidence about a method. `scripts/07_pinn_tuning.py`
varies the choices most plausibly responsible — loss weighting across three orders of
magnitude, two-stage training with the parameters frozen initially, denser
collocation, longer runs, a wider network — and records every outcome on the same
synthetic data. The claim in the paper will rest on that table, not on the single
attempt above.

---

## 2026-07-27 — Bootstrap: the asymptotic intervals were wrong by two orders of magnitude

`scripts/05_bootstrap.py`, 200 parametric replicates per window and model, simulated from
each fitted negative binomial and refitted from scratch.

| Window | Model | R₀ | 95% interval | Population at risk | 95% interval |
|---|---|---|---|---|---|
| National 2013 | climate | 1.599 | [1.50, 1.71] | 0.81 M | [0.43, 1.78] |
| National 2013 | null | 1.372 | [1.32, 1.44] | 0.49 M | [0.25, 0.82] |
| Sindh 2021 | climate | 1.025 | **[0.78, 2.45]** | 0.12 M | **[0.04, 47.9]** |
| Sindh 2021 | null | 1.420 | [1.09, 1.61] | 0.058 M | [0.041, 0.119] |
| KP 2021 | climate | 1.669 | [1.46, 1.88] | 0.045 M | **[0.041, 35.5]** |
| KP 2021 | null | 1.706 | [1.61, 1.80] | 0.044 M | [0.0396, 0.0500] |

### Why the bootstrap was necessary rather than decorative

Comparing its spread against the Jacobian-based standard errors reported earlier:

- Sindh, null model, population at risk: the bootstrap standard deviation is **417 times**
  the asymptotic one. The asymptotic interval was not merely optimistic; it was
  meaningless.
- Sindh, climate model, same parameter: 25 times.
- Several others run 2–3 times too *wide* in the opposite direction.

The asymptotic errors assume a locally quadratic likelihood and well-identified
parameters. Where `corr(β₀, pop_frac)` reached −0.99, neither holds, and the numbers those
formulae produce should not have been believed. They are retained in the tables only for
this comparison.

### What the intervals say

**R₀ is identified; the catchment size often is not.** The national R₀ interval spans
0.21, and KP's null-model interval spans 0.19 — usable. But the population at risk under
the climate model spans a factor of 400 in Sindh and 800 in KP. Adding climate parameters
to a fifteen-week window does not merely fail to help; it destroys the identifiability of
a parameter the null model pins down to ±6%.

**Sindh's climate-model R₀ interval includes values below one.** An outbreak that infected
over two thousand people cannot have had R₀ < 1. That the interval admits it is a further
sign the climate model is over-parameterised for these short windows.

### Model uncertainty exceeds parameter uncertainty

`figures/04_bootstrap_R0.png` shows the two models' bootstrap distributions side by side.
For the national window they **barely overlap**: the climate model concentrates on
[1.50, 1.71] and the null model on [1.32, 1.44]. Choosing between two models that the
data cannot distinguish — ΔAIC of 1.3 — shifts R₀ by roughly 0.2, more than the sampling
uncertainty within either.

A confidence interval conditional on one model therefore understates what is actually
unknown. The paper should report the range spanning both models, approximately
**R₀ = 1.3–1.7 nationally**, rather than a narrow interval from whichever model is
preferred on the day. In Khyber Pakhtunkhwa the two distributions coincide closely, so
the problem is specific to windows where the climate parameters are doing something,
even if that something is not identifiable.

### An inconsistency to resolve before writing up

The bootstrap fits the full window; step 03 fits the first 75% and holds out the rest.
The national climate estimates differ accordingly: R₀ 1.599 against 1.51, and
`a_temp = 0.215` with an interval of [0.12, 0.34] that **excludes zero**, against
`0.023 ± 0.033` on the training subset, which includes it.

So whether the temperature term is distinguishable from zero depends on whether the final
nineteen weeks are included. That instability is itself evidence for the overall
conclusion, but the paper must report headline parameters from one fit — the full window,
using all the data — and reserve the split strictly for the model comparison. Reporting
whichever is more convenient would be indefensible.

---

## 2026-07-27 — Sensitivity, and a label that deleted itself

`scripts/04_sensitivity.py`, 74 fits.

### The reporting fraction: exactly as the algebra requires

| ρ | 0.02 | 0.05 | 0.10 | 0.20 |
|---|---|---|---|---|
| National R₀ | 1.5990 | 1.5990 | 1.5990 | 1.5990 |
| Population at risk | 2.018 M | 0.807 M | 0.404 M | 0.202 M |
| ρ × population | 0.04036 | 0.04037 | 0.04036 | 0.04036 |

R₀ is identical to four decimal places across a tenfold range of ρ, and the product is
constant to 0.01% nationally, 0.1% in Sindh, 0.2% in Khyber Pakhtunkhwa. The structural
identifiability result is confirmed numerically, not merely argued: **the assumed
reporting fraction rescales the estimated catchment and touches nothing else.** Every
conclusion about transmission is therefore free of it.

### The rainfall lag and the climate location: not free at all

| Window | R₀ spread over rainfall lags 3–7 wk | R₀ spread over climate locations |
|---|---|---|
| National 2013 | 1.55–1.62 (0.06) | **1.36–1.68 (0.32)** |
| Sindh 2021 | **0.85–1.69 (0.83)** | **1.02–1.47 (0.45)** |
| KP 2021 | 1.61–1.77 (0.16) | 1.64–1.68 (0.05) |

Nationally, which city stands in for the country's climate moves R₀ by 0.32 — larger than
the bootstrap sampling interval of 0.21 and comparable to the shift between models.
Peshawar gives 1.36 with a thermal exponent driven to zero; Karachi gives 1.68 with an
exponent of 0.59. Sindh is worse: the rainfall lag alone moves R₀ from 0.85 to 1.69,
across the epidemic threshold.

These are choices, not measurements. An arbitrary decision that moves the headline
estimate by more than the sampling uncertainty must be reported as a source of
uncertainty, and the constant-transmission model — which uses no climate data and is
therefore immune to both — is correspondingly more trustworthy.

### The label that deleted itself

The constant-transmission model was labelled `null` in the output tables. Pandas reads
that string back as `NaN`. Re-importing the sensitivity table for analysis therefore
dropped **every one of its 37 rows** — hours of computation, present in the file, absent
from anything that opened it normally, with no error anywhere.

It was caught only because a consistency check expected the constant model's estimates to
be invariant under the rainfall lag (it ignores rainfall) and there were no rows to check.
Had that check not been written, the sensitivity analysis would have been reported for the
climate model alone, and the comparison that matters — that the constant model is immune
to these choices — would have been missed entirely.

The label is now `constant`, and `tests/test_outputs.py` scans every result table for
strings the CSV reader silently converts to missing, verifies row counts survive the round
trip, and pins the expected model labels. All prior tables were regenerated.

`scripts/07_pinn_tuning.py`. All on the same synthetic data, truth R₀ = 1.700,
pop_frac = 0.0025, a_rain = 0.100. Classical estimator on identical data: R₀ error 3.1%.

| Configuration | R₀ | error | pop_frac | a_rain | physics | logL |
|---|---|---|---|---|---|---|
| A baseline (weight 200) | 0.019 | 98.9% | 0.130 | 0.47 | 28.7 | −416 |
| B weight 2,000 | 0.945 | 44.4% | 0.00049 | −0.22 | 60.8 | −330 |
| **B weight 20,000** | **1.702** | **0.1%** | 0.00019 | 0.24 | 82.6 | −433 |
| B weight 200,000 | 0.186 | 89.0% | 0.00056 | 0.60 | 84.2 | −218 |
| C dense grid, 80k epochs | 1.885 | 10.9% | 0.982 | 0.02 | 13.3 | −389 |
| D two-stage | 0.004 | 99.8% | 0.0046 | 0.78 | 78.4 | −429 |
| E wide network | 1.507 | 11.3% | 0.013 | 0.75 | 52.2 | −396 |

Weighting the data term at 20,000 recovers R₀ to **0.1%** — better than the classical
estimator's 3.1%. Taken alone, that is a headline: a PINN matching and beating classical
inversion at a fraction of the compute.

It does not survive contact with a second random seed.

| Same configuration, seed | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| R₀ | 1.702 | 1.565 | 0.848 | 1.001 | 1.953 |
| error | **0.1%** | 8.0% | **50.1%** | 41.1% | 14.9% |

Nothing changed but the network initialisation. The estimate ranges from exact to wrong
by half, with a median error near 15%. The apparent success was a favourable draw, and on
real data — where there is no truth to check against — a practitioner has no way to learn
which draw they received.

### Three further observations

**No configuration recovers the other parameters.** The population at risk, true value
0.0025, came back between 0.00019 and 0.982 — four orders of magnitude. The rainfall
coefficient, true value 0.100, came back between −0.946 and +0.781, with no configuration
even getting the sign reliably. Only R₀, and only sometimes.

**Likelihood is anti-correlated with accuracy.** The configuration that recovered R₀ to
0.1% has the *worst* log-likelihood in the table (−433); the configurations with the best
likelihoods (−218, −225, −232) have R₀ errors of 89%, 50% and 8%. Model selection on fit
quality would systematically pick the wrong configuration.

**A smaller physics residual does not help either.** The lowest residual, 13.3 from the
dense-grid run, came with a population at risk of 0.982 — the entire country in the
transmission catchment. Neither diagnostic available to a practitioner without ground
truth points toward the right answer.

### Conclusion for the paper

Across eleven configurations spanning loss weighting, training schedule, network capacity,
collocation density and initialisation, the PINN does not reliably recover parameters that
the classical estimator recovers to 3% every time. The configuration that appeared to beat
it was seed-dependent, and neither the likelihood nor the physics residual identifies the
good runs.

The mechanism is the one identified earlier: with the equations satisfied only
approximately, the trajectory is not confined to the model's solution manifold, and the
parameters are consequently unconstrained. The classical estimator's cost — thousands of
forward solves — is what buys the constraint, and on this problem the constraint is the
thing of value.

This does not generalise to problems where a forward solve costs hours rather than
milliseconds; that case is untested here and the paper will say so. What it does establish
is that the standard justification for PINNs on inverse problems — comparable accuracy at
lower cost — does not hold on a seven-state epidemic model with weekly surveillance data,
and that the failure is invisible to every diagnostic available without ground truth.

---

## 2026-07-27 — Aggregation bias: the second research question, finally tested

`scripts/09_aggregation_bias.py`. The README has posed this question since the project
began and no analysis had answered it. The 2021 season allows a direct test: districts and
their sum report over the same weeks, so the same epidemic can be fitted at two spatial
resolutions and the answers compared.

The comparison is unusually clean because R₀ does not depend on the assumed population at
all. District census figures affect only the reported catchment, never the quantity being
compared.

| Province | Unit | Cases | R₀ |
|---|---|---|---|
| Sindh | Hyderabad | 929 | 2.21 |
| Sindh | Tharparkar | 1,006 | 1.36 |
| Sindh | Karachi Malir | 168 | 1.65 |
| Sindh | Karachi East | 18 | 1.09 |
| Sindh | **aggregate of the four** | 2,121 | **1.54** |
| KP | Haripur | 1,323 | 1.25 |
| KP | Lakki Marawat | 240 | 1.96 |
| KP | **aggregate of the two** | 1,563 | **1.14** |

**Aggregation biases R₀ downward, and can place the aggregate outside the range of its own
constituents.** In Sindh the bias is mild — −2.2% against the district mean, and the
aggregate sits inside the district range. In Khyber Pakhtunkhwa it is −28.8%, and the
aggregate estimate of 1.14 lies **below every district it is composed of**, 1.25 and 1.96.

The mechanism is visible in `figures/07_aggregation_bias.png`. Haripur peaks in late
October; Lakki Marawat peaks in early October and again in late November. Their sum is
broader and flatter than either, and a homogeneously mixing model can only read a broader
curve as slower transmission. Sindh's districts peak closer together, which is why its
bias is smaller.

### The consequence for the headline estimate

The national 2013 series aggregates outbreaks across the whole country. If summing two
districts of one province depresses R₀ by 29%, summing dozens across four provinces
depresses it by at least as much. **The national R₀ of roughly 1.5 should be read as a
lower bound on transmission in any constituent locality**, not as a characteristic value.
The paper must say so explicitly, because the number invites the opposite reading.

### Two bugs found while doing it, both surfaced by the data rather than by a test

**Duplicate district names.** The OpenDengue Pakistan subset carries the same places under
variant spellings: `HYDERABAD`/`HYDERBAD`, `LAKKI MARAWAT`/`LAKKIMARWAT`/`LAKIMARWAT`,
`KAMBER`/`QAMBER`, `KARACHI WEST`/`KARWEST`. The stray `HYDERBAD` record covers one week
and carries 107 cases; because the comparison requires a window shared by every unit, that
single row collapsed the shared window to one week and produced a silently empty analysis.
Variants are now normalised and listed in the script.

**Incidence over non-contiguous weeks.** `weekly_incidence` differenced consecutive
cumulative values, which equals each observation's own seven days only while the observed
weeks are contiguous. Districts do not report on aligned weeks, so intersecting them
produces gaps — and across a three-week gap the old code attributed three weeks of
incidence to one reported week, over-counting by up to **184×** in the test case. The main
study windows are contiguous by construction and their results are unaffected, which is
exactly why the bug survived: every earlier analysis exercised the one case in which the
two formulations agree. Each observation now spans its own seven days, and
`tests/test_models.py` checks a deliberately gapped series.

Correcting it moved the Khyber Pakhtunkhwa aggregation bias from −16.9% to −28.8%, so the
conclusion strengthened rather than weakened.

The figure needed fixing too: the aggregate was drawn last and thick, so it sat exactly on
top of whichever district dominated, making the sum appear to exceed its parts. That
appearance is what prompted the check that found the incidence bug.

---

## 2026-07-27 — Time-marching, and twenty-one runs pooled

### The remaining PINN remedy does not help

`scripts/10_pinn_timemarching.py` implements the progressive time-window curriculum that
step 07's docstring promised and did not deliver. Note first that the standard
prescription does not transfer unchanged: Krishnapriyan et al. recommend a separate
network per window for *forward* problems, but here the coefficients are what is being
estimated and they are global, so a per-window network with its own parameters would be
estimating a different epidemic in each window. The transferable part is the curriculum.

| Scheme | R₀ error across five seeds | median |
|---|---|---|
| Global (step 07 best) | 0.1% – 50.1% | 14.9% |
| Expanding window | 7.9% – 54.5% | 49.2% |
| Expanding window + optimiser restart | 10.2% – 57.6% | 21.3% |
| Classical estimator | 3.1%, deterministic | 3.1% |

Neither variant improves on global training, and both remain an order of magnitude worse
than the classical estimator.

### Pooling every run: no diagnostic works, and one is actively misleading

`scripts/12_pinn_diagnostics.py` pools all 21 PINN fits — 11 from step 07, 10 from
step 10 — on identical synthetic data with a known truth.

**One run in twenty-one beat the classical estimator.** Median R₀ error 41.1% against
3.1%; estimated population at risk ranged from 0.08× to 398× the true value.

The question that matters for a practitioner is whether anything computable from a
finished run identifies the good ones. Correlations with log R₀ error:

| Diagnostic | Correlation | Verdict |
|---|---|---|
| Log-likelihood on the data | −0.06 | uninformative |
| Physics residual | −0.03 | uninformative |
| Implausibility of the estimated population | −0.04 | uninformative |

None of them. A practitioner with real data has no basis for trusting or rejecting any
particular run.

**And the physics residual is worse than uninformative about the population.**

    corr(log physics residual, error in estimated population) = −0.89

| | Runs | Median population | Median R₀ error |
|---|---|---|---|
| Physics residual < 10 | 6 | **318× the truth** | 31.9% |
| Physics residual ≥ 10 | 15 | 1.3× the truth | 41.1% |

The runs that satisfy the governing equations best are the ones that place most of the
country in the transmission catchment. The mechanism is interpretable: with a very large
population the epidemic barely depletes susceptibles, the nonlinear SI coupling becomes
nearly linear, and a smooth network output can track that regime far more easily than one
in which saturation matters. Driving the residual down therefore pushes the solution into
a regime that fits the equations well and describes the data badly.

This is the sharpest form of the project's central point. It is not that PINNs are
inaccurate. It is that **the diagnostic conventionally reported as evidence a PINN has
worked — a small physics residual — moves in the opposite direction to the accuracy of the
parameters it is being used to justify**, and that no other available diagnostic
compensates.

---

## 2026-07-27 — Profile likelihoods, and a growth rate the model cannot explain

### Practical identifiability

`scripts/11_identifiability.py` profiles each parameter: fix it across a grid, refit
everything else, and look at how fast the likelihood falls away. A profile staying within
1.92 log-likelihood units of the optimum across the whole grid is not bounded by the data.

Most parameters are identified. Three are not:

- **`i0_frac`, national window, both models** — unbounded across a 25-fold range. The
  initial infected fraction is absorbing unobserved history before the series begins,
  exactly as suspected in the first analysis round.
- **`a_rain` in Sindh under the climate model** — 20-fold range within the threshold.
- **`pop_frac` in Sindh under the climate model** — 12-fold, against a 2-fold range under
  the constant model on the same data. Adding climate parameters to a nineteen-week
  window destroys the identifiability of a parameter the simpler model determines.

A caveat on the grids: for the sharply identified parameters the interval collapses onto a
single grid point, meaning the profile rises past the threshold within one step of about
26%. The profiles establish *whether* a parameter is bounded, not a precise interval; the
bootstrap supplies the intervals.

### The epidemics grow faster than the fitted R₀ permits

An independent check, using no compartmental fit at all: estimate the exponential growth
rate r of the early epidemic and bound R₀ from it. The bound depends on how the generation
interval is distributed, and the extremes are far apart — R₀ = exp(rT) for a fixed
interval, R₀ = 1 + rT for an exponentially distributed one. Quoting only the first, as is
common, overstates R₀ badly when rT is large, so both are reported.

| Window | r (per day) | R² | R₀ from growth | Fitted R₀ | |
|---|---|---|---|---|---|
| National 2013 | 0.0956 | 0.94 | 2.91 – 6.77 | **1.37** | below the lower bound |
| Sindh 2021 | 0.0534 | 0.89 | 2.07 – 2.91 | **1.42** | below the lower bound |
| KP 2021 | 0.0397 | 0.88 | 1.79 – 2.21 | **1.71** | marginally below |

**In every window the fitted R₀ is smaller than the observed growth rate can support**,
and by a factor of two nationally even against the most conservative bound.

This is not a fitting failure; it is the model telling us it cannot describe both features
of the data at once. A homogeneously mixing model ties the growth rate to the final size
through a single R₀. The fitted values imply attack rates of 49%, 53% and 69% of the
catchment; the growth-rate bounds would imply 100%, 93% and 85%. Faced with an epidemic
that rises quickly and then infects fewer people than that rise implies, the likelihood —
dominated by the many weeks around the peak and decline — settles on the final size and
under-states transmission.

### This connects directly to the aggregation result

The two findings are the same phenomenon seen from different sides. These series are sums
of asynchronous local outbreaks; the sum is broader and flatter than any component, so a
homogeneous fit reads it as slower transmission. The growth-rate check quantifies the cost
independently of any fit, and reaches the same verdict: **the fitted R₀ is a lower bound,
depressed by spatial aggregation, and should not be read as the transmission rate of any
actual population.**

That, rather than a specific value, is what this analysis can honestly conclude about R₀
in Pakistani dengue from data of this kind.

---

## 2026-07-27 — A third way the climate conclusion is manufactured

Separating the reported estimate from the held-out scoring fit — the full window for the
former, the first 75% for the latter — exposed something the earlier runs concealed by
using one fit for both purposes.

| Window | ΔAIC, full window | ΔAIC, training fit | Held-out deviance |
|---|---|---|---|
| National 2013 | **−31.8** (climate) | **+1.3** (constant) | +205.6 (constant) |
| Sindh 2021 | **+2.9** (constant) | **−1.8** (climate) | +96.3 (constant) |
| KP 2021 | +3.3 (constant) | +2.2 (constant) | +19.8 (constant) |

**The AIC verdict flips in two of three windows**, and in opposite directions, on nothing
more than whether the last quarter of the series is included. Nationally the full window
gives a decisive −31.8 in favour of climate forcing, where the training subset gives +1.3
against it; Sindh does the reverse.

Both are defensible analyses. Fitting all the data is the obvious choice for an estimate;
holding out a tail is the obvious choice for validation. Neither is a mistake, and a paper
reporting either alone would be reporting an artefact of that choice.

### Three choices, three reversals, one consistent answer

This is now the third routine decision found to reverse the in-sample conclusion:

1. **The functional form of the temperature term.** A monotone coefficient fitted to a
   season whose temperature falls while cases rise returns a large negative effect; a
   unimodal response derived from vector thermal biology returns nothing.
2. **The observation model.** Poisson gives a national ΔAIC of −388 in favour of climate;
   a negative binomial, with dispersion of 164 in the data, gives +1.3 against.
3. **How much of the series is fitted.** Full window −31.8; first 75%, +1.3.

Each is a choice an analyst makes without much thought, none is wrong, and each moves the
conclusion across the decision boundary.

**Only one criterion is stable.** Held-out deviance prefers the constant-transmission
model in all three windows, under both observation models, at every rainfall lag and every
climate location tested. Out-of-sample prediction is the only thing here that does not
depend on how the question was set up — which is the paper's central methodological point,
now supported by three independent reversals rather than argued from one.

### Consistency restored

The full-window estimates now match the bootstrap point estimates exactly (national
climate R₀ 1.60 against 1.599; constant 1.37 against 1.372; and likewise for the
provincial windows), and the profile likelihoods are computed from the same fits. The
earlier mismatch — flagged in the bootstrap entry above — is resolved: parameters are
reported from the full window, and the split is used for nothing but scoring.

---

## 2026-07-27 — Which fit is reported, and a third reversal

Until now step 03 fitted the first 75% of each window and reported those parameters, while
the bootstrap and the profile likelihoods used the full window. Two different estimates
were circulating for the same quantity, and the difference was not cosmetic.

The script now performs both fits explicitly and labels them:

* the **full-window fit** is the estimate, uses all the data, and matches what the
  bootstrap and profiles are computed from — this is what the paper quotes;
* the **training fit** exists only to score the held-out weeks, and its parameters are
  never reported.

### The in-sample verdict reverses with the fitting window

| Window | ΔAIC, full window | ΔAIC, first 75% | Held-out deviance |
|---|---|---|---|
| National 2013 | **−31.8** (favours climate) | **+1.3** (favours constant) | +205.6 (favours constant) |
| Sindh 2021 | **+2.9** (favours constant) | **−1.8** (favours climate) | +96.3 (favours constant) |
| KP 2021 | +3.3 (constant) | +2.2 (constant) | +19.8 (constant) |

Two of the three windows change their answer, and in opposite directions, purely from
whether the last quarter of the series is included. Neither analysis is wrong; both are
defensible choices an author might make without comment.

That is now the **third** routine decision found to reverse the in-sample conclusion about
climate forcing, alongside the functional form of the temperature term and the observation
model. Three separate choices, none of them visible in a final model specification, each
capable of flipping the headline.

### What survives all of them

Held-out deviance prefers the constant-transmission model in all three windows, under both
observation models, at every rainfall lag and every climate location tested, on both
fitting windows. It is the only criterion in this study that does not depend on how the
question was set up — which is the practical recommendation the paper can make.

Reported estimates for the record, full-window fits:

| Window | Model | R₀ | Population at risk | a_temp | a_rain | k |
|---|---|---|---|---|---|---|
| National 2013 | climate | 1.60 | 0.81 M | 0.215 | −0.134 | 0.57 |
| National 2013 | constant | 1.37 | 0.49 M | — | — | 0.33 |
| Sindh 2021 | climate | 1.03 | 0.12 M | ~0 | 0.158 | 2.96 |
| Sindh 2021 | constant | 1.42 | 0.06 M | — | — | 2.78 |
| KP 2021 | climate | 1.67 | 0.04 M | ~0 | −0.066 | 28.1 |
| KP 2021 | constant | 1.71 | 0.04 M | — | — | 26.6 |

---

## Results phase closed — 2026-07-27

14 pipeline scripts, 6 modules, 82 tests, 13 tables, 10 figures, 10 dated entries above.
No analysis gap remains that I can identify. Every number the paper will quote is
reproducible from `scripts/` under a fixed seed, and every choice that could have gone
another way is recorded here with what happened when it did.

---

## 2026-07-27 — Testing the explanation, not just the observation

Two results pointed at spatial heterogeneity: aggregation depresses the fitted
transmission rate, and the observed growth rates exceed what the fitted rates permit.
Heterogeneity was offered as the reason for both, and offering a reason is not testing it.

`scripts/13_heterogeneity_test.py` fits a two-patch model to the Khyber Pakhtunkhwa
aggregate — the case with the largest bias — using the aggregate series alone. It is told
nothing about which districts contributed or when they peaked. If heterogeneity is the
explanation, three things should follow.

| | One patch | Two patches |
|---|---|---|
| R₀ | 1.142 | **2.534** and **1.723** |
| Catchment | 6.3% of the province | 0.7% and 1.7% |
| Second patch offset | — | **+46 days** |
| AIC | 188.4 | **179.7** |

1. **Better fit** — ΔAIC = −8.6 for two patches over one, on three extra parameters.
2. **Rates recovered** — the two-patch estimates of 1.72 and 2.53 sit at and above the
   separately fitted district values of 1.25 and 1.96, where the single-patch fit sat
   *below both* at 1.14.
3. **Growth reconciled** — the faster patch at 2.53 clears the growth-rate lower bound of
   1.79 that the single-patch fit violated.

The model also recovered a 46-day offset between the two patches without being given any
timing information. Haripur peaks in late October and Lakki Marawat has its secondary rise
in late November, so the asynchrony it inferred is the asynchrony that is there.

**The explanation survives testing.** The aggregation bias is not merely an observed
artefact of unknown origin; it is what a homogeneous model does when the truth is two
epidemics offset in time, and allowing that structure recovers what aggregation hid.

### A caveat that must be stated

The two-patch estimates do not reproduce the district values exactly — 1.72 and 2.53
against 1.25 and 1.96. They should not be expected to. An aggregate series does not
uniquely determine its decomposition, and the model recovers *a* two-component structure
consistent with the data rather than *the* true one. What the test establishes is that a
heterogeneous structure explains the discrepancy and a homogeneous one cannot, not that
the specific split is correct.

### And a near-miss worth recording

The first version of this test also passed all three predictions — and its fitted first
patch placed **279 million people** in a province of 1.9 million. The population fraction
had been kept positive by a log transform but never bounded above one, so the optimiser
produced a sub-critical patch acting as a slowly decaying background rather than an
epidemic. The fit was better and the verdict was "YES", and both were worthless.

Bounding the fraction to (0, 1) by a logit and refitting gave the result above, which
happens to be stronger. But the lesson is the near-miss, not the outcome: a model
comparison had been won on a parameter set with no physical meaning, and nothing in the
AIC, the likelihood or the three pre-registered checks would have caught it. Only reading
the fitted values did.

---

## Results phase closed — 2026-07-27

15 pipeline scripts, 7 modules, **90 tests passing**, 14 tables, 11 figures, and 12 dated
entries above.

No analysis gap remains that I can identify. Every number the paper will quote is
reproducible from `scripts/` under a fixed seed, and every choice that could have gone
another way is recorded here together with what happened when it did.

### The eight results

1. **R₀ = 1.4–1.7** across three windows, two of them eight years apart — and it is a
   **lower bound**, for two independent reasons (3 and 4).
2. **The reporting fraction and the population at risk are structurally unidentifiable.**
   Only their product. Varying ρ tenfold rescales the population estimate as 1/ρ and
   leaves R₀ at 1.5990 to four decimal places.
3. **Spatial aggregation depresses R₀** — by 2% in Sindh and 29% in Khyber Pakhtunkhwa,
   where the aggregate estimate falls below both constituent districts.
4. **The epidemics grow faster than the fitted R₀ allows.** A model-free bound from the
   early growth rate exceeds the fitted value in all three windows.
5. **Heterogeneity explains 3 and 4, tested not asserted.** A two-patch model fitted to
   the aggregate alone beats one patch by ΔAIC −8.6, recovers rates above both districts,
   clears the growth bound, and infers the real 46-day asynchrony unprompted.
6. **Climate forcing is unsupported, and three routine analysis choices each reverse the
   in-sample verdict** — the functional form of the temperature term, the observation
   model, and how much of the series is fitted. None is visible in a final specification.
7. **Held-out prediction is the only stable criterion.** It prefers constant transmission
   in all three windows under every variation tested.
8. **PINNs do not replace classical inversion here.** One of 21 configurations beat it;
   no available diagnostic identifies the good runs, and the physics residual is
   anti-correlated with the accuracy of the population estimate at −0.89.

### What the log is for

Five bugs in my own work appear above — three caught by tests, two by looking at data or
a figure that seemed wrong. Two of them (a category label parsed as missing, an incidence
window that silently over-counted by 184×) produced no error and no warning: they returned
a plausible number or an empty table. A sixth near-miss passed every pre-registered check
while placing 279 million people in a province of 1.9 million.

That is the argument for keeping a log of this kind rather than only a final method
section. The choices that decide a result are mostly invisible in the result.

---

## 2026-07-27 — The global study: 236 outbreaks, 5,664 fits

`scripts/14_global_windows.py` → `16_global_robustness.py` → `17_robustness_analysis.py`.

The Pakistani result — that three routine analysis choices each reverse the in-sample
verdict on climate-driven transmission — was an anecdote from three windows. This repeats
the test on every usable outbreak in OpenDengue.

**Sample.** 237 windows from 34 countries, extracted by a documented rule: an unbroken
weekly run, split into waves outward from each peak, 30–110 weeks, at least 300 cases and
a peak at least 2.5× the median week with the peak interior to the window. The
peak-prominence threshold is reported across a range of values before one is applied, so
the sample size is not tuned. 236 windows completed the full factorial; one Brazilian
window failed to fit at all and is excluded.

**Design.** Each window fitted under all 24 combinations of four choices — observation
model (Poisson, negative binomial), temperature term (unimodal Brière, conventional
log-linear), rainfall lag (3, 5, 7 weeks), fraction of the series fitted (75%, 100%) —
against both the climate-forced and constant-transmission models. 5,664 fits, 150 minutes.

### Result 1: the verdict is unstable in 83% of outbreaks

| | Windows | Share |
|---|---|---|
| Always favours climate forcing | 41 | 17.4% |
| Never favours climate forcing | 0 | 0.0% |
| **Verdict changes with the choices** | **195** | **82.6%** |

In more than four outbreaks in five, whether the data support climate-driven transmission
depends on decisions the analyst makes before seeing the answer. Among unstable windows
the median split is 18 of 24 combinations in favour — not a marginal wobble around a
boundary, but a substantial disagreement inside a single dataset.

### Result 2: the observation model is the dominant lever

Paired within window, so that differences in series length or case count cannot drive it:

| Choice varied | Fits whose verdict flips | P(climate wins) before → after |
|---|---|---|
| **Negative binomial → Poisson** | **36.4%** | 0.565 → 0.921 |
| 75% → 100% of the series | 21.9% | 0.700 → 0.786 |
| Rainfall lag 3 → 7 weeks | 15.5% | 0.745 → 0.742 |
| Brière → log-linear temperature | 11.9% | 0.694 → 0.791 |
| Rainfall lag 3 → 5 weeks | 10.3% | 0.745 → 0.740 |

Assuming Poisson counts raises the probability of endorsing climate forcing from 57% to
92%. It is the single most consequential decision in the analysis, it is rarely justified
in published work, and it is not a decision about dengue.

Note the rainfall lag: it flips individual verdicts 10–16% of the time while barely moving
the average. That is noise injection rather than a systematic effect — an arbitrary choice
adding variance to a conclusion without shifting it.

### Result 3: out-of-sample validation is not the remedy I claimed

This corrects an earlier entry. On the Pakistani windows, held-out deviance preferred the
constant model in every configuration tested, and the log recorded that "out-of-sample
prediction is the only thing here that does not depend on how the question was set up."
**On the global sample that is false.**

| | Unstable within window |
|---|---|
| In-sample (AIC) | 75.7% |
| Out-of-sample (held-out deviance) | 76.6% |

Held-out validation is *equally* unstable. What it does do is halve the endorsement rate:
climate forcing is preferred in **69.9%** of fits in-sample and **35.0%** of the same fits
out-of-sample. So validation substantially reduces the false endorsement of climate
forcing, but it does not make the conclusion independent of the analyst.

That is a weaker and more useful claim than the one the Pakistani data suggested, and it
had to be checked on more than three windows to find out.

### Result 4: the conventional temperature term reports an arbitrary sign

| Parameterisation | Fits | Mean | Median | Negative in |
|---|---|---|---|---|
| Brière (unimodal, from vector thermal biology) | 2,832 | +0.596 | +0.037 | **0.0%** |
| Log-linear (the conventional choice) | 2,832 | −0.005 | +0.005 | **48.8%** |

The log-linear coefficient is negative in essentially half of all outbreaks and positive
in the other half, averaging to nothing. In Pakistan it came out negative, and read alone
that looked like a systematic artefact; across 236 outbreaks the truth is worse — the sign
it reports is close to a coin flip, determined by whether temperature happens to be rising
or falling during that particular epidemic's growth phase. A unimodal response derived
from thermal biology cannot produce a negative value at all.

### What this establishes

The field's own systematic review of 99 dengue models (PLOS NTD 2023) reports that
validation and methodology are inadequately described. This quantifies the consequence:
across 236 outbreaks worldwide, **the published conclusion about climate-driven
transmission would change in 83% of cases had the analyst made different, equally
defensible, choices** — and the most consequential of those choices is the observation
model, which is a statistical convention rather than an epidemiological hypothesis.

---

## 2026-07-27 — A remedy, and it is nearly free

`scripts/18_stability_remedy.py`. Step 17 diagnosed the problem; a paper that stops there
is a complaint. Four candidate decision rules were evaluated on the same 5,664 fits, each
scored on two things that must be reported together — how often it gives an unstable
answer, and how often it gives any answer at all. A rule that is never wrong because it
never speaks is not a remedy.

### The instability lives entirely in comparisons that were never decisive

| Rule | Unstable | 95% CI | Gives an answer for |
|---|---|---|---|
| Sign of ΔAIC (conventional) | **82.6%** | 77–87% | 100% |
| \|ΔAIC\| > 2 | 72.0% | 67–78% | 100% |
| **\|ΔAIC\| > 4** | **5.5%** | 3.0–8.5% | **99.6%** |
| \|ΔAIC\| > 10 | 4.8% | 2.6–7.5% | 97.9% |
| \|ΔAIC\| > 20 | 3.7% | 1.4–6.2% | 91.9% |
| \|ΔAIC\| > 50 | 4.5% | 1.7–8.2% | 75.4% |

Between a margin of 2 and a margin of 4 the instability collapses from 72% to 5.5%, and
the rule still returns a verdict for 99.6% of windows. Beyond 4 there is nothing further
to gain: the curve is flat to a margin of 50, while the share of windows answered falls
away.

The interpretation is direct. **Almost all of the instability lives in comparisons where
ΔAIC is smaller than 4** — that is, where the two models were never distinguishable to
begin with. Burnham and Anderson's long-standing rule of thumb treats ΔAIC of 4–7 as
"considerably less support" and below 4 as substantial support for both models. The
guidance already exists. What this shows is the cost of not applying it: reporting the
sign of a difference smaller than the threshold at which the difference means anything,
which is what produces a conclusion that moves with the analyst.

### Stricter rules buy nothing more and cost a great deal

| Rule | Unstable | Answers |
|---|---|---|
| Negative binomial only, sign of ΔAIC | 69.9% | 100% |
| Negative binomial, \|ΔAIC\| > 10 | 2.8% | 45.8% |
| All 24 combinations must agree | 0.0% | 17.4% |
| All 24 agree and \|ΔAIC\| > 10 | 0.0% | 2.1% |

Requiring unanimity across the factorial is perfectly stable and useless: it speaks for one
outbreak in six. Fixing the observation model alone — the largest single lever — leaves
70% instability, because the remaining choices still move the verdict. The margin is the
only rule tested that removes the instability while still answering.

### The headline survives the checks a reviewer would demand

**Not a small-sample artefact.** By quartile of reported cases: 89.8%, 84.7%, 74.6%,
81.4% unstable. By series length: 84.9%, 86.3%, 83.9%, 75.0%. Larger outbreaks are
marginally more stable, and nowhere near stable.

**Not an artefact of our own window-selection rule.** This study applied a
peak-prominence threshold when choosing windows, and must hold itself to the standard it
applies to others. Recomputing the headline on progressively stricter subsets: 82.6% at no
threshold, 82.9%, and 85.5% at the strictest (n = 76). Tightening our selection makes the
result slightly stronger, not weaker.

**Interval on the headline.** Bootstrapping over windows, which are the independent units:
82.6% [77.1, 86.9].

### What the paper can now recommend

Not "do not use AIC", and not "climate does not drive dengue". Something narrower and
actionable:

> Model-selection verdicts on climate-driven transmission should not be reported from the
> sign of ΔAIC. On 236 outbreaks, requiring a margin of 4 reduces the share of outbreaks
> whose verdict depends on unremarked analysis choices from 83% to 5%, while still
> returning a verdict for 99.6% of them. Below that margin the correct report is that the
> data do not distinguish the models.

---

## 2026-07-27 — Other criteria, and the mechanism

Two additions, both prompted by objections a referee would raise.

### Is the instability specific to AIC? No — stricter criteria are worse

`scripts/19_criteria_comparison.py`. AIC, BIC and the likelihood-ratio test are all
thresholds on the same quantity, and all three are computable from what was already
stored: AIC = −2 logL + 2k gives logL, and Δk = 2 in every comparison. Nothing was
refitted.

| Criterion | Equivalent to ΔAIC < | Unstable | 95% CI |
|---|---|---|---|
| AIC, sign | 0.00 | 82.6% | 78–87% |
| LRT at p < 0.05 | −1.99 | 88.1% | 84–93% |
| BIC, sign (38-week median) | −3.28 | 91.9% | 88–95% |
| LRT at p < 0.01 | −5.21 | 94.1% | 91–97% |
| **AIC with a margin of 4** | abstain in (−4, +4) | **5.5%** | 3–9% |
| BIC with a margin of 4 | abstain, off-centre | 79.2% | 74–84% |

**Making the criterion stricter makes the instability worse.** The conventional
significance test an epidemiologist would run — a likelihood-ratio test at 5% — is
equivalent to a margin of about 2 and leaves 88% of outbreaks unstable. And a margin of 4
on the BIC scale, which sounds like the same remedy, leaves 79%, because on this problem
BIC's extra penalty shifts the abstention band off centre.

### Why: a hypothesis tested and rejected, then the real explanation

The first explanation attempted was that a threshold flips when comparisons sit near it,
so instability should track the density of comparisons at the boundary. **It does not.**
Scanning the threshold from −14 to +6 gives 95% instability where the local density is 4%
and 59% where it is 16%; the Spearman correlation is −0.43, the wrong sign. Recorded here
because it was tested, and because the correct explanation only became visible once it had
been ruled out.

What is actually happening is two facts that sound contradictory until separated:

**The choices move the evidence enormously.** The range of ΔAIC across the 24 combinations
of a single outbreak has a median of **172 AIC units** (quartiles 53 to 808, 90th
percentile 5,009). The median |ΔAIC| being read as evidence is **9.9**. The analyst's own
choices move the evidence by roughly **seventeen times** the size of the difference being
interpreted.

**Yet the disagreement is marginal.** The margin needed to eliminate all disagreement
within a window has a median of 3.50 and a 95th percentile of 4.31. In **94.5%** of
outbreaks, every combination that dissents from the majority verdict does so on evidence
weaker than ΔAIC = 4. The margin needed is uncorrelated with how far the choices moved the
evidence (ρ = −0.06, p = 0.39).

Both hold because most of the enormous variation stays on one side of zero: it changes how
strong the evidence looks without changing which model it favours. Where the direction does
flip, it flips on comparisons that were weak to begin with.

That is the mechanism, and it explains all three observations at once — why a small margin
works, why the number 4 is not special (it is where the dissent lives, not a property of
the criterion), and why a stricter threshold cannot help, because a stricter threshold is
still a threshold and the choices straddle wherever it is placed.

### Revised recommendation

The recommendation is now sharper than "use a margin of 4":

> Report a verdict only where the evidence exceeds the variation your own analysis choices
> induce in it. Measuring that variation requires refitting under alternative defensible
> choices, which costs a factorial and is what this study did. Absent that, abstaining
> below |ΔAIC| = 4 reproduces the same behaviour on 94.5% of outbreaks. Do not substitute
> a stricter criterion: BIC and the likelihood-ratio test are both less stable than AIC
> here, not more.

---

## 2026-07-27 — False positives: the verdict is not just unstable, it is wrong

`scripts/21_false_positive_rate.py`. Everything until now measured disagreement, which
shows something is wrong without showing *which* answer is wrong — on real data there is
no truth to compare against. So counts were simulated from models whose answer is known,
using each window's own fitted parameters, its own climate series and its own estimated
dispersion, and the full five-factor factorial was run on them.

60 windows, two truths, 48 combinations, 11,520 fits.

| Rule | False positive rate | Power |
|---|---|---|
| Poisson, sign of ΔAIC | **80.8%** | 89.0% |
| Poisson, margin of 4 | 69.8% | 80.2% |
| **Negative binomial, sign of ΔAIC** | **14.2%** | 70.7% |
| Negative binomial, margin of 4 | **3.9%** | 52.6% |

*False positive rate*: data generated from the constant-transmission model, with the real
climate covariates present but exerting no influence. Any endorsement of climate forcing
is wrong. *Power*: data generated with a genuine effect of realistic size.

### Three findings, in order of importance

**Assuming Poisson gives an 81% false-positive rate.** On data containing no climate
signal whatsoever, a Poisson observation model endorses climate forcing four times in
five. Its apparent power of 89% is not a virtue: a rule that says yes to 81% of null data
and 89% of real data is barely distinguishing them at all.

**The observation model is the fix, not the margin.** Switching to a negative binomial
takes the false-positive rate from 80.8% to 14.2% while retaining 71% power — by far the
largest single improvement available, and it costs nothing but a defensible choice.

**A margin cannot rescue Poisson.** Applying the same margin of 4 that works so well
elsewhere leaves Poisson at 69.8% false positives, against 3.9% for the negative binomial.
The two remedies are not interchangeable: the margin removes marginal disagreement, but
Poisson's errors are not marginal — they are large and confident.

### How this connects to the empirical study

The empirical finding was that Poisson endorses climate forcing in 92.1% of real fits and
the negative binomial in 56.5%. The simulation explains the gap: on null data Poisson
endorses 80.8% of the time. Most of the difference between the two observation models on
real data is not detection, it is false positives.

### The recommendation, now quantified

1. **Use a negative binomial observation model.** False positives fall from 81% to 14%
   with power at 71%. This single change matters more than everything else tested.
2. **Add a margin of |ΔAIC| > 4 where a false positive is costlier than a miss.** With the
   negative binomial this gives 3.9% false positives at 53% power.
3. **Do not substitute the margin for the observation model.** On Poisson counts it leaves
   70% false positives, and neither BIC nor a stricter likelihood-ratio test helps either
   (step 19).

### Operating characteristics, and the recommendation in its final form

`scripts/22_operating_characteristics.py`. Reporting either error rate alone is
misleading, so both are given together across a range of margins.

| Observation model | Rule | False positives | Power | Separation |
|---|---|---|---|---|
| Poisson | sign of ΔAIC | 80.8% | 89.0% | **8 pts** |
| Poisson | margin 4 | 69.8% | 80.2% | 10 pts |
| Poisson | margin 16 (its best) | 51.1% | 68.8% | 18 pts |
| **Negative binomial** | **sign of ΔAIC** | 14.2% | 70.7% | **56 pts** |
| Negative binomial | margin 4 | 3.9% | 52.6% | 49 pts |
| Negative binomial | margin 12 | 0.9% | 34.4% | 34 pts |

Separation is power minus false-positive rate: how much the rule actually distinguishes
the two truths, as opposed to how often it says yes.

**Poisson does not discriminate.** Its best achievable separation, over every margin
tested, is 18 points, and at the conventional rule it is **8**. An 89% power figure means
nothing next to an 81% false-positive rate: the procedure endorses climate forcing
whether or not any exists.

**The observation model is worth 48 points of discrimination.** Switching to a negative
binomial moves separation from 8 to 56 at the identical decision rule. Nothing else tested
in this project comes close to that.

**The margin is a specificity dial, not an improvement.** On the negative binomial it
trades power for false positives along the same curve — 14.2%/70.7% at margin 0, 3.9%/52.6%
at margin 4 — rather than shifting the curve outward. It is the right dial to turn when a
false positive costs more than a miss, and it is not a substitute for the observation
model.

### The bug I had already fixed, reintroduced

The two truths in step 21 were labelled `null` and `climate`. Pandas reads `null` back as
`NaN`, so the entire no-effect arm — 2,832 fits — vanished on re-import, exactly as the
`constant` model's rows had earlier in this project.

`tests/test_outputs.py` was written for precisely this and would have caught it; I ran the
analysis before running the tests. The label is now `no_effect`, the script carries a
comment explaining why, and the data were relabelled in place rather than recomputed, the
fits themselves being unaffected.

Worth recording plainly: knowing about a class of bug and having a test for it does not
help if the test is run after the result has already been read.

---

## 2026-07-27 — Model structure as a fifth factor, and two errors it exposed

The strongest objection left to the robustness result was that it might be an artefact of
one model. Every fit used a host–vector formulation; a referee is entitled to ask whether
the instability is a property of dengue analysis or a property of *this* transmission
model. So the factorial gained a fifth factor: host–vector, or a directly transmitted
human-only SEIR, which is what much of the applied literature actually fits to case counts.
The two are matched on R₀ so the comparison is between structures rather than between
transmission intensities.

**11,326 fits, 226 minutes, 236 windows attempted, 234 completing all 48 combinations.**
One Brazilian window was dropped when every start diverged on the anchor fit.

### The answer to the objection

| Factor varied | Verdicts flipped | P(climate wins) before → after |
|---|---|---|
| **Negative binomial → Poisson** | **38.1%** | 0.542 → 0.917 |
| 75% → 100% of the series | 22.4% | 0.682 → 0.778 |
| Rainfall lag 3 → 7 weeks | 15.6% | 0.735 → 0.726 |
| Brière → log-linear temperature | 11.5% | 0.682 → 0.778 |
| Rainfall lag 3 → 5 weeks | 10.4% | 0.735 → 0.728 |
| **Host–vector → human-only SEIR** | **8.3%** | 0.741 → 0.718 |

Model structure is the *weakest* of the five. The one factor that is genuinely about
dengue epidemiology — whether the mosquito is represented at all — moves the verdict less
than any convention surrounding it. That is a better result than a null: it does not merely
fail to undermine the finding, it sharpens it.

Instability rises from 82.6% to **88.0%** (95% CI 82.9–92.3), as it must when a design
gains a factor. 28 windows (12.0%) always favour climate forcing, none never do, and among
unstable windows the median split is 36 of 48.

Out-of-sample and in-sample instability are now **83.7% and 83.7%** — equal to the decimal
on 233 windows. Endorsement still halves out-of-sample, 68.2% → 35.4%. The correction
recorded on 2026-07-26 holds at larger sample size.

The Brière temperature exponent is negative in **0.0%** of 5,616 fits; the log-linear
coefficient in **48.0%**. Unchanged by doubling the sample.

### An error that did not fail: a bound tied to the old design size

Step 18 decided whether a window was unstable with `wins.between(1, 23)` — correct for a
24-combination factorial, wrong for 48. It did not crash. It reported **11.5%** instability
where the true figure was **88.0%**, because a window in which 24 to 47 of the 48
combinations favoured climate was counted as unanimous. The headline printed in section 1
of the same script was right, so the two halves of one output disagreed by a factor of
eight and I nearly wrote the wrong one into the paper.

The fix is not the corrected bound. The definition now lives in `src/dengue_pk/robustness.py`,
reads the factorial size from the data, and is covered by `tests/test_robustness.py`, which
checks every property at both 24 and 48 combinations — so a rule that is right at one size
and wrong at the other fails. One test asserts the old bound's exact failure, to keep it
from being reintroduced quietly.

`scripts/23_paper_numbers.py` now prints every number the paper quotes from the result
tables in one place. The paper had drifted from the pipeline twice; catching it by reading
is not a method.

### A selection bias in the simulation study

The false-positive study generated data under two truths from each window's own fitted
parameters. The alternative was built by switching the climate coefficients on — and that
multiplies β(t) by a factor whose mean is not 1, so the alternative had *higher mean
transmission* than the null, not merely time-varying transmission.

Two consequences. The contrast being measured was confounded: any detected difference could
be a difference in R₀ rather than in climate forcing. And it broke the study quietly —
for **21 of 59 windows the raised R₀ made the generating integration diverge**, so the
climate arm covered 38 windows against the null arm's 59. Power was therefore estimated
only on windows where a strong effect happened not to break the model, which is selection
on something related to the outcome.

β₀ is now rescaled by the mean of the climate multiplier, so both truths have the same
average transmission and differ only in whether it varies. On a four-window check both arms
now complete all 48 combinations. Re-running at **100 windows** rather than 60 for tighter
intervals.

### Is the instability just optimiser noise?

`scripts/24_optimiser_check.py`. The sharpest objection left: the factorial fitted each
model with three multi-starts warm-started from an anchor, chosen for runtime. If three
starts sometimes miss the optimum — and miss it more often for the climate model, which
carries two extra parameters and a rougher surface — then ΔAIC would move for reasons that
have nothing to do with the observation model or the rainfall lag.

Seven outbreaks refitted under all 48 combinations twice: **as run** (3 restarts, warm) and
**thorough** (10 restarts, cold, no warm start). 329 paired comparisons.

| Question | Answer |
|---|---|
| Does thorough find better optima? | Yes — 6.1% of climate fits, 3.3% of constant fits; median gain 0.19 AIC, **max 1,656** |
| Does the verdict change? | **0.6%** of comparisons (weakest genuine factor: 8.3%; strongest: 38.1%) |
| Does the headline change? | No — all seven outbreaks unstable under both settings |

This is the strongest form the answer could take. It is *not* "the optimiser was already
perfect" — it demonstrably was not, and on one fit three warm starts missed by 1,656 AIC
units. It is "the optimiser has room to improve and improving it does not rescue the
conclusion", which is a claim about the evidence rather than about the software.

Reported honestly as a bound rather than an estimate: seven windows, not twenty. The run
was stopped early because it was competing with step 21 for the same machine and the answer
was already unambiguous at that sample size. The script gained an `analyse` mode so the
stored partial table could be reported without refitting — worth having in general, since
every long script here writes after each window precisely so an interrupted run is not a
wasted one.

### Positioning: this is a multiverse analysis, and it should say so

A reviewer would ask how this relates to multiverse analysis (Steegen et al. 2016) and
specification curve analysis (Simonsohn et al. 2020), and the manuscript did not answer.
It does now, in the introduction. The honest positioning is that the *method* is theirs and
the *domain* is not: those literatures grew up around regression, where the forking paths
are covariate sets and exclusion rules. Nothing here has an analogue in a covariate list.

Two differences are real rather than cosmetic. Mechanistic fits can fail to converge, so
the multiverse is not guaranteed complete and has to declare its missing cells — two of 236
windows here. And the comparison is between nested mechanistic hypotheses judged by an
information criterion, which makes the decision rule a lever that can be tested rather than
a fixed part of the setup. That second point is what turns the paper from a complaint into
a recommendation, and it is not available in the regression setting the method came from.

Searched for prior multiverse analyses of compartmental epidemic models and found none.
Also found Khamthong & Phramrung (PLOS NTD 2026), which reaches a compatible conclusion by
a different route — a station with the stronger in-sample climate–dengue association
forecast no better out of sample. Cited: two studies with no design in common finding that
in-sample climate association overstates what the data support.

### The figure the paper was missing

`scripts/25_specification_curve.py`. Having positioned the study as a multiverse analysis,
it was missing that literature's canonical figure — the specification curve. Added, in two
panels because the question has two forms.

**The single-outbreak panel gave the study its clearest sentence.** The most evenly split
window is Argentina, December 2015: exactly 24 of 48 specifications favour climate forcing
and 24 favour the constant model, and ΔAIC runs from **−2050 to +4**. One dataset, one
epidemic. An analyst following one set of defensible conventions reports overwhelming
evidence; another following a different set reports none.

**The across-outbreak panel gave it a new number.** Give each of the 48 specifications the
whole study and ask how often it endorses climate forcing:

| | Specification | Endorses climate in |
|---|---|---|
| Most credulous | Poisson, log-linear, lag 3, host–vector | **96.6%** of outbreaks |
| Least credulous | negative binomial, Brière, lag 7, SEIR | **35.0%** of outbreaks |

**A spread of 62 percentage points.** This is the same finding as the 88% instability
expressed as a rate rather than as a count of flips, and it is the more useful form for
reading the published literature: a paper reporting that climate forcing improved its fit
has told you where in that 62-point range it sat only if it stated all five choices.

Two drawing decisions worth recording, because the first version of the figure was
misleading. On a linear ΔAIC axis the single −2050 bar compresses all 24 constant-favouring
specifications into the axis line, and the panel appears to show unanimity — the exact
opposite of what it exists to show. The scale is now symmetric-log with the linear region
set to ±4, which is also the margin the paper recommends, so the scale change and the
"not distinguishable" band are the same boundary. And every in-axes annotation collided
with a bar at this aspect ratio, so the explanation moved into a legend.

The observation-model row splits the left panel exactly: every Poisson specification on the
climate side, every negative binomial on the other. That is the paper's main result visible
without reading a number.

### Answering "Poisson isn't defensible, so your multiverse is inflated"

A reviewer can reasonably say that a Poisson likelihood on overdispersed counts is simply
wrong, that including it admits a cell that should never have been in the multiverse, and
that the honest headline is the instability among *defensible* analyses only.

The answer was already in the results and was not being used. Restricting the whole study
to negative-binomial fits — discarding half the factorial — still leaves **76.5%** of
outbreaks changing verdict (95% CI 71–82). The observation model is the largest lever, not
the only one, and the finding survives conceding the objection in full.

The second half of the answer is that the field does use Poisson: reviews report that most
dengue studies use Poisson regression or comparably restrictive count models and that
overdispersion-tolerant formulations remain uncommon. A multiverse that excluded it would
not describe the literature it is about.

Both points are now stated in Section 4.2 rather than left implicit in a table. Worth
recording as a general lesson: several of this paper's strongest defensive answers were
already computed and simply not written down where the objection would be raised.

### A concrete recommendation, not just a diagnosis

The discussion now ends with four numbered recommendations ordered by what they are worth
on the evidence: overdispersed observation model (48 points of separation), abstain within
±4 (costs 0.4% of verdicts, removes 92% of the instability), state all five choices (62
points of endorsement rate lie between the extremes), and refit under the alternatives
where affordable (96 fits per outbreak, minutes).

None of it is new statistical advice, and the paper says so. The contribution is the
measured cost of not following it.

---

## 2026-07-27 — The simulation, rerun without its selection bias

`scripts/21_false_positive_rate.py`, rebuilt and rerun at **100 windows, two truths, 9,586
fits, 317 minutes**. The previous run had a bias serious enough to invalidate its power
estimates: the alternative truth was built by switching the climate coefficients on, which
multiplies β(t) by a factor whose mean is not 1, so it differed from the null in *overall
transmission* as well as in whether transmission varies. For 21 of 59 windows the raised R₀
made the generating integration diverge, leaving the climate arm on 38 windows against the
null arm's 59 — power measured only where a strong effect happened not to break the model.

β₀ is now rescaled by the mean of the climate multiplier. Both truths have the same average
transmission and differ only in whether it varies with climate, which is the contrast the
study is actually about. Both arms now cover all 100 windows (4,800 and 4,786 fits).

| Observation model | Rule | False positives | Power | Separation |
|---|---|---|---|---|
| Poisson | sign of ΔAIC | 83.0% | 94.7% | **12 pts** |
| Poisson | margin 4 | 74.6% | 92.0% | 17 pts |
| Poisson | margin 24 (its best) | 48.6% | 85.2% | 37 pts |
| **Negative binomial** | **sign of ΔAIC** | 13.8% | 72.7% | **59 pts** |
| Negative binomial | margin 4 | 3.9% | 61.3% | 57 pts |
| Negative binomial | margin 12 | 1.4% | 45.6% | 44 pts |

The conclusions from the biased run survive and sharpen. Poisson endorses climate forcing on
**five datasets in six containing no signal**; its best separation over every margin tested
is 37 points and reaching that requires abstaining below ΔAIC = 24. The observation model is
worth **47 points** of separation at the identical rule.

### The result I did not expect

The clean design made a new question answerable: **does a factorial manufacture instability
by itself?** Run enough analyses and something disagrees — so a 88% instability rate might
say more about running 48 fits than about the evidence. Here the truth is known and the
factorial is identical in both arms, so the comparison is direct.

* Data with **no** climate effect: **96.0%** of windows unstable.
* Data with a **real** effect of realistic size: **71.4%**.

The 25-point gap, in the expected direction, answers the objection: the design is not
generating disagreement out of nothing, and the truth genuinely moves the verdict.

The 71.4% is the uncomfortable half, and it changes what this paper is claiming. I had been
reading the instability as a symptom of weak evidence — comparisons near zero flip, strong
ones do not. That reading is incomplete. A genuine climate effect, at a size drawn from
these epidemics' own fitted parameters, still has its verdict flipped by the analyst's
choices in seven outbreaks in ten. The methodology does not reliably deliver a stable answer
even when there is a real answer to deliver.

Recorded plainly because it weakens a sentence I had already written and liked. The
instability is not only weak evidence looking weak.

### Does the recommendation survive a real effect?

The 71.4% sat awkwardly against Section 4.7's account, which says the dissent lives in
comparisons too close to call. If that were the whole story, a margin should stabilise the
real-effect arm as well as the null arm. It does not, and the gap is informative.

| Rule | No climate effect | Real climate effect |
|---|---|---|
| sign of ΔAIC | 96.0% unstable | 71.4% unstable |
| abstain within ±4 | **7.4%** unstable | **28.6%** unstable |
| abstain within ±10 | 1.1% unstable | 25.8% unstable |

The margin works in both cases and works less completely where an effect exists, which is
the right direction: a real effect produces large |ΔAIC|, so the dissent is no longer
confined to the abstention band. A margin removes disagreement about comparisons too close
to call. It cannot remove disagreement about a difference the analyses genuinely disagree
about. This is a limitation of the recommendation and is now stated as one.

### The inference I nearly didn't write down

On **real** outbreaks the margin-4 rule leaves **6.9%** unstable. That is close to the
**7.4%** the same rule gives on simulated data with **no** climate effect, and far from the
**28.6%** it gives on simulated data with a real one.

Taken at face value, these 234 real outbreaks behave under this diagnostic more like the
no-effect arm than like the real-effect arm.

It is in the paper, flagged explicitly as suggestion rather than result. The simulation
shares each window's fitted parameters, climate series and dispersion but not its unmodelled
structure — serotype dynamics, reporting-effort changes, the spatial heterogeneity the
worked case demonstrates — and any of those could produce the same resemblance without any
statement about climate being true.

The reason it is written down at all: leaving out an observation because it is inconvenient
to characterise is exactly the class of analytical choice this paper exists to measure.

---

## 2026-07-27 (later) — Answering the two objections the paper had conceded

The five-factor paper named two weaknesses in its own limitations section and left them
there. Both are the first thing a reviewer would reach for, precisely because the paper
handed them over. They are now factors rather than caveats.

### "Your model is too simple"

Every fit held mosquito lifespan constant across seasons whose temperature swings fifteen
degrees. Mechanistic dengue modelling does not do that, so the honest reading of the old
result was "a model simpler than current practice is unstable" — which is not the claim
the paper wanted to make.

A third structure was added: the same host–vector model with adult mosquito mortality
following a quadratic lifespan response, thermal limits fixed from *Aedes aegypti* trait
data (11.7–37.2 °C), normalised so that at the thermal optimum mortality equals the
constant value the simpler structure uses. **It estimates no additional parameter**, so
the structure factor stays a comparison between mechanisms rather than between one model
given more freedom than another.

Implementation detail worth recording: the integrator now carries a mortality grid
evaluated on the same half-step grid as β, rather than branching inside the step. The
constant case fills that grid with a scalar, so there is one code path, and
`test_constant_mu_v_t_reproduces_the_scalar_path` asserts the refactor left every existing
result untouched to 1e-14. Without that test the change would have silently moved every
number in the study for a reason unrelated to the new structure.

**A subtlety I nearly shipped without stating.** Under this structure temperature enters
*both* compared models, because it drives mortality in the structure itself. Its "constant"
model is therefore not climate-free — it is "climate affects mosquito survival but not the
transmission coefficient". That makes it the hardest of the three tests rather than an
unfair one, since its null already contains thermal biology, but a reader who was not told
would reasonably assume a like-for-like comparison. It is now stated in the script, in
`docs/MODEL.md`, and will be stated in the paper.

On a four-window check the new structure completes all 144 combinations like the others,
recovers similar R₀ (1.35 against 1.42 and 1.48), and endorses climate forcing *more often*
(0.91 against 0.81 and 0.80) — consistent with the extra channel through which climate can
matter.

### "You fixed the parameters that matter"

Mosquito lifespan, both incubation periods and the vector-to-host ratio were fixed from the
literature and never varied, while the study criticised others for exactly this class of
unexamined choice. Published ranges are wide — viraemia 4–7 days, incubation 4–7 days,
adult lifespan 8–15 days, vector-to-host ratio from below one to above ten — so choosing
within them is a degree of freedom like the rainfall lag.

A second parameter set at the other end of the same ranges is now a factor
(`config: model.fixed_alt`). It maps β to R₀ differently (13β against 10β), which is
exactly why it belongs in the factorial rather than in a sentence.

**A config edit that broke nine tests.** The new block was first placed between `fixed:`
and the `rho_fixed` key that belonged to it, which silently moved `rho_fixed` into the new
block. Nine tests failed with `KeyError: 'rho_fixed'`. YAML indentation has no closing
brace to get wrong, so a block inserted in the middle of another one looks correct.
The block now sits after the keys it must not capture, with a comment saying why.

### Six factors, and making it run in two hours instead of twenty-two

    observation x2 · temp_form x2 · rain_lag x3 · train_frac x2 · structure x3 · params x2
    = 144 combinations, two models each = 288 fits per window

Single-threaded that is about 22 hours for the full inventory. `scripts/26_*` distributes
windows across twelve processes: the parent does all the I/O — reading the 2.2-million-row
case table once, building each window's Datasets — and workers receive prepared Datasets
and return rows. Workers hold no copy of the case table, which is what makes twelve of them
possible on a laptop. Expected wall clock is under two hours.

Table selection also moved into `dengue_pk.robustness.FACTORIAL_TABLES`, so growing the
design means adding one entry rather than editing six scripts. That is a direct response to
the four-to-five-factor change, which left a hard-coded bound behind in step 18 and
reported 11.5% where the truth was 88.0%.

### A results section nearly deleted by a column that was not carried forward

The six-factor script was written fresh rather than edited from the five-factor one, and in
the process it stopped recording `heldout_delta`. Nothing failed. The run was 60 windows in
and producing a perfectly good table — one that simply could not answer the out-of-sample
question, which is an entire results section and the one showing that validation halves
endorsement without restoring stability.

Caught by running step 17 against the partial table and noticing section 3 had nothing to
work with. Fixed, and the run restarted from scratch at a cost of about an hour of compute,
because resuming would have left the first sixty windows without the column.

The general lesson, which is now the third instance of it in this log: **rewriting a script
loses whatever the old one did that the new author forgot it did.** Each of the three was
silent — a wrong bound, a mislabelled truth, a missing column — and none produced an error.
The defence that has actually worked is running the downstream step early on partial output
rather than waiting for the run to finish.

### A statistic that does not grow with the size of the factorial

Enlarging the design from 48 combinations to 144 raises "the verdict changes in X% of
outbreaks" mechanically: asking whether *any* two cells disagree gets easier the more cells
there are. A reviewer would notice that the headline went up when the only thing that
changed was the enumeration, and would be right to.

`pairwise_disagreement` reports instead the probability that two analyses of the same
outbreak, drawn at random, reach opposite verdicts — `2p(1-p)`, a function of the proportion
alone. Zero for a unanimous outbreak, one half for an evenly split one, and identical for a
quarter-split whether the design has 24 cells or 144. `tests/test_robustness.py` asserts
that invariance exactly.

Drawing *without* replacement is the more literal reading of "two different analyses" and
gives `2w(n-w)/[n(n-1)]`, but that carries a factor `n/(n-1)` and does creep with design
size — 0.391 over 24 cells against 0.378 over 144 at a quarter split. Small, and entirely
an artefact of the enumeration, which is the artefact the statistic exists to remove, so
the independent form is used and the correction deliberately dropped.

On the first 52 windows of the six-factor run it reads **25.6%**: two defensible analyses of
the same outbreak, chosen at random, reach opposite conclusions about a quarter of the time,
against a maximum possible 50%.

### Holding this study to its own standard

A study that keeps adding factors until its headline is large is doing exactly what it
criticises. The defence has to be evidence, not assertion, so step 23 now prints every
design that was run side by side:

| Design | Combinations | Verdict changes | P(two analyses disagree) |
|---|---|---|---|
| 4-factor | 24 | 82.6% | 29.9% |
| 5-factor | 48 | 88.0% | 31.4% |
| 6-factor | 144 | *pending full run* | *pending* |

The middle column climbs, and it should: asking whether *any* two cells disagree gets
easier the more cells there are. The right column, which depends on the proportion rather
than the count, sits near 30% throughout. Reporting both is the only honest way to present
a headline that is partly an artefact of the enumeration.

Two further points belong with that table. Each factor added closed a gap the *previous*
version had named in its own limitations section — model structure and the fixed
epidemiological parameters were both written down as weaknesses before they were measured,
not sought once the result was known. And the direction was not guaranteed: the parameter
set turned out to flip almost nothing, which is a null result for a factor I expected to
matter and is reported as one.

### The result the six-factor design was really for

At 117 of 237 windows the paired flip rates already separate cleanly, and they separate
along a line the paper had not drawn before:

| Factor varied | Verdicts flipped |
|---|---|
| Negative binomial → Poisson | **30.6%** |
| 75% → 100% of the series | 19.4% |
| Rainfall lag 3 → 7 weeks | 14.4% |
| Brière → log-linear temperature | 11.8% |
| Host–vector → thermal-mortality structure | 10.5% |
| Rainfall lag 3 → 5 weeks | 10.0% |
| Host–vector → human-only SEIR | 9.6% |
| **Fixed epidemiological parameters** | **2.8%** |

The factors sort into two groups. The **statistical conventions** — which likelihood, how
much of the series to fit — flip 19–31% of verdicts. The **epidemiological content** —
which transmission mechanism, what the biological parameters are — flips 3–10%. The
question is answered by the statistics, not by the biology.

That is a sharper claim than the five-factor version could make, and it is sharper in a
direction I did not choose: **the fixed parameter set, the factor I was most concerned
about having excluded, turns out to flip 2.8% of verdicts.** I added it expecting it to
matter and it does not. Recorded as a null result rather than dropped, because a factor
included only when it strengthens the headline is the practice this paper exists to
measure.

The pairwise disagreement statistic reads 27.6% at this point, against 29.9% and 31.4% for
the four- and five-factor designs — flat, as it should be, while the "any two disagree"
share climbed from 82.6% to 90.6%.

### An honest caveat on the location test

Building the climate-location study surfaced something worth stating before its result
exists. The alternative point is one degree of latitude away, applied mechanically, with no
regard for terrain. In flat units the two series differ by a few tenths of a degree, which
is the intended comparison. For Costa Rica they differ by **4.3 C** in mean temperature:
the offset has crossed high ground, and that is not a choice between two defensible
stations for the same population but between two climates.

So whatever flip rate it returns is an **upper bound** on what a sensible alternative choice
would do. The alternative — hand-picking a second city per reporting unit — would replace a
rule anyone can reproduce with a gazetteer of our own choices, which is not a trade this
study is in a position to make quietly.

### The conclusion I would have got wrong

Running step 18 against the six-factor table appeared to destroy the paper's remedy. On
48 combinations, abstaining within ±4 cut instability from 88% to **6.9%**. On 144 it cut
91% to **39%**. Read at face value: the remedy does not survive a realistic multiverse, the
paper's one constructive contribution is gone, and what is left is a complaint.

I nearly wrote that. The first check was whether the third structure was responsible, since
its null differs from the others' and it could dissent *confidently* rather than marginally.
Excluding it moved 39.4% to 32.0% — not the explanation.

The actual explanation is the one this log had already written down two entries earlier and
I failed to apply to the remedy: **"some pair still disagrees" gets easier to satisfy the
more cells there are, at every margin, not just at zero.** The instability share is not
comparable across designs, and neither is any rule evaluated on it.

On the design-invariant scale the remedy is intact:

| Margin | P(two speaking analyses disagree) | Windows still answering |
|---|---|---|
| 0 | 29.4% | 100% |
| 2 | 24.6% | 100% |
| **4** | **3.9%** | 100% |
| 10 | 4.2% | 98.9% |
| 20 | 4.6% | 97.2% |

The collapse between 2 and 4 is exactly where the five-factor study found it, at the same
size. Nothing about the remedy changed; the yardstick did.

Two things follow. The paper must report the remedy on the invariant scale, with the
share-of-windows figure alongside and labelled as design-dependent. And the invariant
statistic has now earned its place twice: once against the objection that the headline grows
with the design, and once against my own reading of my own output. It was introduced to
answer a reviewer and its first real use was to stop me publishing a wrong conclusion.

`pairwise_disagreement_at_margin` drops windows where nothing speaks rather than scoring
them zero. Counting silence as agreement is how a remedy flatters itself: push the margin
high enough and every window falls inside the band, which would read as perfect stability.
A test asserts the silent window is absent from the result rather than present as a zero.

---

## 2026-07-27 (final) — The six-factor result

**33,152 fits, 126 minutes on twelve processes.** 237 windows in the inventory; four
Brazilian windows could not be fitted at all (every multi-start diverged on the anchor) and
twelve completed only part of the design, leaving **221 outbreaks in 33 countries**.

| Measure | Value |
|---|---|
| Outbreaks where some pair of the 144 disagrees | 91.9% [87.8–95.5] |
| **P(two analyses drawn at random disagree)** | **29.8%** (median 34.6%) |
| Always favours climate forcing | 18 (8.1%) |
| Never favours climate forcing | 0 |
| Median split among unstable outbreaks | 107 of 144 |
| Per-country range | 85%–100% |

### The factors sort into three groups, and the ordering is the finding

| Factor varied | Verdicts flipped | P(climate) before → after |
|---|---|---|
| *Statistical conventions* | | |
| Negative binomial → Poisson | **34.3%** | 0.570 → 0.900 |
| 75% → 100% of the series | 20.7% | 0.702 → 0.769 |
| *Covariate handling* | | |
| Rainfall lag 3 → 7 weeks | 15.8% | 0.734 → 0.738 |
| Brière → log-linear temperature | 11.5% | 0.696 → 0.774 |
| Rainfall lag 3 → 5 weeks | 11.0% | 0.734 → 0.733 |
| *Epidemiological content* | | |
| Host–vector → human-only SEIR | 9.9% | 0.738 → 0.727 |
| Host–vector → thermal mortality | 9.4% | 0.738 → 0.741 |
| **Central → alternative fixed parameters** | **3.1%** | 0.735 → 0.736 |

The conventions of statistical analysis flip 21–34% of verdicts; the epidemiological
content flips 3–10%. **The question is answered by the statistics, not by the biology.**
That is a sharper claim than any earlier version of this study could make, and it is
sharper because the two factors added to answer the "your model is a toy" objection both
landed near the bottom.

### Everything else, on the full design

* **Specification curve.** Most credulous analysis endorses climate forcing in 95.9% of
  outbreaks, least credulous in 39.8% — a 56-point spread. The single-outbreak panel draws
  Jamaica, July 2023: exactly 72 of 144 for and 72 against, ΔAIC from −337 to +4.
* **Out-of-sample.** Endorsement halves, 70.3% → 36.3%. Stability does not improve: 84.7%
  in-sample against **90.6%** out-of-sample. Smaller designs found these equal; the larger
  one puts out-of-sample *above* in-sample.
* **Temperature sign.** Brière negative in 0.0% of 15,912 fits; log-linear in 48.5%, mean
  zero to three decimals.
* **Criteria.** AIC by sign 91.9%, LRT p<0.05 95.0%, BIC by sign 96.8%, LRT p<0.01 97.3%.
  Strictness still makes it worse, monotonically.
* **Mechanism.** Within-window ΔAIC range median 804 units against a median |ΔAIC| of 14.4
  read as evidence — a factor of **56**.
* **Checks.** By quartile of reported cases 94.6 / 90.9 / 89.1 / 92.7%; by series length
  89.5 / 93.2 / 96.1 / 88.9%. Tightening our own window-selection rule moves the headline
  from 91.9% to 94.1%.

### One claim the larger design overturned

On 48 combinations the margin needed to remove disagreement was uncorrelated with how far
the choices moved the evidence, and the paper said dissent was purely a property of
comparisons too close to call. On 144 that correlation is **ρ = +0.61** (p = 4e-24): in
39.8% of outbreaks at least one dissenter is *confident*, and no band around zero reaches
it. What survives is that such dissenters are few — a median of 3 of the 100 analyses still
speaking after abstaining at ±4.

The script now derives that sentence from the measured correlation instead of asserting it,
because the old wording would have gone on printing after it stopped being true.

The paper and README are rewritten throughout; `tests/test_paper_consistency.py` passes,
and it now reads the current factorial through `latest_factorial` rather than its own list —
it had been silently checking the manuscript against the five-factor table.

### The climate grid cell, measured at last

`scripts/27_*` and `28_*`. The limitation the study had named three times and never
quantified: one representative point supplies the climate for a whole country or province.
Forty outbreaks refitted under all 144 combinations at a second point one degree of latitude
away — about 110 km, well inside the area any of these case series aggregates over.

**5,664 paired comparisons. The verdict changes in 12.0%**, with the endorsement
probability essentially unmoved (0.739 → 0.737). Placed against the factorial:

| Factor | Flipped |
|---|---|
| Negative binomial → Poisson | 34.3% |
| 75% → 100% of the series | 20.7% |
| Rainfall lag 3 → 7 weeks | 15.8% |
| **Climate grid cell, 110 km** | **12.0%** |
| Brière → log-linear temperature | 11.5% |
| Rainfall lag 3 → 5 weeks | 11.0% |
| Host–vector → SEIR | 9.9% |
| Host–vector → thermal mortality | 9.4% |
| Fixed epidemiological parameters | 3.1% |

Where the weather series comes from matters more than which transmission model is fitted
and more than what the mosquito parameters are, and less than any of the statistical
conventions. That is a defensible place in the ranking and not a surprising one, but it had
to be measured to be said.

Reported as an **upper bound**. The offset is mechanical and respects no terrain: for most
units the two series differ by a few tenths of a degree, but for Costa Rica by 4.3 °C,
because the shift crosses high ground. That is a choice between two climates rather than
between two defensible stations for one population. Hand-picking a second city per unit
would have replaced a reproducible rule with a gazetteer of our own choices, which is not a
substitution a paper about undisclosed analytical choices should make quietly.

The script's "for comparison" line now reads the other factors from the stored sensitivity
table instead of carrying them as literals — the first version printed the five-factor
figures next to a six-factor result.

### The simulation, on the six-factor design

`scripts/29_*`, superseding step 21. The paper claims the simulation applies "the full
factorial", and once the factorial reached 144 combinations that had stopped being true.
**80 windows, two truths, 23,000 fits, 69 minutes on twelve processes.**

| Observation model | Rule | False positives | Power | Separation |
|---|---|---|---|---|
| Poisson | sign of ΔAIC | 85.7% | 95.3% | **10 pts** |
| Poisson | margin 24 (its best) | 58.8% | 86.6% | 28 pts |
| **Negative binomial** | **sign of ΔAIC** | 18.7% | 73.4% | **55 pts** |
| Negative binomial | margin 4 | 6.2% | 61.8% | 56 pts |

Poisson endorses climate forcing on **six datasets in seven that contain no signal**. The
observation model is worth **45 points** of separation at the identical rule. Both
conclusions survive the larger design; both are slightly stronger than before.

### The artefact question, answered properly this time

On the design-invariant scale, two analyses of the same simulated outbreak disagree:

| | Sign of ΔAIC | Abstain within ±4 |
|---|---|---|
| No climate effect | **44.5%** | 1.8% |
| A real climate effect | **18.6%** | 4.3% |
| *Real outbreaks, for reference* | *29.8%* | *4.0%* |

44.5% is close to the 50% ceiling: with no signal present the procedure is very nearly
coin-flipping. The truth more than halves that. So the factorial does not manufacture
disagreement — and the 18.6% remains uncomfortable, because a real effect of ordinary size
still leaves two defensible analyses disagreeing almost one time in five.

### An inference withdrawn

The five-factor version reported that real outbreaks behaved like the no-effect arm, on the
grounds that the margin-4 rule left 6.9% of them unstable against 7.4% for no-effect
simulations and 28.6% for real-effect ones. It was flagged as suggestion rather than
result, which turns out to have been the right instinct for the wrong reason.

On the invariant scale the ordering reverses. At the sign rule the real data (29.8%) sit
between the two arms; under a margin of 4 they sit next to the *real-effect* arm (4.0%
against 4.3%, with no-effect at 1.8%). Two ways of making the same comparison point in
opposite directions, which means the diagnostic does not discriminate.

The claim is withdrawn from the paper and the README, and the withdrawal is stated rather
than the paragraph quietly deleted. It was an artefact of measuring on a scale that grows
with the design — the third time in this project that scale has misled me, and the second
time in one day.

`window_verdicts` also had to stop requiring `weeks` and `cases`: the simulation tables
carry neither, so the invariant statistic could not be computed on exactly the data it was
most needed for. Those columns are now carried through when present and ignored when not.

---

## 2026-07-28 — Building rather than measuring the tape

Amna's objection, and it was right: everything so far says *the measuring tape is wrong*.
Nobody gets a paper into a good journal for that alone. A reviewer reads it, agrees, and
rejects. She asked for a house — something that was not known before and that this work
produces, not merely a warning about how others work.

She also gave a method: **published limitations sections do not contain the real
limitations; the methodology section does.** That turned out to be the useful instruction.

### First attempt, and it was already built

The idea: dengue transmission responds to temperature along a unimodal curve, so a season
spent near the flat optimum offers nothing to detect however strong the biology, while a
season crossing the steep flank offers a great deal. Detectability should track the gradient
traversed. A researcher could then tell in advance whether their data can answer the
question at all.

Searched before building, as she has insisted on since 2026-07-27. **Kirk et al. (PLOS
Climate, 2024)** had already found it, by meta-analysing 358 published temperature–dengue
correlations: the association is strongest where temperature varies most.

Her reply to that was the second useful instruction of the day: *a house that already exists
can be built better.* That is how research normally proceeds, and abandoning the idea was
my error, not a conclusion.

### Reading their methodology rather than their limitations

Two problems that no care inside a meta-analysis of published estimates can remove:

* **Analytical confounding.** Their observations are the numbers other people chose to
  report — and this project has just shown those numbers move by tens of percentage points
  with choices nobody states. If analysts in high-variation settings tend toward different
  conventions for any reason, a moderator effect appears that is about analysts.
* **Publication bias.** A correlation enters their sample by having been published, and
  "nothing found" is exactly the outcome their hypothesis predicts in low-variation settings.

Our design removes both: every outbreak fitted the same way, nothing entering or leaving the
sample by having been interesting.

### The house

Before asking *which* property of the setting explains the differences between studies,
bound how much *any* property could. Every outbreak here has been analysed every way, so the
variation splits exactly.

| Design | Analyses | Between outbreaks | Within outbreaks |
|---|---|---|---|
| Four factors | 24 | 21.7% | **78.3%** |
| Five factors | 48 | 20.3% | **79.7%** |
| Six factors | 144 | 23.5% | **76.5%** |

**Roughly three-quarters of the variation in whether an outbreak is judged climate-driven is
variation between analyses of the same outbreak.** Climate, density, wealth, serotype,
health system, surveillance intensity — everything that differs between epidemics competes
for the remaining quarter, and that quarter is an upper bound on what any moderator analysis
can hope to explain.

Crucially the split **does not grow with the design**: 78.3%, 79.7%, 76.5% across 24, 48 and
144 analyses. Adding analyses adds to both the numerator and the denominator. That is the
property the headline instability share lacks, and it is why this can be quoted where the
other cannot.

This does not say local context is irrelevant. It says context has been asked to carry more
than it can, and that a conflicting pair of published findings is more likely to differ
because the analysts differed than because the places did.

### And the hypothesis I started with failed

The mechanistic "thermal traverse" measure — mean absolute slope of the Brière response over
the season's temperatures — was supposed to beat plain temperature variance as a predictor.
It does not: against the estimated temperature exponent it reads rho = +0.13 (p = 0.06)
where the simple SD reads −0.18 (p = 0.006). Reported as a failure, because a study about
undisclosed analytical choices does not get to quietly drop the measure it hoped would work.

Every correlation with the temperature regime is weak (|rho| < 0.25), which is the same
finding from the other direction: the setting explains little.

### The second house: the intervals do not cover

`scripts/30_interval_coverage.py`. If the first contribution says where the disagreement
comes from, this one says what to do about it — and, unusually for this project, it can be
*validated* rather than argued, because the truth is set by us.

Counts simulated from each outbreak's own fitted parameters with the Brière exponent at 1.0
and the rainfall coefficient at 0.30; every Brière combination asked to recover them.
**80 outbreaks, 72 combinations each, 5,760 fits.**

| Interval | Temperature exponent | Rainfall coefficient |
|---|---|---|
| Conventional 95%, one analysis | 82.9% | **71.2%** |
| **Multiverse, Rubin's rules** | **100.0%** | **97.5%** |

**A conventional 95% interval for the rainfall coefficient contains the truth 71% of the
time.** A conventional interval reports sqrt(W) and behaves as though B were zero; for the
rainfall coefficient the analytical SD equals the sampling SD, so half the variance is
simply missing. Rubin's rules — within plus (1 + 1/m) times between, the multiple-imputation
decomposition, which transfers unchanged because the structure is identical — restore
coverage at 1.5 to 2.2 times the width.

The 100% for the temperature exponent is conservative rather than exact and is reported as
such: with m = 72 analyses and a long-tailed estimate distribution, B is generous.

Two details worth recording. Thirty per cent of individual analyses estimate the temperature
exponent at essentially **zero** when the truth is 1.0 — they conclude temperature does
nothing — while others overshoot to ten times the truth. And the temperature exponent's
conventional coverage looks almost respectable (82.9%) only because its intervals are so
wide as to be uninformative: median half-width 2.7 on a parameter whose true value is 1.0.
An interval can keep its promise by saying nothing.

Coverage by factor is flat (82–83% across observation model, structure and parameter set),
with one exception: fitting the whole series rather than three-quarters drops it from 86.9%
to 78.9%. More data, worse interval — because the estimate moves and the interval does not
widen to match.

---

## 2026-07-29 — Reviewing my own paper, and what the review cost

Amna asked for a strict reviewer's report and then asked which of its gaps could actually
be closed. Writing that report found more than reading the paper had.

### The check that saved the coverage result

A referee's first move against Section \ref{sec:coverage} is: *of course a mis-specified
model does not cover — that is the definition of bias, not a finding.* The coverage study
included fits under the wrong structure and the wrong parameter set, so the objection had
to be answered before the section could stand.

| Cells | Temperature exponent | Rainfall coefficient |
|---|---|---|
| Correctly specified only | 82.5% | **71.5%** |
| Mis-specified | 83.0% | 71.1% |

Identical. The coverage failure is present when the model is exactly right; it is the
interval, not the model, that is wrong. With the correct specification *and* the full
series it is worse still — 68.6%. The objection is answered by measurement rather than by
argument, which is the only way it could have been.

### The check that cost me a defence

The decomposition's headline — three-quarters of the variation is analytical — was defended
by its stability across designs of 24, 48 and 144 cells. **That defence is weak and I should
have seen it.** The three designs are nested; all contain the observation model. Showing the
answer survives adding cells to one factor set is not showing it survives choosing another.

Redone properly, dropping each factor in turn and over all 63 subsets:

| Factor list | Within-outbreak share |
|---|---|
| All six | 76.5% [73.1–79.9] |
| Dropping the fixed parameters | 76.7% |
| Dropping the model structure | 73.6% |
| Dropping the rainfall lag | 73.0% |
| Dropping the temperature form | 71.4% |
| Dropping the fraction fitted | 70.5% |
| **Dropping the observation model** | **55.2%** |
| Median over all 63 subsets | 54.5% |
| Adversarial best case (parameters alone) | 2.3% |

Every factor but one costs at most six points. The observation model costs twenty-one. The
claim is **load-bearing on that single choice**, and the paper now says so in a table rather
than waiting to be told.

That makes the question "is the observation model genuinely a free choice?" not a nicety but
the hinge of the paper.

### The literature count, and why I did not do it

I had promised to read 30–50 methods sections and count. On starting, the better instrument
was already available and I had been citing it without reading its numbers: the systematic
review of 99 dengue prediction models reports its own breakdown.

* **Poisson regression 18.3%, linear regression 18.3%** of statistical models. The negative
  binomial does not appear as a category at all.
* **No validation reported 20.2%**, internal only 75.8%, external 5.2%.

Both of the paper's load-bearing factors are therefore documented as live choices by a
published review, not by my impression of one. Mechanistic fits, which that review does not
cover, commonly go further and use ordinary least squares — worse than Poisson, since it
assumes constant variance across an epidemic curve.

**A single-coder ad-hoc count of thirty papers would have been weaker evidence than a
published systematic review of ninety-nine**, and would itself have become a target. Not
doing it was the better call, and the promise is recorded here alongside the reason it was
broken.

One thing the review cannot supply: whether our two or three settings per factor span what
the literature does. They do not, and the limitations say so.

### What the reviewer's report could not fix

Affiliation and a domain co-author are not computable. Validation of the multiverse interval
on real data is impossible, because on real data there is no truth. Both stand as stated
limitations rather than as work outstanding.

---

## 2026-07-29 (later) — A second reviewer pass, and what it found

Amna asked for another strict review before submission. Four things came out of it, and
three were found by testing rather than by reading.

### The circularity objection, answered by generating from outside the family

Step 30 simulated from the model family it then fits, so a referee could dismiss its
coverage result as circular. `scripts/33_out_of_family.py` generates instead from **two**
climate-forced patches sharing one temperature series and one exponent, seeded 46 days
apart — the offset the two-patch model infers for Khyber Pakhtunkhwa in the appendix — and
fits every model with one patch. The mis-specification is not invented; the appendix
demonstrates it in these data.

| | Temperature exponent | Rainfall coefficient |
|---|---|---|
| *In family* — conventional | 82.9% | 71.2% |
| *In family* — multiverse | 100.0% | 97.5% |
| *Out of family* — conventional | 78.9% | **42.0%** |
| *Out of family* — **multiverse** | **95.0%** | **90.0%** |

**A 95% interval covering 42% of the time.** And the recommended fix degrades too: 97.5%
to 90.0%. Combining defensible analyses cannot recover from a structural error that all of
them share. It absorbs most of the damage — a 53-point shortfall becomes 5 — but it is an
improvement, not a solution, and the paper now says the reader should carry the 90%.

Reporting the degradation matters more than reporting the headline. A recommendation that
had only ever been tested where it works is not a tested recommendation.

### An internal contradiction, two paragraphs apart

Section 4.8 said, in consecutive paragraphs, that 39.8% of outbreaks contain a *confident*
dissenter that no band around zero can reach, and then that "where the direction does flip,
it flips on comparisons that were weak to begin with." The second sentence was written for
the five-factor design and survived the rewrite. A referee reading two paragraphs in order
would have found it.

Replaced with the version that holds: in three outbreaks in five all the dissent is weak and
a band reaches it; in the other two something dissents confidently and nothing reaches it,
but there the dissenter is nearly alone. **The surviving disagreement is common across
outbreaks and rare within them** — which is why 4.0% and 39.8% are both true.

### An overstatement, corrected against a better source

The discussion claimed "most studies use Poisson regression or comparably restrictive count
models", cited to a comparative paper. The systematic review gives the actual figures:
18.3% Poisson regression, 18.3% linear regression. That is **more than a third**, not most.
Corrected, and the weaker citation demoted to a supporting role.

Also fixed: two stale simulation numbers in the recommendations (47 points → 45, five in six
→ six in seven), and "negative in half of all outbreaks" where the figure is per fit.

### Two objections tested and dismissed

* **"Your NB standard errors ignore the fixed dispersion, so the under-coverage is
  self-inflicted."** Poisson fits, which do estimate their dispersion and inflate the
  standard error accordingly, cover *worse* (68.6% against 73.7% for the rainfall
  coefficient). Not self-inflicted.
* **"Your sixteen excluded windows are not random."** Mann–Whitney against the retained
  windows: p = 0.98 on length, 0.80 on cases, 0.79 on peak.

### Multiple comparisons in the thermal-regime analysis

Eighteen correlations were being reported with raw p-values and a flag on the small ones —
which is the undisclosed analytical choice this paper exists to measure. The false-discovery
rate is now controlled across the family. Seven reach p < 0.05 against about one expected,
so something is there; **one survives correction**, and its sign is the wrong way round for
the published hypothesis. Both facts are now in the same sentence rather than whichever one
suited.

### The question that found the gap: whose number is it?

Amna asked whether the paper could be improved further. Answering it honestly surfaced
something I had missed for three days: **every coverage result was computed on `a_temp` and
`a_rain`, and no dengue paper reports either.** What this field reports is R₀. The central
result was being measured on a quantity its audience does not look at.

Two additions followed, both from data already on disk.

**1. The variance decomposition, applied to R₀.**

| | Between outbreaks | Within (analyst) |
|---|---|---|
| R₀ under the **constant** model | 71.6% | 28.4% |
| R₀ under the **climate** model | 19.8% | **80.2%** |

The constant model's R₀ is mostly about the epidemic — which is what an estimate is for.
The climate model's is mostly about the analyst. The contrast survives all four treatments
of the long right tail (constant 27–28% throughout; climate 80–98%), so it is not an
artefact of outliers.

On a scale a reader can weigh: the median outbreak's R₀ ranges over **1.31** across the 144
analyses, against a median estimate of **1.29**. The spread the choices induce is about the
size of the thing being estimated.

This is the same fact the appendix reports from the other side. There, adding climate
parameters to a short series destroyed the identifiability of the population at risk that
the simpler model recovered to within a factor of two. Here the same addition destroys the
stability of R₀. **The climate terms are paid for in the precision of everything else, and
that cost is not reported because it is not looked for.**

**2. Coverage of the transmission coefficient, and so of R₀.**

Required recording the generating β₀, which step 30 had not done — so the study was rerun.

| Interval | Coverage |
|---|---|
| Conventional 95% | **78.0%** |
| Multiverse (Rubin) | 100.0% |

And the point estimates: the ratio of estimated to true transmission coefficient has median
**1.27** and interquartile range **0.74–2.59**. In the middle half of analyses the reported
R₀ is between three-quarters and two and a half times the value that generated the data.

Worth recording as a lesson about the reviewer exercise itself. Two rounds of adversarial
review found contradictions, stale numbers and one overstatement — all real, all worth
fixing. But the largest gap was not found by attacking the paper. It was found by asking
who reads it.

---

## 2026-07-29 (third review pass)

### The audit that gave false comfort

I wrote a check that extracts every percentage from the manuscript and looks for it among
every number in every result table. It reported all 117 clean. **It was worthless.** With
forty tables and thousands of values, almost any three-digit percentage matches something
by chance; the check had a false-negative rate near one.

Following it up found the real problem it should have caught: **five results in the paper
were produced by ad-hoc commands and by no pipeline step.** The R₀ decomposition, the R₀
coverage, the leave-one-country-out, the clustered bootstrap, and the how-many-analyses
study. Anyone re-running the repository could not obtain them, which breaks the one claim
this paper cannot afford to get wrong.

All five now live in steps 30 and 32 and write tables (`35_r0_decomposition.csv`,
`36_leave_one_country_out.csv`, `37_analyses_needed.csv`). Re-running them changed several
figures in the third decimal — different random draws — and the manuscript was moved to the
scripts' numbers rather than the other way round.

The general lesson: a verification whose failure mode is "passes anyway" is worse than none,
because it stops you looking.

### The introduction was still describing the old paper

The contributions had moved three times and the introduction had not. It promised "how
often the verdict changes" and "a decision rule that restores stability", neither of which
is now the main result, and the second of which is an overclaim. Replaced with four numbered
findings and one honest qualifier.

A reader forms their expectation in the first page. Leaving it stale is not a cosmetic
error.

### Two robustness questions a referee would ask, answered

**"Your outbreaks are not independent."** 221 outbreaks from 85 reporting units in 33
countries; three countries supply 45%. Resampling whole countries rather than outbreaks
widens the interval from [73.1, 79.9] to [72.6, 81.1] and does not move the point estimate.

**"Is this a fact about Latin America?"** Removing each well-represented country in turn
gives 76.0–78.6%. Removing Bolivia, Nicaragua and Mexico together, leaving 121 outbreaks,
gives **78.5%**. National series alone 79.1%, subnational alone 75.4%. Every subset lands
within three points of the whole.

And one that pushes the other way, now stated in the paper: every combination within an
outbreak is warm-started from one anchor and shares one dispersion estimate, so the 144
analyses are more alike than 144 independent analysts would be. The optimiser check bounds
that at 0.6% of verdicts, but the sign is clear — the analytical share is if anything
understated.

### Eight analyses are enough

The recommendation had been "combine the factorial", which costs 288 fits and would be
admired rather than used. Drawing m analyses at random and combining only those:

| Analyses | 2 | 4 | **8** | 16 | 72 |
|---|---|---|---|---|---|
| Rainfall coefficient | 91.0% | 96.2% | **97.3%** | 97.6% | 97.5% |
| Transmission coefficient | 91.7% | 97.1% | **98.9%** | 99.9% | 100.0% |

**Four reach nominal coverage; eight are indistinguishable from seventy-two.** The interval
also stops narrowing after eight. That turns the recommendation from a research exercise
into something a working epidemiologist can do in a few minutes on one processor — which is
the difference between a paper that is agreed with and one that changes practice.

### The question the paper had never answered

Three review passes and it had not occurred to me that a reader finishing this paper would
ask the obvious thing: *so does climate drive dengue transmission or not?* The paper told
the field its methods were unreliable and stopped. That reads as a refusal, and it invites
a misreading far worse than the refusal — that this work shows climate does not matter. It
shows nothing of the kind.

Applying the paper's own recommendation to the paper's own data — negative binomial only,
abstain where |ΔAIC| < 4:

| | Outbreaks | Share |
|---|---|---|
| **Climate forcing supported** | **118** | **53.4%** |
| No climate forcing | 19 | 8.6% |
| Analyses still disagree | 49 | 22.2% |
| Inconclusive: nothing speaks | 35 | 15.8% |

**A clear verdict for 137 of 221, and 86% of those support climate forcing.** Where these
data can answer the question, the answer is usually yes.

And the other half of it, which is the number to carry away: **in 38% of outbreaks the
question cannot be answered at all** from one season by this method. That is a result, not
a gap in the reporting.

This belongs near the front of the discussion and in the abstract, and it is now in both.
Without it the paper is a complaint; with it, it is a measurement that happens to include a
warning.

### The abstract, three times

It reached 339 words as findings accumulated — over every limit worth targeting. Compressed
to 258 by keeping five results and dropping the factor-by-factor breakdown and the
simulation error rates, both of which are better read in the body than skimmed in the
abstract.

---

## 2026-07-29 (fourth review pass) — what was written but never shown

This pass looked at the parts of the paper that are not prose.

### The main results figure was not in the paper

Twenty figures exist in `results/figures`; the manuscript used four. Among the sixteen
unused were **`12_global_robustness.png`, which is the study's main results figure**, and
**`19_decomposition_robustness.png`, which supports the single most attackable claim in it**.
A paper whose headline is a variance decomposition had no picture of that decomposition's
robustness, and a referee would have had to take the table on trust.

Both are now in. The temperature panel of the first was also redrawn: the Brière exponent
has a long right tail running to twenty, and plotting the full range compressed the part of
the axis that carries the point — whether the estimate falls on the negative side of zero.
Clipped to ±3.

### The model had no equations

Twenty-three pages about how sensitive these fits are to their specification, and the
specification was described only in words. A reader wanting to know what was fitted had to
read the source. The host–vector system, the two alternative structures, the two forms of
β(t) and the observation equation are now stated.

### The inference had no description either

The paper's central claim is about **intervals**, and it did not say how an interval was
computed. It now states the integrator and step, the parameters estimated, the
deviance-residual formulation that lets a least-squares routine perform Poisson or
negative-binomial maximum likelihood, the transforms, the profiled dispersion, the penalty
that stops the optimiser following a divergent proposal downhill, and — the part that
matters — that the standard errors are asymptotic from the Jacobian, inflated by the Pearson
dispersion, delta-method transformed. Those are the intervals whose coverage the paper
measures, so the reader has to be told what they are.

### Two tables no script produces

`13_global_robustness_pilot.csv` and `18_band_placement.csv` were left in `results/tables`
from work that has since been replaced. Nothing in the repository writes them, so anyone
re-running the pipeline would find files that cannot be reproduced. Deleted.

---

## 2026-07-29 (fifth review pass) — an assertion that was backwards

Reading the coverage section as a referee, one sentence stood out because it had no evidence
behind it: *"Which eight matters less than that they differ ... the observation model, being
the largest lever, is the one we would vary first."* Plausible, and asserted.

It is wrong, and wrong in the interesting direction.

| Eight analyses sharing… | Rainfall coef. | Temperature exp. | Transmission coef. |
|---|---|---|---|
| nothing (spanning draw) | 97.5% | 99.9% | 99.2% |
| one observation model | 97.1% | 98.9% | 98.6% |
| one model structure | 97.9% | 99.7% | 96.5% |
| one fixed parameter set | 97.4% | 99.7% | 98.6% |
| one fraction of the series | 97.0% | 99.6% | 99.2% |
| **one rainfall lag** | **88.7%** | 97.7% | 97.8% |

**Holding the observation model fixed costs the interval nothing.** Eight analyses that all
use the same likelihood cover as well as eight that do not. The only choice that cannot be
held fixed is the **rainfall lag**, and only for the **rainfall** coefficient, which loses
nine points.

The rule that fits is narrower and far more useful than the one I assumed: *vary the choices
that construct the covariate whose coefficient you are reporting.* The lag is what makes the
rainfall covariate; fix it and that coefficient's analytical uncertainty never enters the
arithmetic, however much else is varied.

It also separates two things the paper had been treating as one. **The observation model is
the largest lever on the verdict and nearly irrelevant to the interval.** A practitioner
reporting a rainfall effect needs three lags and can fix everything else; one reporting a
verdict needs both observation models. Different quantities, different requirements — which
the paper now says instead of offering a single instruction for both.

Recorded because the sentence had survived four passes. Prose that sounds like a finding and
is not one is harder to see than a wrong number, and this one was load-bearing on the
practical recommendation.

---

## 2026-07-29 (sixth review pass) — the sentences that sounded like findings

Pass five found an untested assertion that turned out to be backwards, so this pass looked
for the rest of them: sentences that read like results and carried no number.

### "A monotone term must return a negative coefficient"

The paper explained the arbitrary sign of the log-linear temperature coefficient by saying
that a monotone term fitted to a season whose temperature falls while cases rise *must*
return a negative value, so the sign is a property of the calendar. Testable against each
outbreak's growth-phase temperature trend, and never tested.

* Correlation of the estimate with the trend: **ρ = +0.25** (p = 1×10⁻⁴) — the predicted
  direction.
* Temperature falling during growth: coefficient negative in **52%**.
* Temperature rising: negative in **32%**.

So the direction holds and the strength does not. "Must" was wrong; the calendar **tilts**
the sign without setting it. That leaves most of the sign unexplained, which is arguably a
worse fact for the log-linear parameterisation than the original claim — a coefficient whose
sign is part calendar and part unaccounted for is not measuring anything. The paper now says
that.

### "Dispersion is two orders of magnitude larger"

The mechanism paragraph said observed dispersion is "far larger" than Poisson; the discussion
said "two orders of magnitude". Neither carried a measurement, and the second is **wrong**.

Measured model-free — residuals about a five-week centred moving average, so the figure does
not depend on any transmission model being right:

| | Index of dispersion |
|---|---|
| median | **5.7** |
| quartiles | 2.6 to 19.1 |
| 90th percentile | 59.1 |
| above 20 | 24% of outbreaks |

Six times at the median, not a hundred. Two orders of magnitude is true only of the extreme
tail. Both sentences now carry the measured values, and `scripts/34_dispersion.py` produces
them.

The concrete version is better anyway: a week of 2,000 cases is treated by Poisson as
accurate to ±45, where the honest figure at the median dispersion is ±107.

### The pattern across passes five and six

Three overstatements have now been found in prose that had survived several readings —
"most studies use Poisson", "must return a negative coefficient", "two orders of magnitude"
— and one assertion that was exactly backwards. All four were in sentences explaining *why*
a result holds rather than reporting the result. Explanations are where a paper is least
guarded, because the numbers around them are right and the sentence sounds like it follows.

---

## 2026-07-30 (seventh review pass) — the paper's central explanation was wrong

Passes five and six found overstatements in explanatory sentences. This pass tested the
biggest explanation in the paper: **why** the observation model dominates.

The paper said — and every draft since the first has said — that Poisson treats ordinary
noise as a large residual, so the fit chases individual weeks, and the climate covariates
supply the flexibility to do it. It is a good story. It makes a prediction: the climate
coefficients should be pushed harder under Poisson.

**They are not.** Paired within outbreak and within every other choice, across 15,912 pairs:

| | median under NB | median under Poisson | larger under Poisson |
|---|---|---|---|
| \|a_temp\| | 0.092 | 0.097 | **40%** of pairs |
| \|a_rain\| | 0.079 | 0.073 | 46% of pairs |

Forty per cent is *less often than chance*. The climate model does the same thing under
both likelihoods.

### What is actually happening

The climate model nests the constant one, so it always fits at least slightly better. The
only question is whether the improvement clears two parameters' worth of penalty. Without a
dispersion parameter to absorb the residual, the same improvement buys far more
log-likelihood:

| | |
|---|---|
| median ΔAIC under NB | **−1.6** |
| median ΔAIC under Poisson | **−68.9** |
| ratio where the sign agrees | **13×** (quartiles 3.8 to 67) |
| \|ΔAIC\| larger under Poisson | 90% of pairs |

And the asymmetry is near-total: climate forcing wins under both in 56.4% of pairs, **under
Poisson alone in 33.6%, and under the negative binomial alone in 0.6%.**

**Poisson does not make the model chase noise. It makes the criterion generous, and generous
in one direction only.** That is a better explanation than the one it replaces — it explains
why the parameters do not move, which the old story could not — and it is measured rather
than told.

### The pattern, now four passes deep

Every error found since pass five has been in a sentence explaining *why* a result holds:

* "most studies use Poisson" — a third, not most
* "must return a negative coefficient" — tilts the sign, does not set it
* "two orders of magnitude" — six times at the median
* "vary the observation model first" — varying it changes nothing for the interval
* "Poisson makes the fit chase noise" — it makes the criterion generous instead

Five for five. The numbers in this paper have been checked to death and the *explanations*
had never been checked at all, because an explanation sitting between two correct numbers
reads as though it were derived from them. It is not, and none of these was.

---

## 2026-08-03 — Pass 8: the sentences get a test of their own

Pass 7 ended by noting that every error found since pass five sat in an explanatory
sentence. Pass 8 took that as the search rule and swept the rest of the prose. It found
five more, one of which was a contradiction the paper had been carrying for two passes,
and one of which turned into the best result in the paper.

### 1. The discussion still carried the mechanism pass 7 refuted

Pass 7 corrected the results section. The discussion, twelve pages later, still read
"rewards any model flexible enough to chase noise, and climate covariates supply exactly
that flexibility" — the claim that had just been measured false. A paper that argues one
mechanism in section 4 and the refuted one in section 5 is worse than one that never
corrected anything, because the reader cannot tell which the authors believe.

The same slip had happened once before, in the other direction: "most studies use Poisson"
was corrected in the discussion in pass 5 and left standing in the results, where it sat
contradicting the corrected version six pages later. Both are now fixed, and both are in a
test — `test_withdrawn_explanations_do_not_return` — that fails if any of the five
withdrawn phrases reappears. Correcting a claim in one place is not correcting it.

### 2. "Most of the gap is not detection" — two errors in one sentence

The paper wrote: *the simulation shows 85.7% of the Poisson endorsements would occur on
data with no effect at all. Most of the gap between the two observation models on real
data is not detection.*

85.7% is the false-positive **rate under a null truth**. It is not a share of endorsements,
and reading it as one is the conditional-probability inversion this paper's own subject
matter is about. The second sentence then rests on the first.

Measured properly, pairing Poisson against NB within outbreak and within every other choice:

| discordant cell (Poisson endorses, NB does not) | share of pairs |
|---|---|
| simulated, no climate effect | **67.1%** |
| simulated, real effect present | **23.3%** |
| real data | **33.6%** |

The cell *is* diagnostic — three times as common under the null. But the real rate sits
nearer the arm that **contains** an effect. A mixture read off the three puts the null
fraction near a quarter, not a majority. So the sentence was not merely imprecise, it was
backwards, and it is withdrawn on the same grounds as the inference already withdrawn two
sections earlier. The paper now reports the three rates and declines the inference.

### 3. "The sign is unambiguous" — the one direction never checked

On warm starts: *the 144 analyses are more alike than 144 independent analysts would be
... but the sign is unambiguous: the analytical share reported here is if anything an
underestimate.*

The optimiser check had refitted seven outbreaks from ten cold starts two passes earlier.
Nobody had asked it this question. Pairwise disagreement is **33.4%** as run and **32.9%**
cold: cold starts give slightly *less* disagreement, not more. Seven outbreaks settle
nothing about the size, and the check varies the warm start but not the shared dispersion
estimate, so half the mechanism is still untested. What it does settle is that the word
"unambiguous" was doing work no measurement supported.

### 4. The limitation named the minority mechanism

*Requiring an unbroken weekly run excludes settings with intermittent reporting, which are
likely to be poorer and more heterogeneous than those retained.*

Step 36 measures it. Of 88 countries in the weekly subset, 34 contribute a usable window.
Of the 54 excluded:

| binding requirement | countries |
|---|---|
| no gap-free run of 30 weeks | 17 (31%) |
| long enough, wave too small | 22 (41%) |
| long and large, rejected on shape | 15 (28%) |

**Sixty-nine per cent report continuously enough.** The dominant selection is on outbreak
size and shape, not on reporting continuity — and it points the *opposite way* from the
concern as written. Large, sharply peaked, single-wave epidemics are exactly what these
models handle best, so the instability measured here is what survives on the most
favourable data available. That is a stronger position than the paragraph claimed, arrived
at by discovering the paragraph was wrong.

Wealth is not recorded in OpenDengue or NASA POWER. That clause is dropped rather than
propped up with an imported classification.

### 5. "Because the analysts differed than because the places did" — and what replaced it

This one was not a loose sentence. It was the paper's reading of its own headline, and it
does not follow from it.

The one-way split reports *between* as the variance of per-outbreak means and *within* as
the mean of per-outbreak variances. Three-quarters within is correct. But "within" is a
remainder, and a remainder is not an attribution: it contains a **main effect of the
analysis**, which standardising conventions would remove, and an **outbreak-by-analysis
interaction**, which it would not. The design is completely crossed and complete by
construction, so the partition is exact and there was never any reason not to do it:

| design | outbreak | analysis | interaction |
|---|---|---|---|
| 24 analyses | 21.7% | 19.6% | 58.7% |
| 48 analyses | 20.3% | 21.2% | 58.5% |
| **144 analyses** | **23.5%** | **15.9%** | **60.6%** |
| 95% CI, countries resampled | 18.6–27.7 | 11.9–22.9 | 56.1–63.8 |

**The largest term belongs to neither.** Three-fifths of the variation in whether an
outbreak is judged climate-driven is specific to the *pairing* of a method with an
epidemic. It is the most stable term across designs, it survives dropping each dominant
country (59.7–62.1%), and with one deterministic fit per cell the only thing it could be
confounded with is optimiser noise, which the cold-start check bounds at 0.6% of verdicts.

Checked in the form the discarded sentence actually claimed:

| two findings that share... | disagree |
|---|---|
| the outbreak, not the analysis | 29.8% |
| the analysis, not the outbreak | **32.7%** |
| neither | 38.9% |

Standardising the analysis removes 6 points; standardising the place removes 9. **The gloss
was the wrong way round.** A variance of means and a mean of variances are not comparable
in that direction, and nobody noticed because the conclusion was the one we wanted.

Three things follow, and they are why this pass produced a result rather than a correction.

**It is what a specification curve cannot show.** The curve orders specifications by their
average result across datasets — that *is* the main effect — and averages the interaction
away. Here the term on display is a quarter of the term it hides. This is a limitation of
the method Simonsohn et al. describe rather than of our use of it. Searching found ANOVA
over analytic *factors* in the multiverse literature, but not the dataset-by-specification
interaction; the neuroimaging multi-pipeline datasets have the crossed design in which the
same quantity could be computed, and we have not found it computed.

**It rescues the headline from its own objection.** "Three-quarters is analytical" invites
"so standardise the analysis". The answer is now quantitative: standardising reaches 16%.

**It reorders the recommendations.** The paper listed four in "descending order of how much
they are worth" and then said the fourth was the one it would keep if only one survived —
an ordering contradicting itself inside one paragraph. Recommendations 1–3 act on the main
effect. Only the fourth, the per-outbreak interval, touches the interaction. That is the
order, and it is not the order we thought of them in.

### The pipeline grew a step whose job is the prose

`35_claims_audit.py` exists because `test_paper_consistency.py` cannot see half of a paper.
It refits nothing, reads the stored tables, and checks the interpretive claims: the two-way
decomposition, the discordant-cell rates, the warm-start direction. `36_selection_audit.py`
does the same for the selection rule. Both write tables the consistency suite now asserts
the paper quotes, and the suite also checks its own quoted size, which had read 157 for
long enough that the real figure was half again as large.

Eight passes. Six errors of arithmetic or code, all found by tests. **Ten errors of
interpretation, every one found by asking a stored table a question the paper had answered
from the armchair.** The tables had the answers the whole time.

---

## 2026-08-03 — Pass 9: the citations, which nothing had ever checked

Eight passes had checked the numbers, and pass 8 started checking the sentences.
Nothing had checked the **references**. Every entry in the bibliography was verified
against the published record, and every direct quotation word for word against its source.

### Two quotations were paraphrases wearing quotation marks

The discussion attributed to Khamthong and Phramrung the sentence *"stronger marginal
climate–dengue associations do not necessarily translate into superior out-of-sample
predictive performance"*. Their paper does not contain it. What it contains is
"forecasting performance is not determined solely by the strength of marginal
climate–dengue associations", and, of the two stations, that their predictive performances
"were statistically equivalent".

The meaning survives; the quotation does not. A referee who checks one quotation and finds
it invented does not check the second — they stop reading. Both are now the source's own
words.

The Leung review was quoted as "...inadequate in many existing prediction models". The
abstract says "in many **of the** existing prediction models". Small, and still a
misquotation.

### One citation had the wrong journal and no authors

`nbreview` was a prose entry — no author list, two works bundled into a sentence, one of
them a Research Square preprint cited only to support a phrase ("remain limited in the
literature") that could not be verified because the preprint is behind a 403. The
peer-reviewed half of it was attributed to **Gene Reports**; it is in the **Journal of
Infection and Public Health** 18(11), 102906 (2025), by Al-Manji, Al Wahaibi, Al-Azri and
Chan.

Replaced with that entry, properly formatted, supporting only what its abstract states: in
a head-to-head comparison on mosquito-borne outbreak data the negative binomial is the
better model. The unverifiable quotation is gone, and the claim it propped up was already
carried by the Leung review's own 18.3% / 18.3% breakdown, which is verified.

Also checked and correct: OpenDengue is CC BY 4.0 for v1.3 as the paper states; the Leung
percentages (18.3, 18.3, 20.2, 75.8, 5.2) are all as published; Kirk et al.'s direction
("associations strongest when temperature variation and population density were high")
matches how the paper describes it; the 2022 correction to Mordecai et al. touches two
supplementary figures and a paragraph of supplementary results, not the thermal trait
values this study fixes its Brière limits from.

`test_every_quotation_has_been_checked_against_its_source` now holds the verified list. Any
new quotation fails the suite until someone checks it. The test also asserts it found at
least as many quotations as it knows about, because a regex that silently matches nothing
would pass forever.

### The paper misdescribed the two ends of its own specification curve

*"The least credulous --- negative binomial, Brière, a seven-week lag --- endorses it in
39.8%."* The least credulous specification uses a **three**-week lag. The description was
left over from the 48-cell design.

Regenerating the figure and reading what it actually plots gives something better than a
correction. The extremes are:

| | observation | temp form | lag | fraction | structure | params |
|---|---|---|---|---|---|---|
| most credulous, 95.9% | poisson | log-linear | 3 | 0.75 | seir | central |
| least credulous, 39.8% | **nb** | **briere** | 3 | 0.75 | seir | central |

**The two ends of a 56-point spread differ in exactly two of the six choices, and both are
statistical conventions rather than biology.** That is the paper's thesis in one row, and
it had been sitting in the output of step 25 unread. The figure caption also still said
"all 48 analyses of a single outbreak" for a panel showing 144.

`mordecai` was in the bibliography and cited nowhere, while the Brière thermal limits it
supplies were stated without attribution. Now cited where they are used.

### The framing paragraph gained the methodological point

A specification curve orders specifications by their average result across datasets, which
is a main effect. With many datasets each analysed every way, the interaction is separately
estimable — and pass 8 found it to be the largest term. That belongs in the paragraph
relating this study to Steegen et al. and Simonsohn et al., not only in the results, because
it is a statement about the method rather than about dengue: a main effect can be legislated
away by convention and an interaction cannot.

Paper is 26 pages, 246 tests, compiles clean.

---

## 2026-08-03 — Pass 10: the interaction is not arbitrary, and mostly it is

Pass 8 left the paper with a large unexplained term. A referee reading "three-fifths of the
variation is an outbreak-by-analysis interaction" asks the obvious next question, so this
pass asked it first.

### Reproducibility, checked rather than claimed

Every analysis-only step — 17, 18, 19, 20, 22, 32 — was rerun against the stored tables and
every regenerated file is **byte-identical**. The claim that the pipeline reproduces the
paper is now a thing that was tested rather than a thing that was said.

Step 25 was rerun too, and its figure is drawn from the 144-cell design as the text
requires; the *caption* had said 48 and the prose describing the least credulous
specification was left over from the smaller design. Both were fixed in pass 9.

### One number was being compared across two designs

The optimiser check reports 0.6% of verdicts changing "against 8.3% for the weakest genuine
factor and 38.1% for the strongest". Those are the **48-combination** figures; the paper's
own main table, four pages earlier, says 3.1% and 34.3%. Both are right — the check was run
on the smaller design and the comparison has to be made within one design — but a reader
sees two pairs of numbers for the same quantity and no explanation. Now stated.

### The question: can an analyst see it coming?

For each outbreak and each factor, the effect is the change in the share of analyses
endorsing climate forcing between two levels of that factor, averaged over every combination
of the rest. That is the outbreak's own slice of the interaction. One hypothesis had a
direction fixed before testing: Poisson misbehaves because it denies overdispersion, so its
effect should be largest where dispersion is worst. Dispersion is the model-free index from
step 34 — a five-week moving average — which enters no cell of the factorial.

**It holds.** ρ = +0.31, p = 2×10⁻⁶. Below median dispersion, switching to Poisson raises
the endorsement rate by 0.25; above it, by 0.42.

This matters beyond the correlation. Pass 7 established the mechanism (Poisson does not
enlarge the coefficients, it enlarges the criterion) from the factorial's own ΔAIC values.
This is the same conclusion reached from a quantity computed on the raw counts, with no
model fitted. Two routes, no shared arithmetic.

**And it stops there.** Regressing each factor's effect on dispersion, total cases, series
length and mean weekly count:

| factor | variation explained | sd of the effect |
|---|---|---|
| **observation model** | **34.3%** | 0.291 |
| fixed parameters | 3.5% | 0.029 |
| temperature form | 2.0% | 0.152 |
| model structure | 1.6% | 0.111 |
| fraction of series fitted | 1.0% | 0.238 |
| rainfall lag | **0.1%** | 0.185 |

Median 1.8% across the six. Of 24 correlations, 5 reach p < 0.05 and 3 survive FDR control,
all for the two factors with a mechanism behind them.

The rainfall lag is the cleanest statement in the pass: it flips one verdict in six, its
effect varies across outbreaks as much as the temperature form's, and **essentially none of
that is predictable from anything visible before fitting**. An analyst cannot know whether
their lag choice mattered without trying the others — which is, independently, exactly what
the coverage study concluded two sections earlier for a different reason (fixing the lag
costs the rainfall coefficient nine points of coverage). Neither result was designed to test
the other.

### A bug in this pass, caught by the result looking too good

The first run printed "5 reach p < 0.05; **24 of 24** survive false-discovery-rate control",
which is not a thing that happens. Benjamini–Hochberg enforces monotonicity by a running
minimum taken from the **largest** p-value downward; taking it forward from the smallest lets
one genuine finding certify every null beside it. With p = 2×10⁻⁶ at the top of the family
that is precisely what happened.

Step 31's implementation was correct, and its result — one of eighteen surviving — was never
affected. But two copies of a correction, one right and one wrong, is the same shape as the
`between(1, 23)` bug that started this log: a definition living in a script instead of in the
package. `benjamini_hochberg` is now in `robustness.py`, both steps call it, step 31's output
is unchanged to the digit, and seven tests cover it — including one that asserts the *wrong*
direction would have passed the whole family, so the failure mode is documented rather than
merely fixed.

Ten passes. The instrument that keeps finding things is the same one every time: take a
sentence that sounds like it follows from a table, and ask the table.

Paper is 27 pages, 258 tests, compiles clean.

---

## 2026-08-03 — Pass 11: the partition, applied where the answer is known

Pass 8 introduced the three-way partition and applied it to one thing. Two places where it
obviously belonged had been left with the old one-way split, which is the same inconsistency
of rigour that pass 8 was created to fix.

### R₀: no convention would fix it

The paper reported that adding climate covariates moves R₀ from mostly-between-outbreak to
mostly-not. "Mostly not" is a remainder, and the practical question — *could the field agree
on a method and be rid of this?* — depends entirely on which part of the remainder it is.

| (log scale) | outbreak | analysis | interaction |
|---|---|---|---|
| R₀, constant model | 71.6% | 5.3% | 23.0% |
| R₀, climate model | 22.5% | **8.4%** | **69.1%** |

**The analysis main effect is about a twentieth under either model.** There is no
systematically high-reading method to stop using. What the climate terms do is move variation
out of the outbreak and into the interaction, 23% → 69%. The instability in the number this
field publishes is not a house style anyone could adopt their way out of.

### The simulation has a fingerprint, and the real data matches one arm

The paper contained a paragraph headed "An inference we decline to draw": two attempts to ask
which simulated arm the real outbreaks resemble had pointed in opposite directions, and one
had already been withdrawn. The reason both failed is now visible — **the arms differ in
three components and each of those statistics compressed them to one.**

| | outbreak | analysis | interaction |
|---|---|---|---|
| simulated, no climate effect | 10.9% | **47.0%** | 42.1% |
| simulated, real effect present | 28.6% | 10.3% | 61.1% |
| **real outbreaks** | **23.5%** | **15.9%** | **60.6%** |

The **analysis main effect** is the discriminating component, and the mechanism is exactly
what it should be: with nothing to detect, the answer is whatever the convention gives, and
the convention does not vary by window — so it becomes the largest term at 47%. Where a real
effect is present it collapses to 10%. The real outbreaks read 15.9%, four times away from
the null arm, and land within five points of the effect arm on the other two components as
well.

Stated as evidence, not proof: the simulated effect size is drawn from these same fits, the
arms are not a calibrated mixture, and matching a fingerprint does not put a signal in every
outbreak. What it does rule out is the reading the withdrawn sentence pointed toward — that
this literature carries the signature of detecting nothing. It agrees with the paper's own
answer under its own rule (137 clear verdicts, 86% supporting climate forcing) from a
completely different direction.

A statistic that behaves this way is also usable by anyone else: a multiverse whose
disagreement is mostly a main effect of the analysis is a multiverse that may have no signal
in it.

### Two guards rather than one

`test_decomposition_under_truth_is_quoted` checks the nine figures. A second test checks the
*property*: the null arm must have the larger analysis main effect, and the real data must
sit nearer the effect arm on it. Without that, a rerun could invert the finding while every
quoted number still matched its table — which is precisely how this project's two worst bugs
behaved.

Paper is 28 pages, 263 tests, compiles clean.

---

## 2026-08-03 — Pass 12: nine references

Eleven passes had gone into whether the paper's claims were true. None had asked whether the
paper looked like a paper. It had **nine references**, four of them methods classics, while
making repeated claims about what "this literature" does.

That is not a stylistic complaint. A manuscript that critiques a field and cites nine works
reads to an editor as one that has not read the field, and it is the cheapest possible
reason for a desk rejection — cheaper than disagreeing with anything in it.

### Twelve added, every one verified before it was written down

Each was checked for author list, journal, volume, issue, pages and year against the
published record, not recalled:

| added | for |
|---|---|
| Silberzahn et al. 2018, AMPPS 1(3) | 29 teams, one dataset, effect sizes 0.89–2.93 |
| Botvinik-Nezer et al. 2020, Nature 582 | 70 teams, one fMRI dataset |
| Wang et al. 2024, J Clin Epidemiol 168 | specification curve in epidemiology — 1,208 analyses of red meat |
| Brière et al. 1999, Environ Entomol 28(1) | the thermal response this model uses, cited nowhere before |
| Mordecai et al. 2019, Ecol Lett 22(10) | the thermal-biology synthesis behind the limits |
| Andraud et al. 2012, PLoS ONE 7(11) | dengue model structures reviewed; no settled convention |
| Lowe et al. 2021, Lancet Planet Health 5(4) | rainfall lags with a mechanism and no agreed value |
| Bhatt et al. 2013, Nature 496 | under-ascertainment, where ρ is fixed |
| Lloyd-Smith et al. 2005, Nature 438 | why counts are overdispersed in the first place |
| McCullagh & Nelder 1989 | the standard remedy for it |
| Rubin 1987 | the rules the interval is built from |
| Raue et al. 2009, Bioinformatics 25(15) | profile-likelihood identifiability, used in the appendix |

Twenty-one references, every one cited, every citation defined — both now checked by
`test_every_reference_is_cited_and_every_citation_defined`, after `mordecai` was found
sitting in the bibliography while the thermal limits it supplies were stated without
attribution.

The three many-analyst citations do more than pad. They put this study in a lineage the
reader may already accept — Silberzahn's 29 teams, NARPS's 70 — and they make the
contribution legible: those designs vary the analyst on **one** dataset and so cannot
estimate a dataset-by-analysis interaction at all. This one can, and it is the largest term.

### The two things a submission still needed

**An Author Summary.** PLOS NTD and PLOS Global Public Health both require 150–200
non-technical words and the paper had none. Written (192 words, `paper/author_summary.tex`),
kept in its own file so the manuscript still compiles for venues that do not want it.

**An honest read of the length.** 29 pages, ~14,400 words, 8 figures, 20 tables. That is
roughly double a typical Epidemics paper. `docs/SUBMISSION.md` now carries the numbers per
venue and an ordered list of what to move to supplementary if an editor pushes back — the
appendix case study first, the 63-subset tables second — together with the four things that
must not be cut. Deciding that in advance is easier than deciding it under a revision
deadline.

The abstract is 294 words against PLOS's 300-word limit, which is fine and has no room in
it. Worth knowing before adding anything.

### Where this leaves the paper

Twelve passes. The first four fixed arithmetic. Passes five to eleven fixed explanations,
which turned out to be where nearly everything wrong actually lived — ten of them — and
which produced the paper's two best results as by-products: the three-way partition came out
of correcting a gloss, and the simulation fingerprint came out of a paragraph that had
admitted defeat. This pass fixed the parts that are not claims at all.

What is left is not something more review can supply: an outside reader, a repository URL,
and a DOI.

### Pass 13 — one last reconciliation

The new three-way $R_0$ table reads 22.5% between outbreaks where the one-way table four
paragraphs above reads 19.8%. Both are right: the two-way version needs a balanced design,
so it drops the two outbreaks whose $R_0$ is not positive and finite under every one of the
144 analyses. A reader comparing the two tables would have found a contradiction and had no
way to resolve it. Now stated in the caption, with the raw-scale figures beside it.

The quotation guard then failed on this very edit — it flagged the scare-quoted word
"outbreak" in the new caption as an unverified quotation. Which is the test doing exactly
its job on the first change made after it was written.

---

## 2026-08-09 — medRxiv declined it, and the reason was not the paper

Submitted 8 Aug (MEDRXIV/2026/360009). Declined the next morning:

> "medRxiv requires authors to have an organizational affiliation. It is necessary for
> submissions to be associated with an organization that provides oversight of research
> activities so that it can adjudicate any ethical issues/disputes that arise."

Every declaration was accepted. The manuscript was never in question. The submission
checklist in this repository had asserted — in writing, without a source — that
"Independent researcher is accepted; medRxiv does not require an institutional address".
That was invented, and four days were spent on the strength of it.

The lesson is the same one this log keeps recording in different clothes: **a plausible
claim sitting next to correct ones is not thereby correct.** It applied to the paper's
explanatory sentences through seven review passes, and it applied here to a claim about a
submission policy. The fix is the same too — check it against the source rather than
against how confident it sounds.

### What the PDF check found before any of this

Two faults reached the compiled PDF and had survived thirteen review passes, both found by
opening the file and looking at it:

* The model's seven differential equations were an `align*` block whose row breaks had lost
  a stroke: `\` became `\` in three places. The rows never broke, the display ran off the
  page, and LaTeX reported neither an error nor an overfull box.
* 138 em dashes in the prose, one every 94 words.

Neither was findable by any test then in the repository. Both now have one:
`test_no_line_ends_in_a_lone_backslash` and the control-character sweep, the second
verified by reintroducing the original corruption and confirming the failure.

Also set the manuscript in a Times-like face. Computer Modern is the LaTeX default and
universal in mathematics; it is rare in epidemiology, which is the readership. A typeface
that reads as "not from this field" is a free signal to give away.

### Where it goes now

MetaArXiv. Verified against the Center for Open Science's own announcement rather than
assumed: the OSF generalist server was suspended in August 2025 and will not return, but
fourteen community-run servers remain operational and MetaArXiv is one. No institutional
affiliation is required.

It may be the better home regardless. The paper's most novel claim is not about dengue: it
is that a specification curve shows the analysis main effect and averages away the
dataset-by-analysis interaction, which here is four times larger. That is metascience, and
MetaArXiv's readership is the multiverse and reproducibility community.

The cost is that epidemiologists will not browse it, and the answer to that is the outreach
emails rather than the server.

### The affiliation problem is now the critical path

It has blocked the work twice: it ended the medRxiv submission, and it is the largest single
factor against journal acceptance. Every route out of it — a co-author, a former department,
a Pakistani dengue researcher, an arXiv endorser — reduces to the same requirement: **one
person who knows the field agreeing to look at this.**

No further review pass changes that, which is worth stating plainly in a log that has
recorded thirteen of them.

---

## 2026-08-12 — MetaArXiv declined it too, and the lesson is about titles

Moderator feedback, in full: "Outside of the scope of this preprint series."

The reasoning that sent it there was that the paper's most novel claim is metascientific.
That reasoning was mine, and it ignored how moderation actually works: a volunteer sees
"climate-driven dengue transmission... 221 outbreaks" in the title and triages it as
epidemiology in the first ten seconds. The metascience sits in the abstract's third
sentence. Nobody doing scope triage reads to the third sentence.

Two venues, two declines, neither about the manuscript: medRxiv on policy (affiliation),
MetaArXiv on scope (title). Both misjudgements were mine and both were of the same kind
as the paper's own subject — a plausible belief acted on without checking how the
decision-maker on the other side actually decides.

The fallback written into SUBMISSION.md on day one now runs: a dedicated Zenodo record
for the manuscript (no moderation exists there, DOI is immediate) and direct journal
submission, which never required a preprint in the first place. The citable-DOI purpose
and the under-review status arrive by separate roads, which is where they always were.
