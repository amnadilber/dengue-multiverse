# Model specification

## Why a host–vector model and not SIR

Dengue does not pass directly between humans. The virus travels from an infectious
person to a feeding *Aedes* mosquito, replicates inside it for one to two weeks, and only
then can be transmitted onward. A single-population SIR model has no way to represent
that delay, and fitting one to dengue data forces the transmission rate to absorb a
mechanism it does not contain. The model here tracks both populations.

**Humans (SEIR).** The exposed compartment carries the intrinsic incubation period,
typically 4–7 days between infection and becoming infectious. Omitting it would bias the
estimated transmission rate upward, because the model would have to reproduce the
observed epidemic growth without the delay that genuinely slows it.

**Mosquitoes (SI).** There is no recovered class: an infected mosquito remains infectious
until it dies. Mosquito births are set to replace deaths so that the population is
constant in the state variables, with seasonal abundance carried by the time-varying
transmission coefficient instead. Modelling abundance as a second time-varying quantity
would add a parameter that human case data cannot separately identify.

The system is given in full in the module docstring of `src/dengue_pk/models.py`.

## Fixed versus estimated parameters

Human case counts alone cannot identify mosquito lifespan or incubation periods; those
quantities trade off against the transmission rate, and estimating them yields whatever
value happens to fit best rather than anything biologically meaningful. They are
therefore fixed at literature values and their influence is checked by sensitivity
analysis rather than by estimation.

| Parameter | Central set | Alternative set | Status | Basis |
|---|---|---|---|---|
| `gamma_h` human recovery rate | 0.2 /day | 0.25 /day | fixed | viraemic period 4–7 days |
| `sigma_h` 1 / intrinsic incubation | 0.2 /day | 0.167 /day | fixed | 4–7 days to infectiousness |
| `mu_v` mosquito mortality | 0.1 /day | 0.071 /day | fixed | adult lifespan 8–15 days |
| `m` vector-to-host ratio | 2.0 | 3.0 | fixed | poorly constrained; <1 to >10 reported |
| `rho` reporting fraction | 0.05 | 0.05 | **fixed** | confounded with population — see below |
| `beta_0` baseline transmission | — | — | estimated | |
| `a_temp` temperature coefficient/exponent | — | — | estimated | |
| `a_rain` lagged-rainfall coefficient | — | — | estimated | |
| `pop_frac` fraction of the population at risk | — | — | estimated | |
| `i0_frac` initial infected fraction | — | — | estimated | seeding is not observed |

`rho` is **fixed, not estimated**. The state variables are fractions, so the population
enters predictions only where incidence is converted to counts, and `rho` and the
population at risk appear purely as the product `rho x N`. No data can separate them —
verified numerically: halving `N` while doubling `rho` changes every predicted week by
exactly zero. The product is attributed to the population, which is where the larger
error lies, and `rho` is fixed from the under-ascertainment literature.

### Why the fixed parameters became a factor

Both parameter sets above are used, as a factor in the robustness study. Every published
range for these quantities is wide, so choosing within it is an analytical degree of
freedom exactly like the rainfall lag. The study criticises the field for leaving such
choices unexamined; leaving its own unexamined would have made its headline a lower bound
for a reason it could have removed. `FixedParams.from_config(cfg, "fixed_alt")` selects
the alternative set, and neither set sits outside the published ranges.

## Three model structures

The robustness study fits each outbreak under three structures, which differ in what
mechanism they represent rather than in how much freedom they have.

| Structure | Description | Extra estimated parameters |
|---|---|---|
| `hostvector` | Humans SEIR, mosquitoes SI, constant mosquito mortality | — |
| `seir` | Directly transmitted human-only SEIR, no vector | — |
| `hostvector_tempmort` | As `hostvector`, but mosquito mortality follows temperature | **none** |

The third exists to answer the obvious criticism of the first: holding mosquito lifespan
constant across a season whose temperature swings by fifteen degrees is not what
mechanistic dengue modelling does. Mortality follows a quadratic lifespan response with
thermal limits fixed from *Aedes aegypti* trait data (11.7–37.2 °C), normalised so that
at the thermal optimum it equals the constant value the simpler structure uses — the two
models therefore coincide exactly there, which `tests/test_models.py` asserts. Outside
the thermal limits the response is floored so that mortality stays finite.

**What the null model means in the third structure.** Temperature enters *both* compared
models there, because it drives mosquito mortality in the structure itself. The constant
model is therefore not climate-free: it is "climate affects mosquito survival but not the
transmission coefficient". The comparison is the same in all three structures — does
letting the transmission coefficient depend on climate improve the fit, holding the
mechanism fixed — but the baseline is richer in the third, which makes it the hardest of
the three tests. This must be stated whenever the structure factor is reported; a
structure whose null is not the same null would otherwise read as a like-for-like
comparison.

## Climate forcing

    β(t) = β₀ · exp( a_T · z_T(t) + a_P · z_P(t) )

with standardised weekly temperature `z_T` and standardised lagged, smoothed rainfall
`z_P`. The log-linear form keeps β positive for any coefficient values, and
standardising puts the two coefficients on a comparable scale so their relative sizes
can be read directly.

### The rainfall lag, and a correction

Rainfall does not affect reported cases immediately. The full chain is:

| Stage | Duration |
|---|---|
| Rain creates habitat; egg to biting adult | 7–10 days |
| Extrinsic incubation inside the mosquito | 7–14 days |
| Intrinsic incubation in the newly infected person | 4–7 days |
| Symptom onset to reporting | ~7 days |
| **Total** | **≈ 25–38 days, or 4–6 weeks** |

The configuration initially specified a 3-week lag, counting only mosquito development
and forgetting the three delays that follow it. Lagged cross-correlation between weekly
cases and rainfall, computed in `scripts/02_explore_data.py`, returns an optimum of 6
weeks nationally in 2013, 5 weeks in Sindh in 2021 and 4 weeks in Khyber Pakhtunkhwa in
2021 — consistent with the biological chain above and not with the assumed 3 weeks.

The lag is therefore set to **5 weeks**, which is supported both by the mechanism and by
all three empirical estimates. Sensitivity to values between 3 and 7 weeks is reported
rather than assumed to be negligible.

This is recorded here rather than quietly corrected because the original value was a
modelling assumption that the data contradicted, and a reader is entitled to know which
choices were made in advance and which were revised.

### Temperature

The empirical temperature lags — 12 weeks nationally, 4 weeks in both 2021 provincial
windows — are not mutually consistent. The national figure is not interpretable as a
biological delay: that window spans two summers, so a lagged correlation picks up the
annual cycle rather than a causal delay. Temperature therefore enters contemporaneously,
and the national estimate is treated as an artefact of window length. This is stated in
the analysis rather than presented as a finding.

## A structural identifiability result

The reporting fraction and the population at risk **cannot both be estimated**. This is
not a limitation of the data or of the optimiser; it follows from the model's structure.

The state variables are fractions, so the population never enters the dynamics. It
appears only in the final step that converts incidence into counts:

    expected reported cases = ρ · (incidence as a fraction) · N

so ρ and N occur exclusively as the product ρN. Any pair with the same product produces
identical predictions for every week. This was verified numerically rather than left as
algebra: halving N while doubling ρ changes the prediction by **exactly zero**, to the
last bit. `tests/test_inference.py` asserts it.

The first fit made the consequence visible. With ρ free and N fixed at Pakistan's census
population of 188 million, the estimate came back at ρ ≈ 10⁻⁴ — one reported case per ten
thousand infections, two to three orders of magnitude below any published dengue
under-ascertainment estimate. The reporting fraction was absorbing an error in the
population denominator, because it was the only parameter free to do so.

**Resolution.** ρ is fixed at 0.05, roughly one reported case in twenty, from the
under-ascertainment literature, and the population at risk is estimated as `pop_frac ×
census population`. The identifiable product is unchanged; what changes is which quantity
carries it. Attributing it to the population is the better choice here because that is
where the larger error lies: the 2013 epidemic was concentrated in a small number of
districts, and treating 188 million people as one well-mixed population is wrong by a far
greater factor than any plausible error in reporting.

Because the product is what is identified, the estimated population at risk scales
exactly inversely with whatever ρ is assumed. Results are therefore reported across
ρ ∈ {0.02, 0.05, 0.10, 0.20}, and the transmission parameters — which are unaffected —
are the quantities the conclusions rest on.

## Observation model

Surveillance reports cases aggregated over a week, so the model is compared against
weekly differences of cumulative incidence, not against an instantaneous rate:

    expected reported cases in week k = ρ · ( C(t_{k+1}) − C(t_k) ) · N

where `C` is the cumulative human infection incidence and `N` the population.

## Null model

A constant-β model, nested inside the climate-forced model at `a_T = a_P = 0`, is fitted
alongside it. If climate forcing cannot beat the constant model on held-out weeks, the
climate terms are not earning their parameters, and the paper reports that outcome.

## Numerical implementation

The system is integrated with fixed-step classical Runge–Kutta. A fixed step is chosen
over an adaptive solver deliberately: inference calls the integrator tens of thousands of
times, and a deterministic cost per call keeps the timing comparison against the PINN
approach meaningful.

### Step size

One day. Under climate forcing, where the driver itself varies weekly, this agrees with
an eightfold finer step to a relative difference below 10⁻⁶ — verified in
`tests/test_models.py` rather than assumed from the constant-coefficient case, since a
time-varying coefficient is the harder problem.

### Stability, and why the guard cannot bias the results

An explicit method has a finite stability region, so a sufficiently large transmission
coefficient makes the integration diverge. During multi-start fitting the optimiser does
propose such values. The integrator therefore checks the state periodically and raises
`IntegrationFailure` rather than returning `NaN`: a silent `NaN` becomes a meaningless
residual, and the optimiser cannot then distinguish an impossible parameter set from a
merely poor one. The fitting routine converts the failure into a large finite penalty, so
the optimiser sees a flat barrier instead of a cliff of arbitrary values it might try to
descend.

A barrier in parameter space is only harmless if it sits far from where the optimum can
plausibly lie, so its location was measured by bisection rather than assumed:

| Step size (days) | Largest stable β | Equivalent R₀ | Cost per solve |
|---|---|---|---|
| 2.0 | 1.08 | 10.8 | 6 ms |
| **1.0** | **1.79** | **17.9** | **12 ms** |
| 0.5 | 3.29 | 32.9 | 25 ms |
| 0.25 | 6.09 | 60.9 | 66 ms |
| 0.125 | 11.73 | 117.3 | 129 ms |

The boundary scales as β_max ≈ 1.5 / Δt, so halving the step doubles the reachable
transmission rate should a future model need it.

Reported dengue outbreaks have R₀ between roughly 1.5 and 6. The one-day step remains
stable to R₀ ≈ 18, a margin of three to twelve times over anything credible, which also
accommodates the transient excursions climate forcing produces. The barrier therefore
cannot constrain the optimum, and `tests/test_models.py` asserts the margin on every run
rather than leaving it to trust.

### Interpreting β

For this parameterisation the next-generation argument gives

    R₀ = β · sqrt( m / (γ_h · μ_v) )

which with the fixed parameters above is **R₀ = 10 β**. Estimated transmission
coefficients are reported as R₀ throughout, since β alone is not comparable across
studies that parameterise the vector component differently.

## Stated limitations

- Dengue has four serotypes with partial and temporary cross-immunity. A single-strain
  model is an approximation, and the estimated parameters should be read as effective
  values for the season modelled rather than as serotype-specific quantities.
- The national series aggregates regional epidemics with different timing. Quantifying
  the resulting bias, by comparison with the separately fitted 2021 provincial windows,
  is one of the paper's research questions rather than an oversight.
- A single city represents the climate of an entire province or country. Sensitivity to
  that choice is tested on the Pakistani windows but is **not** a factor in the global
  factorial, because a second climate series would have to be downloaded for every one of
  the 237 windows. It remains the largest degree of freedom the study names and does not
  measure, and the headline should be read as a lower bound for that reason.
- Mosquito abundance is carried by β(t) rather than modelled as a separate state. A model
  with explicit larval and adult stages would fork further still.
