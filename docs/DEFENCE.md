# What a reviewer will ask, and the answer

Written 2026-08-03. Every answer here is checkable against a numbered script or a stored
table — no answer is "I think so".

The point of this file is not to have arguments ready. It is that **you should be able to
answer these without looking**, because a reviewer's follow-up question is where a paper is
actually judged. If an answer below surprises you, that is the part to go and understand.

Where an answer admits a weakness, it admits it. A reviewer who catches you defending an
indefensible point stops believing the defensible ones.

---

## A. The headline

**Q1. Isn't 144 analyses a straw man? No real analyst would try all of them.**

Correct, and the paper does not claim otherwise. The claim is that a real analyst picks
*one* cell of this grid, usually without saying which, and that the cell they pick decides
the answer. The grid is a way of measuring how much the pick matters — not a description of
what anyone does.

**Q2. You have manufactured instability by running enough analyses that something is bound
to disagree.**

This is the most likely first objection and the paper answers it three ways.

1. The headline statistic is *pairwise disagreement*, `2p(1-p)`, which does not grow with
   the design. It reads 29.9%, 31.4% and 29.8% on designs of 24, 48 and 144 cells
   (Section 4.3 table). The share-of-outbreaks statistic does grow, and the paper says so
   explicitly rather than quoting it alone.
2. The simulation settles it directly. Same windows, same 144 analyses, but a known truth:
   disagreement is 44.5% when no effect exists and 18.6% when a real one does. If the
   design manufactured disagreement, the two arms would not differ.
3. The variance decomposition is a share, not a count. Adding analyses adds to numerator
   and denominator alike.

**Q3. Three-quarters of variation being "within outbreak" is not the same as it being the
analyst's fault.**

Agreed — and this exact objection is why the paper reports the three-way partition
(Section 4.1). The within-outbreak remainder splits into an analysis main effect (15.9%)
and an outbreak-by-analysis interaction (60.6%). An earlier draft glossed the two-way split
as "the analyst", which was wrong, and the paper says so.

**Q4. Your interaction term is confounded with residual error — you have one observation
per cell.**

True by design, and the paper states it. The bound is empirical: each cell is a single
deterministic fit, so the only noise available is the optimiser's, and the cold-start check
(step 24) changes 0.6% of verdicts. That bounds the noise contribution at roughly one part
in a hundred of a term worth three-fifths.

**Q5. The decomposition depends on which factors you chose. Add more and the analytical
share rises.**

Yes, and the paper measures how much rather than defending against it (Section 4.1, second
table). Dropping any single factor except the observation model costs at most six points.
Dropping the observation model costs twenty-one and leaves 55.2%. Across all 63 non-empty
subsets the median is 54.5% and the adversarial minimum — keeping only the fixed
parameters — is **2.3%**. That number is in the paper because a sceptic would find it
otherwise.

---

## B. The model

**Q6. Why a host-vector SEIR rather than [some other structure]?**

Because dengue is vector-borne and a single-population model omits the extrinsic incubation
delay. But the choice is not load-bearing: model structure is one of the six factors, and
it is among the weakest (9.9% of verdicts flipped for host-vector → SEIR, 9.4% for
temperature-dependent mosquito mortality). The paper's result does not depend on the model
being right, which is the point of including structure as a factor at all.

**Q7. Why the Brière form for temperature?**

It is the standard unimodal thermal-performance form for arthropod traits (Brière et al.
1999), with limits taken from *Aedes aegypti* trait data (Mordecai et al. 2017, 2019). It
is written as an exponent, `B(T)^a_temp`, so that `a_temp = 0` removes temperature and
`a_temp = 1` applies the literature response in full — which makes the comparison nested.
The alternative log-linear form is the other level of that factor.

**Q8. Your fixed parameters are guesses.**

They are literature values, and the point is that every published range is wide. That is
exactly why the *choice within* the range is one of the six factors. Moving all of them to
the other end of their published ranges flips **3.1%** of verdicts and moves the
endorsement probability by 0.001 — the weakest factor of the six. We expected it to matter
and it does not.

**Q9. ρ (reporting fraction) and N (population at risk) are not separately identifiable.**

Correct, and the paper says so and demonstrates it: halving N while doubling ρ changes every
predicted week by exactly zero. ρ is therefore fixed from the under-ascertainment
literature and N estimated. Varying ρ over a tenfold range rescales N as 1/ρ and leaves R₀
unchanged to four decimal places.

**Q10. Fixed-step RK4 at one day — is that adequate under climate forcing?**

Verified against an eightfold finer step (`tests/test_models.py`). β and μ_v are evaluated
on a half-step grid so the forcing is not aliased.

---

## C. The statistics

**Q11. AIC is the wrong criterion. Use BIC / a likelihood-ratio test.**

Tested, and both are **worse** (Section 4.7). All three are thresholds on the same
quantity: AIC sign is ΔAIC < 0, the LRT at p<0.05 is ΔAIC < −1.99, BIC at the median 38
weeks is ΔAIC < −3.28. Instability rises monotonically with strictness: 91.9%, 95.0%,
96.8%. BIC with a ±4 band still leaves 90.5% unstable, because its extra penalty shifts the
abstention band off centre — it abstains over (−7.3, +0.7), so one side is still allowed to
speak. **The property that matters is symmetric abstention about zero, not strictness.**

**Q12. Including Poisson inflates your result — nobody should use it on overdispersed
counts.**

Two answers. (a) The field does: the Leung et al. review of 99 dengue prediction models
finds 18.3% using Poisson regression and 18.3% linear regression, with the negative binomial
not appearing as a category. That is a large minority, not a majority, and the paper says
"large minority" rather than "most" — an earlier draft said "most" and was corrected.
(b) The result survives without it: restricting the entire study to negative-binomial fits
still leaves **79.6%** of outbreaks changing verdict (95% CI 75–85).

**Q13. Why does the observation model matter so much?**

Not for the reason usually given. The usual account is that Poisson makes the fit chase
noise via the climate covariates — that predicts larger climate coefficients under Poisson,
and **it is false**: paired within outbreak and within every other choice, median |a_temp|
is 0.092 under NB and 0.097 under Poisson, larger under Poisson in only 40% of 15,912 pairs.
The correct mechanism is that the climate model nests the constant one, so it always fits
slightly better, and without a dispersion parameter to absorb the residual the same
improvement buys far more log-likelihood: median ΔAIC is −1.6 under NB and **−68.9** under
Poisson. Poisson does not change what the model does; it changes what the criterion pays
for it, and it pays in one direction only (climate wins under Poisson alone in 33.6% of
pairs, under NB alone in 0.6%).

**Q14. Is this just optimiser noise?**

Bounded, not assumed. Seven outbreaks were refitted with ten cold restarts and no warm
start. The thorough setting does find better optima — it improves the climate fit in 6.1% of
comparisons, once by 1,656 AIC units — but changes the verdict in **0.6%** of 329 paired
comparisons. The check is small and is reported as a bound, not an estimate.

**Q15. Your multiverse cells are warm-started from a common anchor, so they are more alike
than independent analyses.**

True, and an earlier draft claimed this meant the reported instability was an underestimate.
That direction was never tested; when it was, cold restarts gave slightly *less*
disagreement (32.9% vs 33.4% on seven outbreaks). The paper now says the sign is not known.
Note the check varies the warm start but not the shared dispersion estimate, so one of the
two mechanisms remains untested — also stated.

---

## D. The data

**Q16. 33 countries but three supply 45% of the sample. Is this a fact about Latin America?**

Checked several ways. Removing each well-represented country in turn moves the
within-outbreak share between 76.0% and 78.6%. Removing Bolivia, Nicaragua and Mexico all
at once (leaving 121 outbreaks) gives **78.5%**. National vs subnational: 79.1% vs 75.4%.
The bootstrap CI resamples whole countries rather than outbreaks, because outbreaks within a
country are not independent.

**Q17. Your window-selection rule is itself a researcher degree of freedom.**

Yes, and it is checked (Section 4.11). Recomputing the headline on progressively stricter
peak-prominence subsets gives 91.9%, 91.8%, 94.1% (n=68). Tightening our own selection does
not weaken the result.

**Q18. What does requiring an unbroken weekly run exclude?**

Measured in step 36, and the answer corrected an assumption. Of 88 countries in the weekly
subset, 34 contribute a usable window. Of the 54 excluded, only **17** fail the gap-free run
requirement; **37** report continuously enough and are excluded because no wave is large or
sharply peaked enough to fit. So the dominant selection is on **outbreak size and shape**,
not reporting continuity — and that points the opposite way from the usual concern: large,
single-wave epidemics are what these models handle best, so this instability is what
survives on the most favourable data available.

**Q19. One climate grid cell per country or province is crude.**

Quantified rather than conceded. Forty outbreaks refitted under all 144 combinations at a
second point one degree of latitude away: the verdict changes in **12.0%** of 5,664 paired
comparisons, endorsement probability essentially unmoved (0.739 → 0.737). Read as an upper
bound — the offset is mechanical and respects no terrain, so for Costa Rica the two series
differ by 4.3 °C because the shift crosses high ground.

---

## E. The remedy

**Q20. Rubin's rules are for multiple imputation. What justifies using them here?**

The structure is identical: several analyses of one dataset, each with its own estimate and
its own within-analysis variance. T = W + (1+1/m)B. The justification offered in the paper
is not analogy but **validation** — coverage is measured against a known truth, and it is
the only recommendation in the paper validated that way rather than argued for.

**Q21. Why eight analyses?**

Measured, not chosen. Drawing m analyses at random from an outbreak's 72 and combining only
those: coverage of the rainfall coefficient is 91.0% at m=2, 96.2% at 4, **97.3% at 8**, and
97.5% at 72. The interval also stops narrowing after eight (median half-width 0.545 against
0.542 at 72). Four already reach nominal; eight is indistinguishable from seventy-two.

**Q22. Which eight?**

Also measured, and the answer was **not** what we assumed. We expected the observation
model to be the choice that must be varied, since it is the largest lever on the verdict. It
is close to irrelevant to the interval. What cannot be held fixed is the **rainfall lag**,
and only for the **rainfall** coefficient, which loses nine points of coverage. The rule
that fits: *vary the choices that construct the covariate whose coefficient you are
reporting.* Independently, the lag is also the factor whose effect is least predictable in
advance (0.1% of its variation explained). Neither result was designed to test the other.

**Q23. Your coverage study is circular — you simulate from the family you then fit.**

Answered by repeating it with a generator outside the family, and the mis-specification is
not arbitrary: the appendix demonstrates it in these data. Counts were generated from two
climate-forced patches seeded 46 days apart (the offset the two-patch model infers for
Khyber Pakhtunkhwa) and fitted with one patch. The conventional interval degrades to
**42%**; the multiverse interval degrades too, from 97.5% to **90.0%**. The paper reports
the 90% as the number to carry, not the 97.5%.

**Q24. A margin of 4 is arbitrary.**

It is Burnham and Anderson's long-standing guidance, not ours, and the threshold profile
(step 20) shows the collapse happens between 2 and 4 rather than at a tuned point. The
contribution is not the rule but the **measured cost of not following it**: pairwise
disagreement 29.8% → 4.0% at no cost in outbreaks answered.

---

## F. The uncomfortable ones

**Q25. So is climate driving dengue transmission or not?**

Applying the paper's own rule to its own data: a clear verdict for **137 of 221** outbreaks,
**86% of which support climate forcing**. The paper is not a denial. It says the confidence
with which the question is currently answered is not earned, and that in 38% of outbreaks a
single season cannot answer it by this method.

**Q26. Does your real data look like a literature detecting nothing?**

No, and this is measurable because the simulation has a known truth. The analysis main
effect is 47.0% in the no-effect arm, 10.3% in the real-effect arm, and **15.9%** in the
real data — with the effect arm. Reported as evidence and not proof: the simulated effect
size comes from these same fits and the arms are not a calibrated mixture. Two earlier
one-number attempts at this comparison pointed in opposite directions and one was withdrawn;
the paper says so.

**Q27. This is a critique. What is the positive contribution?**

Three. (a) The three-way partition, which shows the disagreement is neither geographic nor a
ranking of methods but specific to the pairing — and therefore that no house style fixes it.
(b) A validated interval that restores nominal coverage for eight fits. (c) The paper's own
answer under its own rule. The first is also a statement about specification-curve
methodology generally: the curve displays the analysis main effect and averages the
interaction away, and here the hidden term is four times the displayed one.

**Q28. You are a single author with no institutional affiliation and no domain co-author.**

True, and there is no answer that makes it untrue. What is offered instead is verifiability:
38 numbered pipeline steps, 266 tests including automated checks that the manuscript's
quoted numbers still match the stored result tables, and a dated analysis log recording
every formulation that failed and every bug found — including two that produced plausible
wrong numbers rather than crashing. The work does not ask to be taken on authority.

**Q29. Has anyone outside the project read this?**

Not yet, at the time of first submission. Say so if asked. It is a real limitation and
pretending otherwise is worse than admitting it.

**Q30. Did you use AI to do this?**

Yes, and it is disclosed in the manuscript in the form journals require, immediately before
the appendix. A large language model assisted with implementing the code, running the
pipeline, and drafting and revising the text. The research question, study design, the
decisions about which results to pursue and which to discard, and the decision to subject
the manuscript to repeated adversarial review were the author's. The numbers are checked
against the pipeline by the test suite rather than by assertion. The author takes
responsibility for the content.

---

## The three questions to be ready for above all

If you prepare only three, prepare **Q2** (manufactured instability), **Q13** (why the
observation model dominates) and **Q23** (circularity of the coverage study). Those are the
ones that decide whether a reviewer believes the paper.
