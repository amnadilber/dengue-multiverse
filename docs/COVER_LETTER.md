# Cover letter

The editor reads this before the abstract, and often decides on it alone. Three things
have to be on the first screen: **what is new**, **why it belongs in this journal**, and
**that it is not only a complaint**.

Keep it under one page. Do not summarise the paper — the abstract does that.

---

## Version for PLOS Global Public Health

> Dear Editors,
>
> I am submitting *"How much of the evidence for climate-driven dengue transmission
> survives the analyst? A multiverse analysis of 221 outbreaks in 33 countries"* for
> consideration at PLOS Global Public Health.
>
> Mechanistic models fitted to surveillance data are routinely used to test whether dengue
> transmission is climate-driven, and reviews explain their conflicting findings by
> appealing to local context. Whether that explanation is right has not been checkable,
> because each published estimate comes from one analysis of one place and the two cannot
> be separated after the fact. This study separates them: every usable outbreak in a global
> surveillance compilation is fitted under all 144 combinations of six analysis choices
> that published models rarely state — 33,152 fits.
>
> Roughly three-quarters of the variation in whether climate forcing is endorsed lies
> between analyses of the same outbreak rather than between outbreaks. Partitioned three
> ways, the largest term is neither: 61% is an outbreak-by-analysis interaction against 16%
> for the analysis itself. That distinction decides which remedy is available — a main
> effect can be settled by convention and an interaction cannot — and it is invisible to a
> specification curve, which displays the main effect and averages the interaction away.
>
> The paper does not stop at the diagnosis. On simulated data whose answer is known,
> conventional 95% intervals contain the truth 71–78% of the time and 42% when the epidemic
> is a sum of two asynchronous ones, which these data demonstrably are. Combining **eight**
> analyses by Rubin's rules restores nominal coverage at 1.5–2.2 times the width — a few
> minutes of computation, not a factorial. This is the only recommendation in the paper
> validated against a known truth rather than argued for, and it is the one we would keep if
> only one survived.
>
> Applying that rule to the data gives a clear verdict for 137 of 221 outbreaks, 86% of
> which support climate forcing. The claim is not that the effect is absent. It is that the
> confidence with which the question is currently answered is not earned, and that in 38% of
> outbreaks a single season cannot settle it.
>
> The work is relevant to PLOS Global Public Health because the surveillance data on which
> this literature rests comes overwhelmingly from countries where dengue is endemic, and the
> recommendation is deliberately cheap enough to be adopted without specialist computing.
>
> All analyses are reproducible from a public repository containing 38 numbered pipeline
> steps, 268 tests — including automated checks that the manuscript's quoted figures still
> match the stored result tables — and a dated analysis log recording every formulation that
> failed and every error found, including those that produced plausible wrong numbers rather
> than crashing.
>
> This manuscript is not under consideration elsewhere. A preprint is posted on medRxiv
> [DOI]. I have no competing interests and received no funding. I am an independent
> researcher without institutional affiliation, and I would be glad of reviewers willing to
> be blunt.
>
> Yours sincerely,
> Amna Dilber

---

## Changes for other venues

**Epidemics** — replace the penultimate-but-one paragraph with:

> The work is relevant to *Epidemics* because it is a study of inference in epidemic models
> rather than a forecasting paper: the object of study is the estimator and the decision
> rule, not the epidemic. The main text is accompanied by a supplement holding the factor-
> list robustness analysis, the criterion comparison, further checks and a worked
> identifiability case; I am happy to move more material there if the editors prefer a
> shorter main text.

**PLOS Neglected Tropical Diseases** — add after the second paragraph:

> The study takes as its starting point the systematic review of 99 dengue prediction models
> published in this journal (Leung et al., PLOS NTD 2023), which found that "the reporting
> of methodology and model performance measures were inadequate in many of the existing
> prediction models". This paper measures the consequence of that inadequacy.

**Royal Society Open Science** — lead with the methodological point instead:

> The contribution is partly to multiverse and specification-curve methodology itself. Those
> designs vary the analysis on one dataset, or the analyst on one dataset, and so cannot
> estimate a dataset-by-specification interaction. With many datasets each analysed every
> way it becomes estimable, and here it is four times the term a specification curve
> displays.

---

## What not to write

- Do not apologise for being unaffiliated. State it once, plainly, at the end.
- Do not claim the paper is important. Say what it measured and let the number speak.
- Do not oversell "first ever". "To our knowledge" is in the paper; that is enough.
- Do not list the limitations here. They are in the paper, which is where an editor
  expects to find them.
