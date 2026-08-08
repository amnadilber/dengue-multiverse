# Emails to send once the preprint is live

**Send these only after the medRxiv DOI exists.** Without a link, the email is a request for
a favour. With one, it is a colleague showing work.

## The rules that make these work

1. **Ask for criticism, not co-authorship.** Nobody agrees to co-author an email. People do
   answer "what did I get wrong?"
2. **Under 200 words.** These people get thirty emails a day.
3. **Say what it bears on *their* work** in one sentence, specifically. Generic praise is
   the fastest route to being ignored.
4. **No attachment.** A link.
5. **Send from your own address**, subject line factual, no "URGENT" or "Request".
6. **Do not follow up more than once**, after three weeks, in one line.

Expect eight of ten not to reply. That is normal and is not about the work.

---

## 1. Rachel Lowe (LSHTM) — climate-dengue modelling

**Subject:** Preprint: how much of a climate–dengue verdict survives the analyst's choices?

> Dear Professor Lowe,
>
> I have posted a preprint that may bear on your work on hydrometeorological drivers of
> dengue in Brazil, and I would value your view on where it is wrong.
>
> I fitted every usable outbreak in OpenDengue — 221 in 33 countries — under all 144
> combinations of six analysis choices that papers rarely state: the observation model, the
> temperature functional form, the rainfall lag, how much of the series is fitted, the
> compartmental structure, and the fixed entomological parameters. Roughly three-quarters of
> the variation in whether climate forcing is endorsed sits between analyses of the same
> outbreak rather than between outbreaks, and most of that is an outbreak-by-analysis
> interaction rather than some methods being systematically more permissive.
>
> The part I am least sure of is the treatment of the rainfall lag. Your 2021 Lancet
> Planetary Health paper shows drought and extreme rainfall acting at different delays; my
> design treats the lag as a nuisance choice with three levels, which may be too crude to be
> fair to the mechanism.
>
> Preprint: [medRxiv DOI]. Code and data: github.com/amnadilber/dengue-multiverse
>
> I would be grateful for any objection, however brief.
>
> Amna Dilber
> Independent researcher, Lahore

---

## 2. Erin Mordecai (Stanford) — thermal biology of transmission

**Subject:** Preprint using your *Aedes aegypti* thermal limits — a question about the exponent

> Dear Professor Mordecai,
>
> I have posted a preprint that uses the Brière thermal response with limits from your 2017
> PLOS NTD paper, and I would like to check I have not misused it.
>
> I write temperature into the transmission coefficient as B(T) raised to an estimated
> exponent, so that 0 removes temperature and 1 applies the literature response in full,
> which makes the climate and constant models nested. Across 221 outbreaks I find that this
> unimodal form and the conventional log-linear one disagree about the verdict in 11.5% of
> cases — and that the log-linear coefficient is negative in essentially half of all fits,
> which is a property of the calendar rather than of the biology.
>
> My question is whether the exponent formulation is a defensible way to nest the hypothesis,
> or whether it distorts the thermal response in a way I have not noticed.
>
> Preprint: [medRxiv DOI]. Code: github.com/amnadilber/dengue-multiverse
>
> Any correction would be welcome.
>
> Amna Dilber
> Independent researcher, Lahore

---

## 3. The OpenDengue team (Oliver Brady / Joseph Clarke, LSHTM) — the data

**Subject:** OpenDengue v1.3 used across 221 outbreaks — a note on what the weekly subset supports

> Dear Dr Brady,
>
> I have posted a preprint built entirely on OpenDengue v1.3, and I thought the selection
> statistics might be of some use to you.
>
> Of the 88 countries in the weekly-resolution subset, 34 yield at least one window that can
> carry a single-epidemic fit under my criteria. Of the 54 that do not, only 17 fail because
> the weekly record has gaps; the other 37 report continuously enough and are excluded
> because no wave is large or sharply peaked enough to fit. I had assumed the binding
> constraint was reporting continuity and it is not — it is outbreak size and shape.
>
> If that breakdown is useful, the script producing it is step 36 in the repository and it
> runs on the released extract unmodified.
>
> Preprint: [medRxiv DOI]. Code: github.com/amnadilber/dengue-multiverse
>
> With thanks for making the compilation available.
>
> Amna Dilber
> Independent researcher, Lahore

---

## 4. Xing Yu Leung / Md Nazmul Karim — authors of the 99-model review

**Subject:** Following up your 2023 review — measuring the consequence of the reporting gap

> Dear Dr Leung,
>
> Your 2023 PLOS NTD review found that "the reporting of methodology and model performance
> measures were inadequate in many of the existing prediction models". I have posted a
> preprint that tries to measure what that inadequacy costs.
>
> Fitting 221 outbreaks under all 144 combinations of six unreported choices, the verdict on
> climate forcing changes in 92% of them; between the most and least credulous specification
> lies 56 percentage points of endorsement rate. Two of your review's own figures do a lot of
> work in the paper — the 18.3% using Poisson regression and 18.3% linear regression, and the
> 20.2% reporting no validation at all — as evidence that the factor list describes practice
> rather than my imagination.
>
> I would value your view on whether the six choices I varied are the right ones, and what a
> reader of your review would say I have missed.
>
> Preprint: [medRxiv DOI]. Code: github.com/amnadilber/dengue-multiverse
>
> Amna Dilber
> Independent researcher, Lahore

---

## Two more worth sending

**Uri Simonsohn / Sara Steegen** — the multiverse methodologists. Angle: their designs vary
the analysis on one dataset and so cannot estimate a dataset-by-specification interaction;
with many datasets each analysed every way it becomes estimable, and here it is four times
the term a specification curve displays. This is the audience most likely to find the
methods contribution interesting, and least likely to care about your affiliation.

**Anyone who replies at all** — reply within a day, answer the specific point, and do not
ask for anything. If a conversation continues over two or three exchanges, then it is
reasonable to ask whether they would be willing to look at the manuscript properly. Not
before.

---

## If someone offers to collaborate

Say yes, and be clear about what you did and what you want. A domain co-author who joins at
this stage would normally: review the epidemiological framing, help target the journal,
and be second author. That is a good trade — a paper with a co-author at a good journal is
worth more than a sole-authored preprint nobody reviewed.

Do not give away first authorship. You did the work; the analysis log records exactly what
was done and when, which is unusually strong evidence if it ever matters.
