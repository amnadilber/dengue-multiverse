# Emails to send once the preprint is live

Rewritten 2026-08-09. Every person below was found by searching the published literature
and every affiliation and address was verified against the paper it appears in. Nobody here
is a guess.

**Send nothing before the preprint has a DOI.** Without a link this is a request for a
favour; with one it is a researcher showing work.

---

## Why this file is now the critical path

The affiliation problem has stopped the work twice: medRxiv declined the submission for it,
and it is the largest single factor against journal acceptance. Every route out of it
reduces to the same requirement — **one person who knows the field agreeing to look at
this.** That is what these emails are for.

Do not ask for co-authorship. Nobody agrees to that in a first email. Ask what is wrong
with the paper. The people who answer that question are the people who might later say yes
to more.

## The rules

1. **Under 200 words.** These people get thirty emails a day.
2. **One specific sentence about their work**, showing you read it. Generic praise is the
   fastest route to deletion.
3. **A link, never an attachment.**
4. **Ask for criticism.** It is the only request that costs the recipient nothing to grant
   and is genuinely useful to you.
5. Send from your own address. A Gmail address is not the problem; nobody looks at it. The
   repository is the credential.
6. **One follow-up only**, after three weeks, one line.

Expect most not to reply. That is normal and it is not about the work.

---

# Tier 1 — Pakistani climate–dengue researchers

These two are the highest-value contacts in this file: the subject is theirs, the data
overlap is real, and a Pakistani independent researcher writing to them is not an odd
event.

## 1. Waqas Shabbir — Institute of Statistics, Alpen-Adria-Universität Klagenfurt

`waqashabbir@outlook.com`
*Shabbir W, Pilz J, Naeem A. "A spatial-temporal study for the spread of dengue depending
on climate factors in Pakistan (2006–2017)." BMC Public Health 20:995 (2020).*

**Why him first.** He is a statistician, he fitted a GLM of Pakistani dengue on rainfall
and maximum temperature, and that is precisely the class of analysis this paper measures.
He is Pakistani, he has an institutional affiliation, and his 2006–2017 window overlaps the
2013 national wave in the appendix.

> **Subject:** Your 2020 BMC Public Health dengue–climate model, and a question about how much of it survives the analyst
>
> Dear Dr Shabbir,
>
> I have posted a preprint that bears directly on your 2020 BMC Public Health paper, and I
> would value your view on where it is wrong.
>
> I fitted every usable dengue outbreak in OpenDengue — 221 in 33 countries — under all 144
> combinations of six analysis choices that papers rarely state: the observation model, the
> temperature functional form, the rainfall lag, how much of the series is fitted, the
> compartmental structure, and the fixed entomological parameters. The verdict on climate
> forcing changes in 92% of outbreaks. The observation model alone flips a third of them.
>
> Your GLM used average maximum temperature and rainfall on Pakistani data covering
> 2006–2017; my appendix analyses the 2013 national wave and the 2021 provincial ones. The
> question I cannot answer from outside is whether the six choices I varied are the ones a
> statistician working on this data would actually consider defensible, or whether I have
> included one nobody would make and missed one everybody does.
>
> Preprint: [DOI]
> Code and result tables: github.com/amnadilber/dengue-multiverse
>
> Any objection, however short, would be useful.
>
> Amna Dilber
> Independent researcher, Pakistan
> ORCID 0009-0008-5684-4516

## 2. Muhammad Nasar-U-Minallah — Institute of Geography, University of the Punjab, Lahore

`Nasarbhalli@gmail.com`
*Rehman W, Nasar-U-Minallah M, Butt I. "Spatial mapping of dengue fever prevalence and its
association with geo-climatic factors in Lahore, Pakistan." Environmental Monitoring and
Assessment 196:812 (2024).*

**Why him.** He works on Lahore dengue and climate, at **the university that awarded your
degree**. That is the affiliation route and the domain expert in the same person. His 2021
Lahore study is the same season as two of the provincial windows in the appendix.

> **Subject:** Dengue and climate in Lahore — a preprint from a Punjab University graduate, and a question
>
> Dear Dr Nasar-U-Minallah,
>
> I took my BS in Mathematics at the University of the Punjab and have since been working
> independently on statistical inference in epidemic models. I have posted a preprint that
> bears on your 2024 study of dengue and geo-climatic factors in Lahore, and I would be
> grateful for your criticism of it.
>
> The question is how much of a climate–dengue conclusion survives the analyst's own
> choices. Fitting 221 outbreaks worldwide under all 144 combinations of six routine
> choices, the verdict changes in 92% of them, and most of that instability is not some
> methods being consistently more permissive: it depends on which method meets which
> epidemic. The appendix analyses Pakistani surveillance, including the 2021 provincial
> waves your study covers for Lahore.
>
> I have no institutional affiliation, which I mention because it is the honest position
> and because it has already cost me a preprint server. What I have is the full pipeline,
> 267 tests, and a log of every formulation that failed.
>
> Preprint: [DOI]
> Code: github.com/amnadilber/dengue-multiverse
>
> If any of it is useful to your group I would be glad to hear.
>
> Amna Dilber
> ORCID 0009-0008-5684-4516

---

# Tier 2 — authors cited in the paper

Lower reply probability, higher value if they answer.

## 3. Rachel Lowe — LSHTM / Barcelona Supercomputing Center

Cited as `lowe2021` for drought and extreme rainfall acting at different delays.

> **Subject:** Preprint: how much of a climate–dengue verdict survives the analyst's choices?
>
> Dear Professor Lowe,
>
> I have posted a preprint that may bear on your work on hydrometeorological drivers of
> dengue in Brazil, and I would value your view on where it is wrong.
>
> Every usable outbreak in OpenDengue — 221 in 33 countries — was fitted under all 144
> combinations of six analysis choices that papers rarely state. Roughly three-quarters of
> the variation in whether climate forcing is endorsed sits between analyses of the same
> outbreak rather than between outbreaks, and partitioned three ways the largest term is
> neither: 61% is an outbreak-by-analysis interaction.
>
> The part I am least sure of is the treatment of the rainfall lag. Your 2021 Lancet
> Planetary Health paper shows drought and extreme rainfall acting at different delays; my
> design treats the lag as a nuisance choice with three levels, which may be too crude to
> be fair to the mechanism.
>
> Preprint: [DOI]. Code: github.com/amnadilber/dengue-multiverse
>
> Amna Dilber
> Independent researcher, Pakistan

## 4. Erin Mordecai — Stanford

Cited as `mordecai` and `mordecai2019` for the *Aedes aegypti* thermal limits.

> **Subject:** Preprint using your Aedes aegypti thermal limits — a question about the exponent
>
> Dear Professor Mordecai,
>
> I have posted a preprint that uses the Brière thermal response with limits from your 2017
> PLOS NTD paper, and I would like to check I have not misused it.
>
> I write temperature into the transmission coefficient as B(T) raised to an estimated
> exponent, so that 0 removes temperature and 1 applies the literature response in full,
> which makes the climate and constant models nested. Across 221 outbreaks this unimodal
> form and the conventional log-linear one disagree about the verdict in 11.5% of cases, and
> the log-linear coefficient is negative in essentially half of all fits — a property of the
> calendar rather than of the biology.
>
> My question is whether the exponent formulation is a defensible way to nest the
> hypothesis, or whether it distorts the thermal response in a way I have not noticed.
>
> Preprint: [DOI]. Code: github.com/amnadilber/dengue-multiverse
>
> Amna Dilber
> Independent researcher, Pakistan

## 5. The OpenDengue team — Oliver Brady / Joseph Clarke, LSHTM

You have something to give them rather than only something to ask, which is the best
position to write from.

> **Subject:** OpenDengue v1.3 across 221 outbreaks — a note on what the weekly subset supports
>
> Dear Dr Brady,
>
> I have posted a preprint built entirely on OpenDengue v1.3, and the selection statistics
> may be of some use to you.
>
> Of the 88 countries in the weekly-resolution subset, 34 yield at least one window that can
> carry a single-epidemic fit under my criteria. Of the 54 that do not, only 17 fail because
> the weekly record has gaps; the other 37 report continuously enough and are excluded
> because no wave is large or sharply peaked enough to fit. I had assumed the binding
> constraint was reporting continuity. It is not — it is outbreak size and shape.
>
> If that breakdown is useful, the script producing it is step 36 in the repository and runs
> on the released extract unmodified.
>
> Preprint: [DOI]. Code: github.com/amnadilber/dengue-multiverse
>
> With thanks for making the compilation available.
>
> Amna Dilber
> Independent researcher, Pakistan

---

# Tier 3 — the metascience community, and arXiv endorsement

The preprint is going to MetaArXiv, whose readership is this group. They are unusually
willing to engage with work from outside a university, because their subject is what goes
wrong inside one.

**Uri Simonsohn, Sara Steegen, Wolf Vanpaemel, Brian Nosek, Andrew Gelman.**

The angle is the one thing in the paper that is about their method rather than about
dengue:

> Specification-curve and many-analyst designs vary the analysis on one dataset, or the
> analyst on one dataset, and so cannot estimate a dataset-by-specification interaction.
> With many datasets each analysed every way it becomes estimable. Here it is 61% of the
> variance, four times the analysis main effect a specification curve displays — which means
> the curve shows the smaller term and averages the larger one away.

**arXiv endorsement.** If you want the paper on arXiv (q-bio.PE or stat.AP) you need an
endorser in that category, and since January 2026 an institutional email no longer
qualifies anyone automatically. Ask only someone who has already replied to you once. An
endorsement request to a stranger is the least likely email in this file to succeed.

---

# When someone replies

Reply within a day. Answer their specific point and ask nothing further. If it runs to two
or three exchanges, then it is reasonable to ask whether they would look at the manuscript
properly — and only after that, whether they would consider joining it.

**A domain co-author joining now would normally**: review the epidemiological framing, help
target the journal, supply the institutional affiliation, and be second author. That is a
good trade. A paper with a co-author at a good journal is worth more than a sole-authored
preprint nobody reviewed.

Do not give away first authorship. You did the work, and `docs/ANALYSIS_LOG.md` records
exactly what was done and when, which is unusually strong evidence if it ever matters.
