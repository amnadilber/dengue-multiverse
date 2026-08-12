# Posting the preprint

Rewritten 2026-08-09 after medRxiv declined the submission. The earlier version of this
file asserted that "Independent researcher is accepted; medRxiv does not require an
institutional address". **That was wrong and was never verified.** It cost four days.

---

## Update 2026-08-12: MetaArXiv also declined it

Moderator feedback, in full: *"Outside of the scope of this preprint series."* A
dengue-titled manuscript on a metascience server was triaged by its title, which in
retrospect is exactly how a volunteer moderator has to work. The methodological framing
that made the fit plausible sits in the abstract's third sentence; moderation decisions
are made on the first.

Two venues, two rejections, neither about the manuscript: one policy (affiliation), one
scope. The fallback that was written into SUBMISSION.md on day one is now the plan:
**Zenodo record for the manuscript (no moderation, instant DOI) + direct journal
submission.** See ZENODO_PREPRINT.md for the paste-ready record.

---

## What happened, so it is not repeated

Submitted to medRxiv on 8 Aug 2026 (MEDRXIV/2026/360009). Declined on 9 Aug:

> *"medRxiv requires authors to have an organizational affiliation. It is necessary for
> submissions to be associated with an organization that provides oversight of research
> activities so that it can adjudicate any ethical issues/disputes that arise."*

Nothing about the manuscript was at issue. The submission was complete and every
declaration was accepted. The affiliation was the sole reason.

**Rule learned: verify a venue's eligibility criteria before preparing a submission, not
after.** The paper was ready four days earlier than it needed to be.

---

## Where it went second (declined 12 Aug): MetaArXiv

<https://osf.io/preprints/metaarxiv/>

Verified 2026-08-09 against the Center for Open Science's own announcement: the OSF
*generalist* server was suspended in August 2025 and will **not** return, but fourteen
community-run servers remain fully operational, MetaArXiv among them.

**Why this one, and why it may be the better home anyway.** The paper's most novel claim is
not about dengue. It is that a specification curve displays the analysis main effect and
averages away the dataset-by-analysis interaction, which here is four times larger. That is
metascience: research about how research is done. MetaArXiv is run by the Berkeley
Initiative for Transparency in the Social Sciences and its readership is the multiverse and
reproducibility community, which is precisely the audience for that claim.

What it costs: epidemiologists will not browse it. That is a real loss, and the answer to it
is the outreach emails, not the server.

**No institutional affiliation is required** (OSF preprint FAQ: affiliation is optional and
attaches only through an email domain if you have one). Submissions are moderated before
posting. A DOI and a persistent URL are assigned on acceptance.

---

## Form fields, ready to paste

**Title**

```
How much of the evidence for climate-driven dengue transmission survives the analyst? A multiverse analysis of 221 outbreaks in 33 countries
```

**Abstract** — paste from `paper/abstract_for_form.txt` (287 words, ASCII only, no special
characters that a web form can mangle). Do **not** copy from the PDF; the line breaks come
with it.

**File** — `paper/paper.pdf` (24 pages).
Add `paper/supplementary.pdf` as a second file if the form allows supplementary material;
if it does not, the supplement is on GitHub and Zenodo and the data-availability statement
points there.

**Contributors**

| Field | Value |
|---|---|
| First name | Amna |
| Last name | Dilber |
| Email | amnadilber.bi@gmail.com |
| ORCID | 0009-0008-5684-4516 |
| Affiliation | leave blank, or `Independent researcher` if the field is free text |

**Licence** — `CC-By Attribution 4.0 International`. Matches the code (MIT) and the source
data (OpenDengue CC BY 4.0). Not CC-0: attribution is the one thing worth keeping.

**Disciplines / subjects** — the form requires at least one top-level discipline.

- Primary: **Social and Behavioral Sciences → Methodology** (or the nearest equivalent
  offered; MetaArXiv's taxonomy is social-science shaped)
- Add: **Life Sciences → Epidemiology** and **Physical Sciences and Mathematics →
  Statistics and Probability** if multiple subjects are permitted

**Author assertions** — these increase the chance a moderator accepts it quickly, and every
one is true here.

| Assertion | Answer |
|---|---|
| Public data available? | **Yes** — `https://github.com/amnadilber/dengue-multiverse` and `https://doi.org/10.5281/zenodo.21854302` |
| Preregistration? | **No** |
| Conflicts of interest? | **None** |

**Supplementary / data links**

```
https://github.com/amnadilber/dengue-multiverse
https://doi.org/10.5281/zenodo.21854302
https://opendengue.org
https://power.larc.nasa.gov
```

**Conflict of interest statement**

```
The author declares no competing interests. This research received no funding.
```

---

## After it is accepted

You get a DOI and a persistent OSF URL. Then, in this order:

- [ ] Add the DOI to `CITATION.cff` (`preferred-citation`), to `README.md`, and to the CV
      under **Preprints**:
      *Dilber, A. (2026). How much of the evidence for climate-driven dengue transmission
      survives the analyst? MetaArXiv. https://doi.org/…*
- [ ] Send the outreach emails in `OUTREACH_EMAILS.md`. Every one now carries two asks:
      a criticism, and — where the recipient publishes on arXiv — an **endorsement** for
      q-bio.PE or stat.AP.
- [ ] Only then approach a journal (`SUBMISSION.md`).

---

## The affiliation problem, stated plainly

It has now blocked the work twice: medRxiv declined it, and it is the largest single factor
against journal acceptance. It is no longer a background weakness. It is the critical path.

Options, with what each actually costs:

| Route | Cost | Realistic? |
|---|---|---|
| **A domain co-author with an affiliation** | one person saying yes | **the real answer** — solves the preprint server, the journal, and the credibility question at once |
| Former department at University of the Punjab | an email; possibly a visiting/affiliate status | plausible, and the degree is from there |
| A Pakistani dengue researcher | an email | plausible; dengue is endemic there and the appendix is Pakistani data |
| Current employer (Synaptrix) | their **written** permission | only with permission. They do not oversee this research, and some employers claim IP over employees' work. Do not list them silently. |
| Ronin Institute | application windows, attend an event; without a doctorate, a published track record | slow |
| IGDORE | applications suspended since Aug 2024 | closed |
| arXiv | an endorser in the category; since Jan 2026 an institutional email no longer qualifies anyone | needs the same person as the co-author route |

Every row except the last two comes down to the same thing: **one person who knows the
field agreeing to look at this.** That is what the emails are for.
