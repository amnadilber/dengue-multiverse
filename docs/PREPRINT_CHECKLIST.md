# medRxiv submission — everything the form asks for, filled in

Written 2026-08-03. Copy these fields straight into the submission form at
<https://www.medrxiv.org/submit-a-manuscript>. Screening takes 2–4 days.

---

## Before you start

- [x] Manuscript PDF — `paper/paper.pdf` (25 pages)
- [x] Supplement PDF — `paper/supplementary.pdf` (6 pages, S1–S6)
- [x] Repository public — <https://github.com/amnadilber/dengue-multiverse>
- [x] Repository URL in the manuscript's data-availability statement
- [x] AI disclosure section in the manuscript
- [x] Zenodo DOI — 10.5281/zenodo.21854302, in the manuscript, CITATION.cff and README

---

## Form fields, ready to paste

**Article type:** New Results

**Subject area:** Epidemiology
*(Second choice if asked: Health Informatics)*

**Title**

```
How much of the evidence for climate-driven dengue transmission survives the
analyst? A multiverse analysis of 221 outbreaks in 33 countries
```

**Author**

| Field | Value |
|---|---|
| First name | Amna |
| Last name | Dilber |
| Email | amnadilber.bi@gmail.com |
| Institution | Independent researcher |
| Corresponding author | Yes |
| ORCID | `0009-0008-5684-4516` |

> "Independent researcher" is accepted. medRxiv does not require an institutional address,
> and it does not verify affiliation. Do not invent one.

**Abstract** — paste from `paper/paper.tex` (294 words, under the 300-word cap).

---

## The declarations — these are what screening actually checks

**Competing interest statement**

```
The author declares no competing interests.
```

**Funding statement**

```
This research received no specific grant from any funding agency in the public,
commercial, or not-for-profit sectors.
```

**Author declarations / ethics** — this is the field most often filled in wrongly.

```
All data analysed in this study are publicly available, aggregated, and
de-identified: weekly reported dengue case counts from OpenDengue v1.3
(CC BY 4.0) and meteorological reanalysis from NASA POWER (public domain). No
human subjects were involved, no individual-level data were accessed, and no
ethical approval or consent was required.
```

**Data availability statement**

```
All analysis code, result tables and figures are available at
https://github.com/amnadilber/dengue-multiverse and archived at
https://doi.org/10.5281/zenodo.21854302. Raw data are not redistributed;
download scripts record source URLs and checksums. Case data are from
OpenDengue v1.3 (CC BY 4.0); climate data from NASA POWER (public domain).
```

**Clinical trial?** No.
**Prospective registration?** Not applicable.

---

## Zenodo DOI — do this first, it takes ten minutes

A DOI is a permanent address. The GitHub link can change; a DOI cannot. Journals increasingly
ask for one, and a reviewer who clicks a dead link stops trusting the rest.

1. Go to <https://zenodo.org> → **Log in with GitHub** (use Amna's GitHub account)
2. Zenodo asks for permission to see your repositories → allow
3. Go to **Zenodo → your profile → GitHub**
4. Find `amnadilber/dengue-multiverse` in the list and switch the toggle **ON**
   *(nothing happens yet — Zenodo now watches for releases)*
5. Back on GitHub: **Releases → Create a new release**
   - Tag: `v1.0.0`
   - Title: `v1.0.0 — preprint submission`
   - Description:
     ```
     Code, result tables and figures accompanying the manuscript submitted as a
     preprint. 38 numbered pipeline steps, 266 tests, and a dated analysis log.
     ```
   - **Publish release**
6. Within a few minutes Zenodo mints a DOI and shows it on your Zenodo dashboard
7. Put that DOI in three places:
   - `CITATION.cff` — add a `doi:` field
   - the manuscript's data-availability statement
   - the medRxiv data-availability field above

**Order matters:** Zenodo DOI → update manuscript → recompile PDF → submit to medRxiv. Doing
it the other way round means the preprint points at nothing.

---

## After the preprint is live

You get a DOI like `10.1101/2026.MM.DD.XXXXXXXX`. Then:

- [ ] Add it to `CITATION.cff` (`preferred-citation`) and to the README
- [ ] Add "Preprints" to your CV:
      *Dilber, A. (2026). How much of the evidence for climate-driven dengue
      transmission survives the analyst? medRxiv. DOI: …*
- [ ] Email it to two or three authors you cite — Rachel Lowe, Erin Mordecai, the
      OpenDengue team. Short message: what the paper does, why it bears on their work,
      an invitation to point out what is wrong. Most will not reply. **One might, and
      one is what you need** — a domain co-author is the single change that would most
      raise the chance of journal acceptance.
- [ ] Only then submit to a journal (see `SUBMISSION.md` for order and costs)

A preprint does **not** block journal submission. PLOS, Elsevier and the Royal Society all
accept previously preprinted work; the submission form asks you to declare it, and you
declare it.

---

## What screening rejects, so you can avoid it

medRxiv rejects for: incomplete author or affiliation details; missing competing-interest or
funding statements; material already published elsewhere; content that could pose a public
health or biosecurity risk; article types that are not research (opinion, review without
methods and data).

None of these applies here — the manuscript has methods, data and results, and the
declarations above are complete. The one field that trips people up is the ethics
declaration, which is why the exact wording is given above rather than left to be written
under time pressure.
