# Where to send this, and in what order

Written 2026-07-27. Verify every fee and policy before acting — APCs and waiver
schemes change yearly, and this file will go stale.

## The constraint that actually decides things

Two things are true at once and they pull in opposite directions:

* **Peer review takes 4–10 months.** A paper submitted in September 2026 is
  unlikely to be accepted before mid-2027. Scholarship applications for the 2027
  intake close between December 2026 and February 2027.
* **A citable preprint takes days.** It carries a DOI, it can be listed on a CV
  under "Preprints", and it can be read by a professor evaluating an application.

So the order is not "pick a journal". It is **preprint first, journal in
parallel**. Waiting for acceptance means the work contributes nothing to the
2027 cycle, which is the reason it exists.

## Step 1 — preprint

| Server | Fee | Fit | Catch |
|---|---|---|---|
| **medRxiv** | free | epidemiology, health methodology | Screening; requires a plausible affiliation and a completed competing-interests statement. Independent researchers are accepted. |
| **arXiv** (q-bio.PE) | free | quantitative biology, populations and evolution | **First-time submitters to q-bio need endorsement** from an existing arXiv author in that category. This is a real hurdle without an institution. |
| **Zenodo** | free | anything | No screening, instant DOI, but carries no venue signal at all. Use as a fallback and for the code archive regardless. |

**Recommended:** medRxiv for the manuscript, Zenodo for the code and data-derived
result tables (Zenodo mints a DOI per GitHub release, so the repository becomes
citable at a fixed version). If medRxiv screening rejects on affiliation grounds,
fall back to Zenodo and cite that.

## Step 2 — journal

Ordered by what this paper actually is: a methodological robustness study of
epidemic-model inference, not a dengue-forecasting paper.

| Journal | Cost to publish | Fit | Notes |
|---|---|---|---|
| **Epidemics** (Elsevier) | **free** under the subscription route | strong — publishes methods and inference work on epidemic models | Hybrid: open access is optional and expensive; the default subscription route costs nothing. Best first choice. |
| **Mathematical Biosciences** (Elsevier) | **free** under subscription | decent — model identifiability and estimation | Slower, more mathematical readership. |
| **PLOS Global Public Health** | likely **free** — Research4Life countries, incl. Pakistan, are covered by PLOS's Global Equity agreement | reasonable — global health methodology | Verify eligibility with the journal before submitting; it keys on the corresponding author's institution, and "independent researcher" is an untested case. |
| **PLOS Neglected Tropical Diseases** | APC, waivers by application | **best topical fit** — it published the 99-model review this paper answers | The waiver is discretionary, not automatic. Apply at submission, not after. |
| **Royal Society Open Science** | APC, waivers available | good — broad, methods-friendly | |
| **Infectious Disease Modelling** (KeAi) | **USD 1,100** | strong topical fit | Fully open access, no subscription route. Only if funded. |
| **BMC Infectious Diseases** | APC | acceptable | Waiver policy exists but is tight. |

**Recommended order:** Epidemics → PLOS NTD (with a waiver request) → Infectious
Disease Modelling. Epidemics first because it costs nothing, the fit is genuine,
and a rejection there still leaves the preprint doing its work.

## Do not submit to

Any journal that solicits by email, promises review in under two weeks, or is not
listed in DOAJ and indexed in Scopus or PubMed. A predatory venue on a CV is worse
than no publication: reviewers of scholarship applications recognise the names,
and it converts a strength into a question about judgement.

Check any candidate against:

* DOAJ (`doaj.org`) — is it a listed open-access journal?
* Scopus source list / PubMed — is it actually indexed?
* Whether the editorial board members exist and list the journal on their own
  university pages.

**Note on PUJM** (Punjab University Journal of Mathematics): it introduced APCs in
January 2025, so it is no longer the free local option it used to be, and its fit
for an epidemiological methods paper was always weak.

## Before submitting anywhere

- [ ] Fill in the repository URL in `paper/paper.tex` (currently the one
      deliberate placeholder — `tests/test_paper_consistency.py` enforces that it
      is the only one).
- [ ] Tag a release and mint the Zenodo DOI, so the paper can cite the exact code
      that produced its numbers.
- [ ] Re-run `scripts/23_paper_numbers.py` and the test suite; the paper has
      drifted from the pipeline twice already.
- [ ] For a PLOS venue, include `paper/author_summary.tex` (192 words, written
      2026-08-03) and confirm the current word limit.
- [ ] Have someone read the manuscript who has not been near it. Every
      collaborator on this project is the same person.

## Manuscript size against venue (added 2026-08-03)

The paper is now **29 pages, ~14,400 words of body text, 8 figures, 20 tables, 21
references**, with a 294-word abstract. That is long, and length is the first thing an
editor weighs against fit.

| Venue | Limit | Where this stands |
|---|---|---|
| PLOS NTD / PLOS Global Public Health | no hard word limit; abstract **300 words** | abstract fits at 294 with nothing to spare. Both need a separate non-technical **Author Summary** (~150–200 words) that does not yet exist. |
| Epidemics | no hard limit, but typical papers are 6,000–8,000 words | **roughly double**. Expect to move material to supplementary. |
| Royal Society Open Science | no limit | fine as is. |

**What to move to supplementary if a venue pushes back**, in the order it should go —
each of these supports the paper without being load-bearing:

1. The appendix case study (identifiability, aggregation bias, two-patch model) — already
   separated for a different reason, and the natural first cut.
2. The alternative-factor-list tables (all 63 subsets) — keep the leave-one-out table and
   the sceptic's bound in the main text, move the rest.
3. The criteria-comparison table — the sentence "stricter is worse, monotonically" carries
   it; the table is confirmation.
4. The climate-grid-cell and optimiser checks — both are bounds, both belong in a robustness
   supplement.

Do **not** cut: the three-way partition, the coverage study, the operating characteristics,
or the paper's own answer under its own rule. Those are the four things the paper is for.

## What changed after the review passes (2026-08-03)

Worth knowing before writing a cover letter, because the contribution is no longer only the
critique:

* The headline is now a **three-way partition** — outbreak 23.5%, analysis 15.9%,
  **interaction 60.6%** — not a two-way split. The interaction is the finding: it is what
  no standardised convention can reach, and it is what a specification curve, which shows
  the analysis main effect, cannot display. That last point is a statement about the method
  of Steegen et al. and Simonsohn et al., which makes the paper partly a methods
  contribution and widens the plausible venue list.
* The same partition applied to the **simulation** gives a fingerprint that discriminates:
  the analysis main effect is 47% where no effect exists and 10% where one does. The real
  data reads 15.9%. This recovers, on firmer ground, an inference an earlier draft had to
  withdraw.
* The **interaction is largely unpredictable** — 34% of the observation model's variation is
  explained by pre-fit outbreak descriptors, a median of 1.8% for the other five factors.
* Twelve references were added and one was corrected (wrong journal, unverifiable
  quotation). Two direct quotations were paraphrases wearing quotation marks and are now the
  sources' own words. `test_every_quotation_has_been_checked_against_its_source` holds the
  verified list.

## Still to do, unchanged

The four items in the checklist above are all still open, and the last one — an outside
reader — is the one that matters most and the one no amount of further self-review can
substitute for.
