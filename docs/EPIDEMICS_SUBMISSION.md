# Epidemics (Elsevier) — submission package

Prepared 2026-08-17, the day the preprint went live. Everything below is paste-ready for
Editorial Manager: <https://www.editorialmanager.com/epidemics/>

**Why Epidemics first:** free under the subscription route, and the fit is genuine — the
journal publishes inference and methods work on epidemic models, which is what this is.
A rejection costs nothing and the preprint keeps working meanwhile.

---

## What Editorial Manager will ask, in order

### 1. Article type
`Research Paper`

### 2. Title
```
How much of the evidence for climate-driven dengue transmission survives the analyst? A multiverse analysis of 221 outbreaks in 33 countries
```

### 3. Author
One author: **Amna Dilber**, `amnadilber.bi@gmail.com`, ORCID `0009-0008-5684-4516`.
Affiliation field: `Independent researcher, Lahore, Pakistan`.

### 4. Abstract
Paste from `paper/abstract_for_form.txt` (287 words, ASCII-safe).

### 5. Keywords
```
dengue; multiverse analysis; specification curve; model selection; researcher degrees of freedom; compartmental models; reproducibility
```

### 6. Highlights — Elsevier wants 3–5 bullets, **max 85 characters each** (counted)

```
221 dengue outbreaks fitted under all 144 combinations of six analysis choices
The verdict on climate forcing changes with the analyst in 92% of outbreaks
61% of the variation is an outbreak-by-analysis interaction no convention removes
Nominal 95% intervals cover the truth 71-78% of the time on simulated data
Combining eight analyses by Rubin's rules restores nominal coverage
```

### 7. Files to upload

| File | Item type in EM |
|---|---|
| `paper/paper.pdf` | Manuscript |
| `paper/supplementary.pdf` | Supplementary material |

(Initial submission par PDF kaafi hota hai; source `.tex` revision par maanga jayega.)

### 8. Cover letter

> Dear Editors,
>
> I am submitting "How much of the evidence for climate-driven dengue transmission
> survives the analyst? A multiverse analysis of 221 outbreaks in 33 countries" for
> consideration at Epidemics.
>
> Mechanistic models fitted to surveillance data are routinely used to test whether
> dengue transmission is climate-driven, and reviews explain their conflicting findings
> by appealing to local context. Whether that explanation is right has not been
> checkable, because each published estimate comes from one analysis of one place and
> the two cannot be separated after the fact. This study separates them: every usable
> outbreak in a global surveillance compilation is fitted under all 144 combinations of
> six analysis choices that published models rarely state — 33,152 fits.
>
> Roughly three-quarters of the variation in whether climate forcing is endorsed lies
> between analyses of the same outbreak rather than between outbreaks. Partitioned three
> ways, the largest term is neither: 61% is an outbreak-by-analysis interaction against
> 16% for the analysis itself. That distinction decides which remedy is available — a
> main effect can be settled by convention and an interaction cannot — and it is
> invisible to a specification curve, which displays the main effect and averages the
> interaction away.
>
> The paper does not stop at the diagnosis. On simulated data whose answer is known,
> conventional 95% intervals contain the truth 71–78% of the time and 42% when the
> epidemic is a sum of two asynchronous ones, which these data demonstrably are.
> Combining eight analyses by Rubin's rules restores nominal coverage at 1.5–2.2 times
> the width — a few minutes of computation, not a factorial. Applying that rule to the
> data gives a clear verdict for 137 of 221 outbreaks, 86% of which support climate
> forcing: the claim is not that the effect is absent, but that in 38% of outbreaks a
> single season cannot settle it.
>
> The work is relevant to Epidemics because it is a study of inference in epidemic
> models rather than a forecasting paper: the object of study is the estimator and the
> decision rule, not the epidemic. The main text is accompanied by a supplement holding
> the factor-list robustness analysis, the criterion comparison, further checks and a
> worked identifiability case study.
>
> All analyses are reproducible from a public repository containing 38 numbered pipeline
> steps, 269 tests — including automated checks that the manuscript's quoted figures
> still match the stored result tables — and a dated analysis log recording every
> formulation that failed and every error found. The manuscript is available as a
> preprint on Zenodo (https://doi.org/10.5281/zenodo.21984527); it is not under
> consideration elsewhere. I have no competing interests and received no funding. I am
> an independent researcher without institutional affiliation, and I would be glad of
> reviewers willing to be blunt.
>
> Yours sincerely,
> Amna Dilber
> ORCID 0009-0008-5684-4516

### 9. Declarations (EM ke sawal)

| Sawal | Jawab |
|---|---|
| Competing interests | None |
| Funding | This research received no specific grant from any funding agency |
| Data availability | statement paper mein hai; form mein: `Code and result tables at https://github.com/amnadilber/dengue-multiverse (archived: 10.5281/zenodo.21854302); raw data are public (OpenDengue v1.3, NASA POWER) and fetched by scripted download` |
| Generative AI declaration | **Yes, used** — manuscript mein section mojood hai; form mein wohi repeat: `A large language model (Anthropic Claude) assisted with implementing the analysis code, running the pipeline, and drafting and revising the manuscript. The author reviewed and verified all content and takes full responsibility for it.` |
| Preprint? | Yes — Zenodo, 10.5281/zenodo.21984527 |

### 10. Suggested reviewers (EM aksar 3 maangta hai)

Wohi log jo is literature mein hain aur jinka kaam paper cite karta hai — sab verified:

| Naam | Idara | Kyun munasib |
|---|---|---|
| Waqas Shabbir | Institute of Statistics, Alpen-Adria-Universität Klagenfurt | Pakistani dengue–climate GLM (BMC Public Health 2020) |
| Devin Kirk | (Kirk et al., PLOS Climate 2024) | temperature–dengue meta-analysis jise paper test karta hai |
| Kris Murray ya Oliver Brady | LSHTM | OpenDengue / dengue mapping |

*(Conflict nahi hai — kisi ke sath kabhi kaam nahi kiya. Ye normal practice hai.)*

### 11. Submit ke baad

- Status "With Editor" → weeks; "Under Review" → mahine. **"Under review at Epidemics"
  kehne ka haq submit hote hi hai.**
- Desk reject aaye to agla: **PLOS Global Public Health** (waiver check ke sath) —
  cover letter variant `COVER_LETTER.md` mein tayyar hai.

---

## CV entry (abhi likhne layak)

```
PREPRINTS
Dilber, A. (2026). How much of the evidence for climate-driven dengue transmission
survives the analyst? A multiverse analysis of 221 outbreaks in 33 countries.
Zenodo. https://doi.org/10.5281/zenodo.21984527
```

Journal submit hone ke baad us mein jor do: `(under review at Epidemics)`.

---

## SUPERSEDED — 2026-08-18

The submission questionnaire itself disclosed what this file had wrong: **Epidemics is
now a fully open access journal.** There is no subscription route. The APC is **USD
2,860**; Elsevier applies automatic Research4Life discounts, but Pakistan is a Group B
country (discount, not waiver), leaving roughly USD 1,430 — not payable here, and the
form requires the author to confirm payment upon acceptance, which would be false.

The submission was abandoned at the Additional Information step. Nothing false was
confirmed. The draft sits in Incomplete Submissions and can be deleted.

This was the third venue whose terms this project asserted without verification
(medRxiv's affiliation rule, MetaArXiv's scope, now Epidemics' fee model — the note
"free under the subscription route" was written 2026-07-27 and the journal had
converted). The files prepared here (cover letter, declaration, highlights) carry over
to the next venue with minor edits.

**Next venue: PLOS Global Public Health** — topical fit already rated good, PLOS's own
policy publishes Research4Life Group B authors free, and PLOS Publication Fee
Assistance can be requested at submission by authors without funds. Both are to be
verified ON PLOS'S OWN PAGES during package preparation, before anything is filled in.
