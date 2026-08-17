# Zenodo manuscript record — har field paste-ready

Likha gaya 2026-08-12, MetaArXiv ke scope-reject ke baad. Zenodo par koi moderation
nahi hai: publish karte hi DOI milta hai, usi din.

⚠️ **Do alag records hain — confuse mat karna:**

| Record | DOI | Kya hai |
|---|---|---|
| Pehle se mojood | `10.5281/zenodo.21854302` | **Code** (GitHub release v1.1.0 ka archive) |
| **Ab banana hai** | naya milega | **Manuscript** (paper.pdf + supplementary.pdf) |

GitHub wale record ko haath nahi lagana. Naya upload banana hai.

⚠️ **Publish permanent hai** — file delete nahi hoti, sirf naya version aa sakta hai.
Is liye upload se pehle confirm: `paper.pdf` 24 safhat wala Times-font version hai
(equations theek, koi em-dash nahi).

---

## Files (dono upload karni hain)

```
D:\A Scholarships\dengue-pakistan\paper\paper.pdf
D:\A Scholarships\dengue-pakistan\paper\supplementary.pdf
```

## Form fields

**Resource type:** `Publication` → **`Preprint`**

**Title:**
```
How much of the evidence for climate-driven dengue transmission survives the analyst? A multiverse analysis of 221 outbreaks in 33 countries
```

**Publication date:** aaj ki tareekh

**Creators:**
| Field | Value |
|---|---|
| Name | `Dilber, Amna` |
| ORCID | `0009-0008-5684-4516` |
| Affiliation | *khali chhor do* |

**Description:** (abstract — `paper/abstract_for_form.txt` se, ya neeche se copy karo)
```
Compartmental models fitted to surveillance data are routinely used to test whether dengue transmission is climate-driven, and reviews explain their conflicting findings by appealing to local context. We fit every usable outbreak in a global surveillance compilation under all 144 combinations of six analysis choices that published models rarely state, so that no choice varies with the setting because none varies at all: 33,152 fits, of which the 221 outbreaks in 33 countries completing every combination are analysed. Roughly three-quarters of the variation in whether climate forcing is endorsed lies between analyses of the same outbreak, not between outbreaks (76.5%; 95% CI 72.6-81.1 resampling whole countries). Partitioned three ways the largest term is neither: 61% is an outbreak-by-analysis interaction, against 16% for the analysis itself, so no standardised convention can reach it and a specification curve, which displays the main effect, cannot show it. Adding climate covariates moves estimated R0 the same way, its interaction share rising from 23% to 69%. The intervals do not keep their promise. A nominal 95% interval contains the true transmission coefficient 78% of the time and the true rainfall coefficient 71% (42% when the epidemic is a sum of two asynchronous ones, as these data demonstrably are). Combining eight analyses by Rubin's rules restores nominal coverage at 1.5-2.2 times the width, reaching 90% rather than 95% under structural mis-specification. Applying our own rule gives a clear verdict for 137 of 221 outbreaks, 86% of which support climate forcing, and the three-way partition of the real data matches simulated outbreaks containing a real effect rather than simulated outbreaks containing none. The claim is not that the effect is absent, but that in 38% of outbreaks a single season cannot settle the question.

All analysis code, result tables and figures: https://github.com/amnadilber/dengue-multiverse (archived at https://doi.org/10.5281/zenodo.21854302). The supplementary PDF (Sections S1-S6) is included in this record.
```

**License:** `Creative Commons Attribution 4.0 International`

**Keywords** (ek ek kar ke):
```
dengue
multiverse analysis
specification curve
researcher degrees of freedom
reproducibility
model selection
compartmental models
epidemiology
climate and health
```

**Related works / Related identifiers** (agar section mile):
| Relation | Identifier |
|---|---|
| `Is supplemented by` | `10.5281/zenodo.21854302` (DOI) |
| `Is supplemented by` | `https://github.com/amnadilber/dengue-multiverse` (URL) |

**Language:** English

**Communities, Funding, Conference:** sab khali chhor do.

---

## ✅ HO GAYA — 17 Aug 2026

| | DOI |
|---|---|
| **Version v1** (citations mein yehi) | `10.5281/zenodo.21984527` |
| All-versions (hamesha latest par jata hai) | `10.5281/zenodo.21984526` |

## Publish ke baad — DOI mujhe dena (mukammal)

Main phir:
1. `CITATION.cff` mein `preferred-citation` add karungi
2. README mein preprint DOI + badge
3. `OUTREACH_EMAILS.md` ke har `[DOI]` placeholder mein daal doongi
4. CV entry likh doongi
5. Journal cover letter mein daal kar Epidemics ka package tayyar karungi
