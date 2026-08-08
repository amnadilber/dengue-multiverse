"""The paper must quote the numbers the pipeline actually produced.

The manuscript has drifted from the results twice. Once when the factorial grew
from four factors to five and the prose kept the old figures; once when a bound
tied to the old design size made a script report 11.5% where the truth was 88.0%,
and the wrong figure nearly reached the text. Both were caught by reading, which
does not scale and does not repeat.

These tests recompute the headline quantities from the stored result tables and
assert that the manuscript contains them. They are deliberately literal: each
checks for a specific string, so a number that changes in the pipeline fails here
with the exact substring the paper is missing.

The tests skip rather than fail when a result table is absent, so a fresh clone
without `results/` still passes its suite. They fail loudly when the tables exist
and disagree with the paper, which is the case that matters.
"""

from __future__ import annotations

from pathlib import Path

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve
from dengue_pk.robustness import (complete_windows, latest_factorial,
                                  window_verdicts)

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "paper.tex"
README = ROOT / "README.md"

cfg = load_config()
TABLES = resolve(cfg, "tables")


def _robustness() -> pd.DataFrame:
    """The same factorial the analysis steps read.

    Deliberately routed through `latest_factorial` rather than naming files
    here. An earlier version of this helper carried its own list, so when the
    six-factor table appeared these tests went on checking the paper against the
    five-factor one — passing where they should have failed, and failing where
    the paper had correctly been updated. A consistency test with its own idea of
    which results are current is worse than no test.
    """
    try:
        return pd.read_csv(latest_factorial(TABLES))
    except FileNotFoundError:
        pytest.skip("no global robustness table present")


def _table(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    return pd.read_csv(path)


def _text(path: Path) -> str:
    if not path.exists():
        pytest.skip(f"{path.name} not present")
    return path.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def verdicts():
    d, n_combos = complete_windows(_robustness())
    return window_verdicts(d, n_combos), n_combos


@pytest.mark.parametrize("doc", [PAPER, README])
def test_headline_instability_is_quoted(verdicts, doc):
    """The share of windows whose verdict changes, to one decimal place."""
    per, _ = verdicts
    pct = f"{per['unstable'].mean() * 100:.1f}"
    assert pct in _text(doc), f"{doc.name} does not contain {pct}% instability"


@pytest.mark.parametrize("doc", [PAPER, README])
def test_factorial_size_is_quoted(verdicts, doc):
    """The number of combinations, which changed once and broke a script."""
    _, n_combos = verdicts
    assert str(n_combos) in _text(doc)


def test_fit_count_is_quoted():
    """The total number of fits, in the paper's LaTeX thousands form."""
    n = len(_robustness())
    latex = f"{n // 1000}{{,}}{n % 1000:03d}"
    assert latex in _text(PAPER), f"paper does not contain {latex} fits"


def test_window_counts_are_quoted(verdicts):
    """Windows always favouring climate, and windows that split."""
    per, _ = verdicts
    text = _text(PAPER)
    assert str(int(per["always_climate"].sum())) in text
    assert str(int(per["unstable"].sum())) in text


@pytest.mark.parametrize("doc", [PAPER, README])
def test_every_factor_flip_rate_is_quoted(doc):
    """Each of the five choices, at the rate the sensitivity table gives."""
    sens = _table("15_choice_sensitivity.csv")
    text = _text(doc)
    missing = [f"{r['factor']}/{r['comparison']} = {r['flipped'] * 100:.1f}%"
               for _, r in sens.iterrows()
               if f"{r['flipped'] * 100:.1f}" not in text]
    assert not missing, f"{doc.name} is missing flip rates: {missing}"


def test_the_remedy_is_quoted():
    """The margin-4 rule: its instability and how often it still answers."""
    rules = _table("16_stability_rules.csv")
    row = rules[rules["rule"] == "|dAIC| > 4"]
    if row.empty:
        pytest.skip("margin-4 rule not in the stability table")
    row = row.iloc[0]
    text = _text(PAPER)
    assert f"{row['unstable_pct']:.1f}" in text
    assert f"{row['answers_pct']:.1f}" in text


def test_criterion_comparison_is_quoted():
    """Every criterion's instability, so 'stricter is worse' stays true."""
    crit = _table("17_criteria_comparison.csv")
    text = _text(PAPER)
    missing = [r["criterion"] for _, r in crit.iterrows()
               if f"{r['unstable_pct']:.1f}" not in text]
    assert not missing, f"paper is missing criteria: {missing}"


def test_simulation_error_rates_are_quoted():
    """False positives and power, for the rules the paper tabulates.

    Added after the simulation was rerun at a different sample size with its
    selection bias removed: the section quoting it is the one most likely to be
    left stale, because it is the furthest downstream.
    """
    oc = _table("29_operating_characteristics_6factor.csv")
    text = _text(PAPER)
    missing = []
    for _, r in oc[oc["margin"].isin([0, 4])].iterrows():
        for col in ("false_positive_pct", "power_pct"):
            if f"{r[col]:.1f}" not in text:
                missing.append(f"{r['observation']} margin {r['margin']:.0f} "
                               f"{col}={r[col]:.1f}")
    assert not missing, f"paper is missing simulation rates: {missing}"


def test_simulation_scope_is_quoted():
    """The number of simulated fits, in the paper's LaTeX thousands form."""
    path = TABLES / "28_false_positive_6factor.csv"
    if not path.exists():
        pytest.skip("simulation table not present")
    n = len(pd.read_csv(path))
    latex = f"{n // 1000}{{,}}{n % 1000:03d}"
    assert latex in _text(PAPER), f"paper does not contain {latex} simulated fits"


def test_two_way_decomposition_is_quoted():
    """Outbreak, analysis and interaction, at every design size the paper shows.

    The interaction is the term the paper's own recommendation now rests on, and
    it is computed nowhere else, so a drift here would be silent.
    """
    tw = _table("43_two_way_decomposition.csv")
    text = _text(PAPER)
    missing = []
    for _, r in tw.iterrows():
        if pd.isna(r.get("n_combos")):
            continue
        for col in ("outbreak", "analysis", "interaction"):
            if f"{r[col] * 100:.1f}" not in text:
                missing.append(f"{r['n_combos']:.0f} cells {col}="
                               f"{r[col] * 100:.1f}")
    assert not missing, f"paper is missing decomposition shares: {missing}"


def test_the_three_shares_still_sum_to_one():
    """The partition is exact; a check column drifting from 1 means the design
    stopped being complete and every share above is then meaningless."""
    tw = _table("43_two_way_decomposition.csv")
    rows = tw[tw["n_combos"].notna()] if "n_combos" in tw.columns else tw
    for _, r in rows.iterrows():
        assert r["check"] == pytest.approx(1.0, abs=1e-9)


def test_withdrawn_explanations_do_not_return():
    """Five interpretive claims were measured and found wrong. Each was replaced
    by a sentence that sounds less tidy, which is exactly the kind of edit a
    later revision undoes by accident.
    """
    text = _text(PAPER)
    for phrase, why in [
        ("flexible enough to chase noise",
         "Poisson does not enlarge the coefficients; it enlarges the criterion"),
        ("most studies use Poisson",
         "the review reports 18.3%, a large minority"),
        ("more likely to differ because the analysts differed",
         "measured the other way round: 29.8% against 32.7%"),
        ("the sign\nis unambiguous",
         "cold starts gave slightly less disagreement, not more"),
        ("two orders of magnitude",
         "the median index of dispersion is 5.7, not 100"),
    ]:
        assert phrase not in text, f"withdrawn claim is back: {phrase!r} — {why}"


def test_interaction_predictability_is_quoted():
    """The share of each factor's effect that outbreak descriptors explain.

    This is the result that turns the interaction from an unexplained remainder
    into a measured one, and it lives in a single table with no other check on
    it.
    """
    res = _table("46_interaction_structure.csv")
    text = _text(PAPER)
    top = res.sort_values("p").iloc[0]
    assert top["factor"] == "observation", (
        "the strongest association should still be the observation model against "
        f"dispersion, not {top['factor']}/{top['descriptor']}")
    assert f"{abs(top['rho']):.2f}" in text, (
        f"paper does not quote rho = {top['rho']:.2f}")
    n_fdr = int((res["q"] < 0.05).sum())
    assert str(n_fdr) in text, f"paper does not state that {n_fdr} survive FDR"


def test_decomposition_under_truth_is_quoted():
    """Both simulated arms and the real data, all three components each."""
    tw = _table("47_decomposition_under_truth.csv")
    text = _text(PAPER)
    missing = [f"{r['source']} {c}={r[c] * 100:.1f}"
               for _, r in tw.iterrows()
               if "raw" not in str(r["source"])
               for c in ("outbreak", "analysis", "interaction")
               if f"{r[c] * 100:.1f}" not in text]
    assert not missing, f"paper is missing: {missing}"


def test_the_null_arm_is_still_the_one_with_the_large_analysis_effect():
    """The claim the fingerprint rests on, asserted as a property not a number.

    If a rerun ever inverted this, every sentence built on it would be wrong
    while every quoted figure still matched its table.
    """
    tw = _table("47_decomposition_under_truth.csv").set_index("source")
    null = tw.loc["simulated, no climate effect"]
    effect = tw.loc["simulated, real effect"]
    real = tw.loc["real outbreaks"]
    assert null["analysis"] > effect["analysis"], (
        "the no-effect arm must have the larger analysis main effect")
    assert abs(real["analysis"] - effect["analysis"]) < abs(
        real["analysis"] - null["analysis"]), (
        "the real data must sit nearer the effect arm on the discriminating term")


def test_the_quoted_test_count_is_the_real_one(request):
    """The data-availability paragraph quotes the size of this suite.

    It said 157 for long enough that the suite had grown by half again, which is
    a small dishonesty of exactly the kind this project keeps finding in itself.
    `request.session.items` is the collected suite, so the claim now fails when
    it goes stale instead of quietly overstating or understating the work.

    Skipped when the run is a subset, since a subset's count is not the suite's.
    """
    collected = len(request.session.items)
    total = sum(1 for _ in (ROOT / "tests").glob("test_*.py"))
    if collected < 100 or total == 0:
        pytest.skip("partial collection; the count is only meaningful for a full run")
    assert str(collected) in _text(PAPER), (
        f"paper does not quote the suite size ({collected} tests collected)")


def test_the_quoted_script_count_is_the_real_one():
    """Same again for the number of numbered pipeline steps."""
    steps = sorted(p.name for p in (ROOT / "scripts").glob("*.py")
                   if p.name[:2].isdigit())
    highest = max(int(name[:2]) for name in steps)
    assert f"{highest} numbered pipeline" in _text(PAPER), (
        f"paper does not quote {highest} numbered pipeline scripts")


#: Every direct quotation of a source in the manuscript, checked word for word
#: against the source itself on 2026-08-03. Two of the four were wrong when this
#: list was first built: one attributed to Khamthong and Phramrung a sentence
#: they do not contain, and one dropped "of the" from the Leung review. Both were
#: paraphrases that had acquired quotation marks, which is the easiest kind of
#: citation error to commit and the hardest for a co-author to catch.
VERIFIED_QUOTATIONS = {
    "[t]he reporting of methodology and model performance measures were inadequate"
    "\nin many of the existing prediction models":
        "Leung et al. 2023, PLOS NTD 17(2) e0010631, abstract",
    "forecasting performance is not determined solely by the\nstrength of marginal "
    "climate--dengue associations":
        "Khamthong & Phramrung 2026, PLOS NTD 20(4) e0014270, abstract",
    "statistically\nequivalent":
        "Khamthong & Phramrung 2026, discussion section 4.3",
}

#: Quotation marks used for the paper's own terms rather than a source. These
#: need no verification because no one else is being credited with them.
OWN_PHRASES = {"Between", "within", "climate affects", "outbreak"}


def test_every_quotation_has_been_checked_against_its_source():
    """A quotation mark is a factual claim about what someone else wrote.

    Nothing else in this suite can catch a misquotation, and a fabricated or
    drifted quote is the kind of error that ends a submission rather than
    delaying it. Any new quotation fails here until it is verified and listed.
    """
    import re
    text = _text(PAPER)
    found = re.findall(r"``(.*?)''", text, flags=re.DOTALL)
    # A regex that silently matches nothing would pass this test forever.
    assert len(found) >= len(VERIFIED_QUOTATIONS), (
        f"found only {len(found)} quotations; the manuscript should contain at "
        f"least the {len(VERIFIED_QUOTATIONS)} verified ones")
    unverified = [q for q in found
                  if q not in VERIFIED_QUOTATIONS
                  and not any(q.startswith(p) for p in OWN_PHRASES)]
    assert not unverified, (
        "unverified quotations in the paper — check each against its source and "
        f"add it to VERIFIED_QUOTATIONS: {unverified}")


def test_every_reference_is_cited_and_every_citation_defined():
    """A dangling \\cite compiles to a bold [?] and an uncited entry is padding.

    Both had happened: `mordecai` sat in the bibliography while the thermal
    limits it supplies were stated without attribution.
    """
    import re
    text = _text(PAPER)
    defined = set(re.findall(r"\\bibitem\{(\w+)\}", text))
    cited = set()
    for group in re.findall(r"\\cite\{([^}]*)\}", text):
        cited.update(k.strip() for k in group.split(","))
    assert not (defined - cited), f"never cited: {sorted(defined - cited)}"
    assert not (cited - defined), f"no bibliography entry: {sorted(cited - defined)}"


def test_the_bibliography_is_wide_enough_to_be_credible():
    """Not a style rule — a submission signal.

    A paper making claims about a literature while citing nine works reads as
    one that has not read it. The floor is deliberately low; the point is to
    fail if entries are ever deleted back toward the original nine.
    """
    import re
    n = len(re.findall(r"\\bibitem\{", _text(PAPER)))
    assert n >= 18, f"only {n} references"


def test_no_placeholders_remain():
    """Nothing in the manuscript is still waiting to be filled in."""
    text = _text(PAPER)
    for marker in ("TODO", "XXX", "FIXME", "??", "[number]", "\\todo"):
        assert marker not in text, f"paper still contains {marker!r}"


def test_no_prose_placeholders_remain():
    """Nothing bracketed is still waiting to be filled in.

    Until the repository was created this test allowed exactly one placeholder,
    `[repository URL]`, and failed on any other. That one is now filled, so the
    allowance is gone: any bracketed prose is an oversight.
    """
    text = _text(PAPER)
    brackets = [seg.split("]")[0] for seg in text.split("[")[1:] if "]" in seg]
    # LaTeX optional arguments (htbp, margin=2.5cm) and maths intervals
    # ([0.01, 20]) are bracketed too. A prose placeholder is words and spaces
    # and nothing else, which none of those are.
    prose = [b for b in brackets
             if " " in b and all(ch.isalpha() or ch == " " for ch in b)]
    assert not prose, f"unfilled placeholders: {prose}"


def test_the_repository_is_named_and_reachable_in_the_paper():
    """The data-availability statement must point somewhere real."""
    text = _text(PAPER)
    assert "github.com/amnadilber/dengue-multiverse" in text, (
        "the paper does not give the repository URL")
