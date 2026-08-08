"""Tests for the instability definition shared by the analysis steps.

The first test here is the one that would have caught the bug this module was
written to prevent: an instability rule with the factorial size hard-coded to 24
reports the right answer on a 24-combination design and the wrong one on a
48-combination design, without failing. Every test below therefore checks the
same property at two different design sizes.
"""

from __future__ import annotations

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk.robustness import (
    benjamini_hochberg,
    complete_windows,
    instability_share,
    is_unstable,
    pairwise_disagreement,
    pairwise_disagreement_at_margin,
    two_way_decomposition,
    variance_decomposition,
    window_verdicts,
)

import numpy as np
import pandas as pd
import pytest


def make_frame(wins_per_window: list[int], n_combos: int) -> pd.DataFrame:
    """A factorial result table in which window i has wins_per_window[i] fits
    favouring the climate model, out of n_combos."""
    rows = []
    for i, wins in enumerate(wins_per_window):
        for j in range(n_combos):
            rows.append(dict(country="X", unit=f"u{i}", window_start="2015-01-01",
                             climate_wins=int(j < wins), weeks=40, cases=5000))
    return pd.DataFrame(rows)


@pytest.mark.parametrize("n_combos", [24, 48])
def test_unanimous_windows_are_stable(n_combos):
    """All wins or no wins is unanimous at any design size."""
    assert not is_unstable(np.array([0, n_combos]), n_combos).any()


@pytest.mark.parametrize("n_combos", [24, 48])
def test_any_split_is_unstable(n_combos):
    """Every split verdict counts, including ones near unanimity.

    A rule of the form `between(1, 23)` passes this at n_combos=24 and fails at
    48, where wins of 24 to 47 are split verdicts that it would call unanimous.
    """
    wins = np.arange(1, n_combos)
    assert is_unstable(wins, n_combos).all()


@pytest.mark.parametrize("n_combos", [24, 48])
def test_share_matches_a_direct_count(n_combos):
    """The share agrees with counting split windows by hand."""
    wins = [0, n_combos, 1, n_combos - 1, n_combos // 2, n_combos - 2]
    frame = make_frame(wins, n_combos)
    expected = sum(0 < w < n_combos for w in wins) / len(wins)
    assert instability_share(frame, n_combos) == pytest.approx(expected)


@pytest.mark.parametrize("n_combos", [24, 48])
def test_verdict_categories_partition_the_windows(n_combos):
    """Always, never and unstable are mutually exclusive and cover everything."""
    frame = make_frame([0, n_combos, 3, n_combos - 3], n_combos)
    per = window_verdicts(frame, n_combos)
    total = per["unstable"] + per["always_climate"] + per["never_climate"]
    assert (total == 1).all()


def test_complete_windows_drops_partial_ones_and_reads_the_size():
    """The factorial size comes from the data, and short windows are dropped."""
    frame = make_frame([10, 20], 48)
    # Remove nine fits from the second window, leaving it incomplete.
    drop = frame.index[frame["unit"] == "u1"][:9]
    frame = frame.drop(index=drop)

    complete, n_combos = complete_windows(frame)
    assert n_combos == 48
    assert set(complete["unit"]) == {"u0"}


def test_a_hardcoded_bound_of_24_would_fail_here():
    """Documents the exact failure the module exists to prevent.

    On a 48-combination design, the old rule counts 36-of-48 as unanimous. This
    test asserts the two rules disagree, so that if anyone reintroduces the
    hard-coded bound the difference is visible rather than silent.
    """
    n_combos = 48
    wins = pd.Series([36, 40, 47])
    old_rule = wins.between(1, 23)
    assert not old_rule.any(), "the old rule sees these as unanimous"
    assert is_unstable(wins, n_combos).all(), "they are in fact all split"


# --- pairwise disagreement, the design-size-invariant statistic --------------

@pytest.mark.parametrize("n_combos", [24, 48, 144])
def test_unanimous_windows_never_disagree(n_combos):
    frame = make_frame([0, n_combos], n_combos)
    assert pairwise_disagreement(frame, n_combos).eq(0.0).all()


@pytest.mark.parametrize("n_combos", [24, 48, 144])
def test_even_split_approaches_one_half(n_combos):
    """An evenly split outbreak is the maximum: two analyses drawn from it
    agree half the time by chance alone."""
    frame = make_frame([n_combos // 2], n_combos)
    value = float(pairwise_disagreement(frame, n_combos).iloc[0])
    assert value == pytest.approx(0.5)


def test_it_does_not_move_with_the_size_of_the_design():
    """The point of the statistic: the same proportion split gives the same
    answer whether the factorial has 24 cells or 144.

    The `unstable` flag does not have this property — it reads as 1.0 for both
    designs while saying nothing about how lopsided the split is — which is why
    a headline built on it alone invites the objection this statistic answers.
    """
    quarter_24 = pairwise_disagreement(make_frame([6], 24), 24).iloc[0]
    quarter_144 = pairwise_disagreement(make_frame([36], 144), 144).iloc[0]
    assert quarter_24 == pytest.approx(quarter_144, abs=1e-12)
    assert quarter_24 == pytest.approx(2 * 0.25 * 0.75)


def test_it_is_ordered_by_how_lopsided_the_split_is():
    """More lopsided means less disagreement, monotonically."""
    n = 144
    values = [float(pairwise_disagreement(make_frame([w], n), n).iloc[0])
              for w in (1, 12, 36, 72)]
    assert values == sorted(values), "should rise toward an even split"


# --- the same statistic under an abstention rule -----------------------------

def make_delta_frame(deltas_per_window: list[list[float]]) -> pd.DataFrame:
    """A factorial table with given delta_aic values per window."""
    rows = []
    for i, deltas in enumerate(deltas_per_window):
        for dv in deltas:
            rows.append(dict(country="X", unit=f"u{i}", window_start="2015-01-01",
                             delta_aic=float(dv)))
    return pd.DataFrame(rows)


def test_margin_ignores_the_analyses_that_abstain():
    """Only combinations outside the band count toward disagreement."""
    # Two decisive for climate, two decisive against, four inside the band.
    frame = make_delta_frame([[-10, -20, 10, 20, 1, -1, 2, -2]])
    at_zero = pairwise_disagreement_at_margin(frame, 0.0).iloc[0]
    at_four = pairwise_disagreement_at_margin(frame, 4.0).iloc[0]
    # At margin 0 all eight speak, four each way: an even split.
    assert at_zero == pytest.approx(0.5)
    # At margin 4 only the four decisive ones speak, still even.
    assert at_four == pytest.approx(0.5)


def test_a_margin_that_removes_dissent_drives_it_to_zero():
    """Where dissent lives only inside the band, abstaining removes it."""
    frame = make_delta_frame([[-10, -20, -30, 1, 2, 3]])
    assert pairwise_disagreement_at_margin(frame, 0.0).iloc[0] > 0
    assert pairwise_disagreement_at_margin(frame, 4.0).iloc[0] == pytest.approx(0.0)


def test_silent_windows_are_dropped_not_counted_as_agreeing():
    """A rule that says nothing has not agreed with itself.

    Counting silence as agreement is how a remedy flatters itself: push the
    margin high enough and every window falls inside the band, which would read
    as perfect stability.
    """
    frame = make_delta_frame([[-10, 10], [1, -1]])   # second window all inside
    result = pairwise_disagreement_at_margin(frame, 4.0)
    assert len(result) == 1, "the silent window should be absent, not zero"


@pytest.mark.parametrize("n_combos", [24, 144])
def test_margin_version_is_also_design_invariant(n_combos):
    """Same proportions among speakers give the same answer at any size."""
    half = n_combos // 2
    deltas = [-10.0] * (half // 2) + [10.0] * (half // 2) + [0.5] * half
    value = pairwise_disagreement_at_margin(make_delta_frame([deltas]), 4.0).iloc[0]
    assert value == pytest.approx(0.5)


# --- variance decomposition: setting versus analyst ---------------------------

def test_all_variation_within_when_outbreaks_are_identical():
    """If every outbreak has the same split, nothing varies between them."""
    frame = make_frame([12, 12, 12], 24)
    r = variance_decomposition(frame)
    assert r["between_share"] == pytest.approx(0.0, abs=1e-9)
    assert r["within_share"] == pytest.approx(1.0, abs=1e-9)


def test_all_variation_between_when_each_outbreak_is_unanimous():
    """If every outbreak agrees with itself, nothing varies within them."""
    frame = make_frame([0, 24, 0, 24], 24)
    r = variance_decomposition(frame)
    assert r["within_share"] == pytest.approx(0.0, abs=1e-9)
    assert r["between_share"] == pytest.approx(1.0, abs=1e-9)


def test_shares_sum_to_one():
    """The decomposition is exact for a binary outcome."""
    frame = make_frame([3, 12, 20, 24, 0], 24)
    r = variance_decomposition(frame)
    assert r["between_share"] + r["within_share"] == pytest.approx(1.0)


@pytest.mark.parametrize("n_combos", [24, 144])
def test_the_split_does_not_move_with_the_design_size(n_combos):
    """The property that makes this quotable.

    Doubling the number of analyses at the same proportions adds to both the
    numerator and the denominator, so the share is unchanged — unlike "some pair
    disagrees", which rises mechanically. Measured on this study's three designs
    the figure reads 78.3%, 79.7% and 76.5%.
    """
    wins = [int(n_combos * f) for f in (0.25, 0.5, 0.75, 1.0, 0.0)]
    r = variance_decomposition(make_frame(wins, n_combos))
    # Worked by hand for these proportions: within = mean p(1-p) = 0.125,
    # between = var([.25,.5,.75,1,0]) = 0.125, total = 0.25, so the split is
    # exactly half and half at any design size.
    assert r["within_share"] == pytest.approx(0.5, abs=1e-9)
    assert r["total"] == pytest.approx(0.25, abs=1e-9)


def test_a_constant_outcome_does_not_divide_by_zero():
    """No variation at all must return cleanly rather than raise."""
    frame = make_frame([0, 0], 24)
    r = variance_decomposition(frame)
    assert r["total"] == 0.0


# --- two-way: outbreak, analysis, and the interaction -------------------------

CELLS = ["observation", "temp_form"]


def make_crossed(values: dict[tuple[str, str], int]) -> pd.DataFrame:
    """A complete crossed table. Keys are (outbreak, analysis) pairs."""
    rows = [dict(country="X", unit=u, window_start="2015-01-01",
                 observation=a.split("/")[0], temp_form=a.split("/")[1],
                 climate_wins=v)
            for (u, a), v in values.items()]
    return pd.DataFrame(rows)


def grid(fn) -> pd.DataFrame:
    """Two outbreaks by four analyses, filled by fn(outbreak_i, analysis_j)."""
    obs, forms = ["nb", "poisson"], ["briere", "loglin"]
    return make_crossed({(f"u{i}", f"{o}/{f}"): fn(i, 2 * oi + fi)
                         for i in range(2)
                         for oi, o in enumerate(obs)
                         for fi, f in enumerate(forms)})


def test_the_three_shares_sum_to_one():
    """The partition is exact on a complete crossed design."""
    r = two_way_decomposition(grid(lambda i, j: (i + j) % 2), CELLS)
    assert r["check"] == pytest.approx(1.0)


def test_variation_only_between_outbreaks():
    """Every analysis agrees; the outbreaks differ."""
    r = two_way_decomposition(grid(lambda i, j: i), CELLS)
    assert r["outbreak"] == pytest.approx(1.0)
    assert r["analysis"] == pytest.approx(0.0, abs=1e-12)
    assert r["interaction"] == pytest.approx(0.0, abs=1e-12)


def test_variation_only_between_analyses():
    """Every outbreak agrees; the analyses differ. A house style would fix this."""
    r = two_way_decomposition(grid(lambda i, j: int(j >= 2)), CELLS)
    assert r["analysis"] == pytest.approx(1.0)
    assert r["outbreak"] == pytest.approx(0.0, abs=1e-12)
    assert r["interaction"] == pytest.approx(0.0, abs=1e-12)


def test_pure_interaction_has_no_main_effects():
    """The case the one-way split cannot distinguish, and the reason this
    function exists.

    Each analysis endorses exactly half the outbreaks and each outbreak is
    endorsed by exactly half the analyses, so neither margin carries any signal
    at all — yet which analysis endorses which outbreak is fully determined.
    Standardising the analysis would remove none of this. The one-way
    decomposition reads it as 100% "within outbreak", which an earlier draft of
    the paper glossed as attributable to the analyst; it is not.
    """
    r = two_way_decomposition(grid(lambda i, j: (i + j) % 2), CELLS)
    assert r["interaction"] == pytest.approx(1.0)
    assert r["outbreak"] == pytest.approx(0.0, abs=1e-12)
    assert r["analysis"] == pytest.approx(0.0, abs=1e-12)

    one_way = variance_decomposition(grid(lambda i, j: (i + j) % 2),
                                     key=["country", "unit", "window_start"])
    assert one_way["within_share"] == pytest.approx(1.0)


def test_a_constant_result_returns_cleanly():
    """No variation must not divide by zero."""
    r = two_way_decomposition(grid(lambda i, j: 1), CELLS)
    assert np.isnan(r["outbreak"])


# --- Benjamini-Hochberg, and the wrong version that passes everything --------

def test_one_strong_result_does_not_carry_the_whole_family():
    """The failure that motivated putting this in the package.

    Taking the running minimum in the wrong direction lets the single smallest
    q-value propagate to every larger p-value, so one genuine finding certifies
    all its neighbours. Step 37 printed exactly that: one correlation at
    p = 2e-6 among 23 nulls, and 24 of 24 "surviving" correction.
    """
    p = np.array([2e-6] + [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95])
    q = benjamini_hochberg(p)
    assert q[0] < 0.05, "the real finding should survive"
    assert (q[1:] > 0.05).all(), "the nulls must not ride along with it"

    wrong = pd.Series(np.sort(p) * len(p) / np.arange(1, len(p) + 1)).cummin()
    assert (wrong < 0.05).all(), "the wrong direction passes the whole family"


def test_a_family_of_pure_noise_survives_nothing():
    rng = np.random.default_rng(0)
    assert (benjamini_hochberg(rng.uniform(size=200)) < 0.05).sum() == 0


def test_q_values_are_never_smaller_than_their_p_values():
    """Correction can only make a result less significant, never more."""
    p = np.array([0.001, 0.01, 0.02, 0.04, 0.3, 0.9])
    assert (benjamini_hochberg(p) >= p - 1e-12).all()


def test_it_is_monotone_in_the_p_values():
    """A larger p-value can never earn a smaller q-value."""
    p = np.array([0.001, 0.004, 0.01, 0.02, 0.2, 0.8])
    q = benjamini_hochberg(p)
    assert list(q) == sorted(q)


def test_the_answer_does_not_depend_on_the_input_order():
    """q-values come back aligned to the p-values as given."""
    p = np.array([0.2, 0.001, 0.9, 0.01])
    q = benjamini_hochberg(p)
    shuffled = np.array([1, 3, 0, 2])
    assert np.allclose(benjamini_hochberg(p[shuffled]), q[shuffled])


def test_a_single_test_is_left_alone():
    """With one hypothesis there is nothing to correct for."""
    assert benjamini_hochberg([0.03])[0] == pytest.approx(0.03)


def test_it_matches_a_worked_example():
    """Benjamini and Hochberg's own arithmetic, m = 4."""
    p = np.array([0.01, 0.02, 0.03, 0.04])
    # raw p*m/rank = 0.04, 0.04, 0.04, 0.04 -> all 0.04 after monotonicity
    assert benjamini_hochberg(p) == pytest.approx([0.04, 0.04, 0.04, 0.04])


def test_an_empty_family_returns_empty():
    assert len(benjamini_hochberg([])) == 0
