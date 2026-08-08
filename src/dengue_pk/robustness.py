"""Deciding whether a window's verdict is stable across the factorial.

This lives in the package rather than in a script because the definition is
shared by four analysis steps, and because getting it wrong is quiet. An earlier
version of step 18 asked whether the number of combinations favouring the
climate model fell ``between(1, 23)`` — correct for the 24-combination design it
was written against, and wrong the moment a fifth factor took the design to 48.
It did not fail; it reported 11.5% instability where the true figure was 88.0%,
because windows in which 24 to 47 of the 48 combinations favoured climate were
counted as unanimous. A wrong number that looks plausible is worse than a crash.

The definition itself is simple: a window is unstable when the verdict is not
unanimous — when the climate model wins under at least one combination and loses
under at least one other. Everything here is written so that the factorial size
is read from the data rather than assumed.
"""

from __future__ import annotations

import pathlib
from typing import Sequence

import numpy as np
import pandas as pd

#: The columns that jointly identify one outbreak window.
WINDOW_KEY: tuple[str, ...] = ("country", "unit", "window_start")

#: Factorial result tables, richest design first. Each analysis step reads the
#: first one present rather than naming a file, so growing the design means
#: adding one entry here instead of editing every downstream script — which is
#: how the four-to-five-factor change came to leave a hard-coded bound behind in
#: one of them.
FACTORIAL_TABLES: tuple[str, ...] = (
    "25_global_robustness_6factor.csv",
    "20_global_robustness_5factor.csv",
    "13_global_robustness.csv",
)


def latest_factorial(tables_dir) -> pathlib.Path:
    """Path to the richest factorial table available.

    Raises rather than returning None: an analysis step that silently ran on no
    data would print an empty summary and look like a result.
    """
    for name in FACTORIAL_TABLES:
        path = tables_dir / name
        if path.exists():
            return path
    raise FileNotFoundError(
        f"no factorial table in {tables_dir}; expected one of "
        f"{', '.join(FACTORIAL_TABLES)}")


def complete_windows(frame: pd.DataFrame,
                     key: Sequence[str] = WINDOW_KEY) -> tuple[pd.DataFrame, int]:
    """Restrict to windows that completed the full factorial.

    Returns the restricted frame and the factorial size, which is taken to be
    the largest number of fits any window achieved. Windows with fewer fits are
    dropped: a window that completed 30 of 48 combinations would look unstable
    or stable for reasons that have nothing to do with the analysis choices.
    """
    key = list(key)
    counts = frame.groupby(key).size()
    n_combos = int(counts.max())
    complete = counts[counts == n_combos].index
    return frame.set_index(key).loc[complete].reset_index(), n_combos


def is_unstable(wins: np.ndarray | pd.Series, n_combos: int) -> np.ndarray:
    """Whether each window's verdict changes across the factorial.

    ``wins`` counts the combinations favouring the climate model. Unstable means
    strictly between none and all of them, so the bound follows the design size
    instead of being written down.
    """
    w = np.asarray(wins)
    return (w >= 1) & (w <= n_combos - 1)


def window_verdicts(frame: pd.DataFrame, n_combos: int,
                    key: Sequence[str] = WINDOW_KEY) -> pd.DataFrame:
    """One row per window: how many combinations favoured climate, and whether
    that verdict is unanimous.

    ``climate_wins`` is expected to be a 0/1 column, one row per fit.
    """
    key = list(key)
    # `weeks` and `cases` are carried through when present because callers cut
    # the result by them, but they are descriptive rather than required: the
    # simulation tables have no such columns, and demanding them made this
    # function unusable on exactly the data the invariant statistic was needed
    # for.
    agg = {"wins": ("climate_wins", "sum")}
    for col in ("weeks", "cases"):
        if col in frame.columns:
            agg[col] = (col, "first")
    per = frame.groupby(key).agg(**agg).reset_index()
    per["unstable"] = is_unstable(per["wins"], n_combos)
    per["always_climate"] = per["wins"] == n_combos
    per["never_climate"] = per["wins"] == 0
    return per


def instability_share(frame: pd.DataFrame, n_combos: int,
                      key: Sequence[str] = WINDOW_KEY) -> float:
    """Share of windows whose verdict changes across the factorial."""
    return float(window_verdicts(frame, n_combos, key)["unstable"].mean())


def pairwise_disagreement(frame: pd.DataFrame, n_combos: int,
                          key: Sequence[str] = WINDOW_KEY) -> pd.Series:
    """Probability that two analyses of the same outbreak, drawn at random
    without replacement, reach opposite verdicts.

    Why this exists. The headline "the verdict changes in X% of outbreaks" is
    read off whether *any* two cells disagree, and that is easier to satisfy the
    more cells the design has: enlarging the factorial from 48 combinations to
    144 raises the figure even if nothing about the evidence changed. A reviewer
    is entitled to point that out, and the right answer is a statistic that does
    not move with the size of the design rather than a defence of one that does.

    With ``p`` the proportion of analyses favouring climate forcing, two drawn
    independently disagree with probability

        2 p (1 - p)

    which depends on the proportion alone and not at all on how finely the space
    was enumerated. It is zero for a unanimous outbreak and one half for an
    evenly split one — the maximum, since two analyses drawn from an even split
    agree half the time by chance.

    Drawing *without* replacement instead gives ``2 w (n - w) / [n (n - 1)]``,
    which is the more literal reading of "two different analyses" but carries a
    factor ``n / (n - 1)`` and so does creep with the design size: at a
    one-quarter split it reads 0.391 over 24 cells and 0.378 over 144. The
    difference is small and entirely an artefact of the enumeration, which is
    exactly the artefact this statistic exists to avoid, so the independent form
    is used and the finite-design correction deliberately dropped.

    Returned per window, so callers can average, bootstrap or cut it by
    covariate as they would any other per-window quantity.
    """
    per = window_verdicts(frame, n_combos, key)
    p = per["wins"].to_numpy(dtype=float) / float(n_combos)
    return pd.Series(2.0 * p * (1.0 - p), index=per.index)


def pairwise_disagreement_at_margin(frame: pd.DataFrame, margin: float,
                                    key: Sequence[str] = WINDOW_KEY,
                                    column: str = "delta_aic") -> pd.Series:
    """Pairwise disagreement among the analyses that do not abstain.

    The counterpart to :func:`pairwise_disagreement` for a decision rule that
    declines to answer when the evidence is weak. Only combinations with
    ``|delta| > margin`` speak; among those, this is the probability that two
    drawn at random reach opposite verdicts.

    Why it is needed. Comparing a remedy across designs of different sizes runs
    into the same trap as the headline: "some pair disagrees even after
    abstaining" gets easier the more cells there are, so a margin can look as
    though it stopped working when the only thing that changed was the
    enumeration. This statistic is a function of the proportions among the
    speaking analyses and so is comparable across designs.

    Windows where nothing speaks are dropped rather than scored zero — a rule
    that is silent is not a rule that agrees with itself, and counting silence
    as agreement is how a remedy flatters itself.
    """
    key = list(key)
    delta = frame[column].to_numpy(dtype=float)
    speaks = np.abs(delta) > margin
    sub = frame.loc[speaks, key].copy()
    sub["_climate"] = (delta[speaks] < 0).astype(float)
    grouped = sub.groupby(key)["_climate"].agg(["sum", "size"])
    grouped = grouped[grouped["size"] > 0]
    p = grouped["sum"] / grouped["size"]
    return 2.0 * p * (1.0 - p)


def variance_decomposition(frame: pd.DataFrame, column: str = "climate_wins",
                           key: Sequence[str] = WINDOW_KEY) -> dict:
    """Split the variation in a result into "between outbreaks" and "within".

    The question this answers is not how often analyses disagree but *where the
    disagreement in a literature comes from*. Reviews of climate--dengue studies
    routinely explain their conflicting findings by appealing to local context —
    this setting is more variable, that one denser, another poorer. That
    explanation requires the variation to sit between places. Here it can be
    checked, because every outbreak has been analysed every way:

        between = variance of the per-outbreak means
        within  = mean of the per-outbreak variances
        total   = variance over all fits

    ``between`` is the share attributable to everything that differs between
    epidemics — climate, health system, serotype, population, surveillance — and
    is therefore a ceiling on what any moderator analysis can explain.

    ``within`` is *not* the share attributable to the analysis, and an earlier
    version of this docstring said it was. It contains two different things: a
    main effect of the analysis, which standardising conventions would remove,
    and an outbreak-by-analysis interaction, which it would not. On this study's
    six-factor design the interaction is nearly four times the main effect, so
    the two readings recommend different remedies. Use
    :func:`two_way_decomposition` when the distinction matters.

    Unlike the share of outbreaks in which some pair disagrees, this quantity
    does not inflate as the factorial grows: adding analyses adds to both the
    numerator and the denominator. Measured on the 24-, 48- and 144-combination
    designs of this study it reads 78.3%, 79.7% and 76.5%.
    """
    key = list(key)
    y = frame[column].astype(float)
    total = float(y.var(ddof=0))
    if total <= 0:
        return dict(total=total, between=0.0, within=0.0,
                    between_share=float("nan"), within_share=float("nan"))
    between = float(frame.groupby(key)[column].mean().var(ddof=0))
    within = float(frame.groupby(key)[column]
                   .apply(lambda s: s.astype(float).var(ddof=0)).mean())
    return dict(total=total, between=between, within=within,
                between_share=between / total, within_share=within / total)


def benjamini_hochberg(pvalues) -> np.ndarray:
    """Benjamini-Hochberg q-values, in the order the p-values were given.

    Written out here because the obvious one-liner is wrong in a way that does
    not look wrong. Sorting the p-values ascending and taking a *forward*
    running minimum of ``p * m / rank`` drives every q-value down to the
    smallest one in the family, so everything "survives" — which is what a first
    draft of step 37 reported: 5 correlations below p = 0.05 and 24 of 24
    surviving correction. The monotonicity constraint runs the other way: q for
    the i-th smallest p is the minimum over all *larger* p-values, so the
    running minimum is taken from the largest downward.

    A study whose subject is undisclosed analytical freedom cannot afford a
    multiple-comparison correction that silently passes everything.
    """
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    if m == 0:
        return np.empty(0)
    order = np.argsort(p)
    ranks = np.arange(1, m + 1)
    q_ascending = p[order] * m / ranks
    q_ascending = np.minimum.accumulate(q_ascending[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q_ascending, 1.0)
    return out


def two_way_decomposition(frame: pd.DataFrame, cells: Sequence[str],
                          key: Sequence[str] = WINDOW_KEY,
                          column: str = "climate_wins") -> dict:
    """Split the variation three ways: outbreak, analysis, and the interaction.

    :func:`variance_decomposition` answers "is it the place?" and returns a
    remainder. The remainder was then read as though it were the analysis, which
    it is not. The design here is completely crossed and complete by
    construction — every retained outbreak was fitted under every combination —
    so the sum of squares partitions exactly:

        outbreak    variance of the per-outbreak means
        analysis    variance of the per-analysis means
        interaction what is left, which is neither

    The three shares sum to one. The distinction is not academic. A main effect
    of the analysis says some conventions are systematically more credulous than
    others, which stating and standardising the choices would fix. An interaction
    says *which* analysis endorses depends on *which* outbreak, so no fixed
    convention is the right one and the only honest report is the spread within
    each outbreak. On this study's data the interaction is the largest of the
    three, which is why the recommendation that survives is the per-outbreak
    interval rather than a house style.

    A specification curve shows the analysis main effect and nothing else: it
    orders specifications by their average result, which averages the interaction
    away. Here that hidden term is about four times the one on display.

    ``cells`` names the columns identifying an analysis. Passing an incomplete
    or unbalanced design gives shares that no longer sum to one, so the caller
    should restrict to complete windows first; the returned ``check`` reports the
    sum for exactly this reason.
    """
    key, cells = list(key), list(cells)
    y = frame[column].astype(float)
    if len(y) == 0 or float(y.var(ddof=0)) <= 0:
        return dict(outbreak=float("nan"), analysis=float("nan"),
                    interaction=float("nan"), check=float("nan"))
    grand = float(y.mean())
    a = frame.groupby(key)[column].transform("mean").astype(float)
    b = frame.groupby(cells)[column].transform("mean").astype(float)
    ss_outbreak = float(((a - grand) ** 2).sum())
    ss_analysis = float(((b - grand) ** 2).sum())
    ss_interaction = float(((y - a - b + grand) ** 2).sum())
    ss_total = float(((y - grand) ** 2).sum())
    return dict(outbreak=ss_outbreak / ss_total,
                analysis=ss_analysis / ss_total,
                interaction=ss_interaction / ss_total,
                check=(ss_outbreak + ss_analysis + ss_interaction) / ss_total)
