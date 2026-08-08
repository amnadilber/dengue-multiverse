"""Pipeline step 35 — test the sentences, not the numbers.

Every figure in this paper is checked against a stored table by
``tests/test_paper_consistency.py``, and has been since a manuscript drifted from
the pipeline twice. That suite cannot see the other half of a paper. Between two
correct numbers sits a sentence saying what they mean, and four review passes
found five such sentences wrong --- each one plausible, each one sitting between
figures that were themselves right, and none of them checkable by any test then
in the repository.

This step exists so that the interpretive claims have somewhere to be checked
too. It refits nothing; it reads the stored tables and asks whether the claims
built on them survive.

**1. Outbreak, analysis, and the interaction.** The paper reported a one-way
split --- three-quarters of the variation in the verdict lies within outbreaks
rather than between them --- and glossed it as "two conflicting published
findings are more likely to differ because the analysts differed than because
the places did". The gloss does not follow from the split, because the two terms
are not the same kind of quantity: "between" is the variance of per-outbreak
means, "within" the mean of per-outbreak variances. The design is completely
crossed, so the honest version is available exactly: outbreak main effect,
analysis main effect, and their interaction, partitioning the total sum of
squares with nothing left over.

**2. Does the discordant observation-model cell diagnose a false positive?**
Poisson endorses climate forcing where the negative binomial does not in a third
of real paired fits. An earlier draft read that as mostly false positives. The
simulation can be asked directly, since both truths were run through the same
factorial.

**3. Do warm starts understate the analytical share?** Every combination within
an outbreak is warm-started from one anchor fit, which the paper argued makes the
144 analyses more alike than 144 independent analysts, so that the reported
analytical share is if anything an underestimate. The optimiser check refitted
seven outbreaks from cold starts, so the direction is measurable rather than
arguable --- on seven outbreaks, which is worth stating alongside the answer.

Writes ``43_two_way_decomposition.csv`` and ``44_claim_checks.csv``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402
from dengue_pk.robustness import (FACTORIAL_TABLES, complete_windows,  # noqa: E402
                                  latest_factorial, two_way_decomposition)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

cfg = load_config()
tables = resolve(cfg, "tables")
KEY = ["country", "unit", "window_start"]
FACTORS = ["observation", "temp_form", "rain_lag", "train_frac", "structure",
           "params"]
N_BOOT = 300


def cells_in(frame: pd.DataFrame) -> list[str]:
    return [f for f in FACTORS if f in frame.columns]


def main() -> None:
    d, n_combos = complete_windows(pd.read_csv(latest_factorial(tables)), KEY)
    cells = cells_in(d)
    print(f"{len(d):,} fits, {d.groupby(KEY).ngroups} outbreaks, "
          f"{n_combos} combinations\n")

    # -----------------------------------------------------------------------
    print("=" * 78)
    print("1. OUTBREAK, ANALYSIS, AND THE INTERACTION")
    print("=" * 78)
    print("  The one-way split says three-quarters of the variation is not")
    print("  between outbreaks. It does not say the analysis explains it.\n")

    rows = []
    for name in FACTORIAL_TABLES[::-1]:
        path = tables / name
        if not path.exists():
            continue
        raw = pd.read_csv(path)
        sub, n = complete_windows(raw, KEY)
        r = two_way_decomposition(sub, cells_in(raw), KEY)
        rows.append(dict(design=name, n_combos=n,
                         outbreaks=sub.groupby(KEY).ngroups, **r))
        print(f"  {n:3d} analyses   outbreak {r['outbreak'] * 100:5.1f}%   "
              f"analysis {r['analysis'] * 100:5.1f}%   "
              f"interaction {r['interaction'] * 100:5.1f}%")
    print("\n  The interaction is the largest term and it is the one that does")
    print("  not move with the size of the design.")

    # An interaction this size is only interesting if it is not an artefact of
    # the three countries that dominate the sample.
    print("\n  leaving out each dominant country:")
    for c in d["country"].value_counts().head(3).index:
        r = two_way_decomposition(d[d["country"] != c], cells, KEY)
        print(f"    without {c:12s} outbreak {r['outbreak'] * 100:5.1f}%   "
              f"analysis {r['analysis'] * 100:5.1f}%   "
              f"interaction {r['interaction'] * 100:5.1f}%")

    rng = np.random.default_rng(cfg["seed"])
    groups = [g for _, g in d.groupby("country")]
    boot = []
    for _ in range(N_BOOT):
        pick = rng.integers(0, len(groups), len(groups))
        frames = []
        for j, i in enumerate(pick):
            g = groups[i].copy()
            g["country"] = f"boot{j}"
            frames.append(g)
        boot.append(two_way_decomposition(pd.concat(frames, ignore_index=True),
                                          cells, KEY))
    bo = pd.DataFrame(boot)
    print(f"\n  resampling whole countries, {N_BOOT} draws:")
    base = two_way_decomposition(d, cells, KEY)
    for k in ("outbreak", "analysis", "interaction"):
        lo, hi = np.percentile(bo[k], [2.5, 97.5])
        rows.append(dict(design="six-factor CI", term=k,
                         point=round(base[k], 4), lo=round(float(lo), 4),
                         hi=round(float(hi), 4)))
        print(f"    {k:12s} {base[k] * 100:5.1f}%  "
              f"[95% CI {lo * 100:.1f}–{hi * 100:.1f}]")

    pd.DataFrame(rows).to_csv(tables / "43_two_way_decomposition.csv",
                              index=False)

    # The pairwise form of the same fact, which is what the discarded gloss was
    # really claiming and the form a reader can weigh against a pair of papers.
    print("\n  the same fact as two published findings would meet it:")
    p = d.groupby(KEY)["climate_wins"].mean()
    q = d.groupby(cells)["climate_wins"].mean()
    r = float(d["climate_wins"].mean())
    same_outbreak = float((2 * p * (1 - p)).mean())
    same_analysis = float((2 * q * (1 - q)).mean())
    neither = 2 * r * (1 - r)
    print(f"    same outbreak, two analyses      {same_outbreak * 100:5.1f}%")
    print(f"    same analysis, two outbreaks     {same_analysis * 100:5.1f}%")
    print(f"    different outbreak and analysis  {neither * 100:5.1f}%")
    print("    Standardising the analysis removes "
          f"{(neither - same_analysis) * 100:.1f} points of disagreement;")
    print("    standardising the place removes "
          f"{(neither - same_outbreak) * 100:.1f}. The gloss had it backwards.")

    checks = [dict(claim="P(disagree | same outbreak)", value=round(same_outbreak, 4)),
              dict(claim="P(disagree | same analysis)", value=round(same_analysis, 4)),
              dict(claim="P(disagree | neither shared)", value=round(neither, 4))]

    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("2. IS THE POISSON-ONLY ENDORSEMENT A FALSE POSITIVE?")
    print("=" * 78)
    sim_path = tables / "28_false_positive_6factor.csv"
    other = [c for c in cells if c != "observation"]
    if sim_path.exists():
        sim = pd.read_csv(sim_path)
        print(f"  {'data':12s} {'both':>7} {'Poisson only':>13} {'NB only':>9}")
        for truth in ("no_effect", "climate"):
            g = sim[sim["truth"] == truth]
            w = g.pivot_table(index=KEY + other, columns="observation",
                              values="climate_wins", aggfunc="first").dropna()
            if not {"nb", "poisson"} <= set(w.columns):
                continue
            only_p = float(((w["poisson"] == 1) & (w["nb"] == 0)).mean())
            print(f"  {truth:12s} "
                  f"{float(((w['poisson'] == 1) & (w['nb'] == 1)).mean()) * 100:6.1f}% "
                  f"{only_p * 100:12.1f}% "
                  f"{float(((w['poisson'] == 0) & (w['nb'] == 1)).mean()) * 100:8.1f}%")
            checks.append(dict(claim=f"Poisson-only, {truth}",
                               value=round(only_p, 4)))

        w = d.pivot_table(index=KEY + other, columns="observation",
                          values="climate_wins", aggfunc="first").dropna()
        real_only_p = float(((w["poisson"] == 1) & (w["nb"] == 0)).mean())
        print(f"  {'real':12s} "
              f"{float(((w['poisson'] == 1) & (w['nb'] == 1)).mean()) * 100:6.1f}% "
              f"{real_only_p * 100:12.1f}% "
              f"{float(((w['poisson'] == 0) & (w['nb'] == 1)).mean()) * 100:8.1f}%")
        checks.append(dict(claim="Poisson-only, real", value=round(real_only_p, 4)))
        print("\n  The discordant cell is diagnostic: it is nearly three times as")
        print("  common where no effect exists. But the real rate sits nearer the")
        print("  arm that contains an effect, not the arm that does not, so the")
        print("  claim that most of the gap is false positives is not supported")
        print("  by this table. It is the claim we withdraw.")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("3. DO WARM STARTS UNDERSTATE THE ANALYTICAL SHARE?")
    print("=" * 78)
    opt_path = tables / "23_optimiser_check.csv"
    if opt_path.exists():
        o = pd.read_csv(opt_path)
        for arm in ("as_run", "thorough"):
            g = o[o["arm"] == arm]
            if g.empty:
                continue
            pp = g.groupby(KEY)["climate_wins"].mean()
            val = float((2 * pp * (1 - pp)).mean())
            checks.append(dict(claim=f"pairwise disagreement, {arm}",
                               value=round(val, 4)))
            print(f"  {arm:9s} {len(pp)} outbreaks, pairwise disagreement "
                  f"{val * 100:5.1f}%")
        print("\n  Cold restarts give slightly less disagreement, not more. The")
        print("  paper had asserted the opposite direction and called it")
        print("  unambiguous. On seven outbreaks this settles nothing; what it")
        print("  does settle is that the direction was not known.")

    pd.DataFrame(checks).to_csv(tables / "44_claim_checks.csv", index=False)
    print(f"\nTables: {tables / '43_two_way_decomposition.csv'}")
    print(f"        {tables / '44_claim_checks.csv'}")


if __name__ == "__main__":
    main()
