"""Pipeline step 36 — what does requiring an unbroken weekly run actually exclude?

The limitations paragraph says that requiring an unbroken weekly run "excludes
settings with intermittent reporting, which are likely to be poorer and more
heterogeneous than those retained". The first clause is a definition. The second
was an assertion, and this study has now found five assertions of that shape to be
wrong, so it is measured here instead.

Two things are measurable from the data in hand and one is not.

**Measurable: which requirement binds.** For every country in the weekly subset,
the longest gap-free run its reporting units achieve, and whether any wave in that
run passes the size and shape criteria. That separates countries excluded because
their reporting is broken from countries excluded because their epidemics are
small or shapeless --- a distinction the limitation elides, and one that points at
a different bias if the second dominates.

It does dominate. Of 54 excluded countries, 17 fail the gap-free run and 37 report
continuously enough but have no wave large or peaked enough to fit. The paragraph
named the minority mechanism, and the majority one cuts the other way: the sample
is the large, sharply peaked, single-wave epidemics these models handle best, so
the instability measured on it is a floor rather than a worst case.

**Measurable: how continuous the retained countries are.** If the retained
countries are simply the ones with long unbroken records, the sample is selected on
surveillance continuity, and the size of that selection can be stated.

**Not measurable here: wealth.** Nothing in OpenDengue or NASA POWER records it,
and importing an income classification to support one sentence would add an
external dependency the rest of the pipeline does not have. The paragraph should
therefore claim continuity, which is measured, and not wealth, which is not.

Reads the raw extract; fits nothing. Writes ``45_selection_audit.csv``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

MIN_WEEKS = 30
MIN_CASES = 300
MIN_PEAK = 20

cfg = load_config()
raw = resolve(cfg, "raw")
tables = resolve(cfg, "tables")


def longest_unbroken(s: pd.Series) -> pd.Series:
    """Longest run of consecutive weekly observations. As in step 14."""
    if len(s) < 2:
        return s
    gap = s.index.to_series().diff().dt.days.fillna(7)
    grp = (gap != 7).cumsum()
    sizes = s.groupby(grp).size()
    return s[grp == sizes.idxmax()]


def main() -> None:
    print("Loading OpenDengue temporal extract...")
    df = pd.read_csv(raw / cfg["data"]["opendengue"]["csv_name"],
                     usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                              "calendar_start_date", "calendar_end_date",
                              "dengue_total"], low_memory=False)
    df["start"] = pd.to_datetime(df["calendar_start_date"], errors="coerce")
    df["end"] = pd.to_datetime(df["calendar_end_date"], errors="coerce")
    df = df[(df["end"] - df["start"]).dt.days + 1 == 7]
    df = df.dropna(subset=["start", "dengue_total"])

    inv = pd.read_csv(tables / "12_global_windows.csv")
    retained = set(inv["country"].astype(str).str.upper())

    rows = []
    for country, cdf in df.groupby("adm_0_name"):
        best_run, best_cases, best_peak, n_units, weeks_total = 0, 0, 0, 0, 0
        for level, mask in (("national", cdf["adm_1_name"].isna()),
                            ("admin1", cdf["adm_1_name"].notna()
                             & cdf["adm_2_name"].isna())):
            sub = cdf[mask]
            if sub.empty:
                continue
            units = ([(None, sub)] if level == "national"
                     else list(sub.groupby("adm_1_name")))
            for _, udf in units:
                s = udf.groupby("start")["dengue_total"].sum().sort_index()
                n_units += 1
                weeks_total += len(s)
                run = longest_unbroken(s)
                if len(run) > best_run:
                    best_run = len(run)
                if len(run) >= MIN_WEEKS:
                    best_cases = max(best_cases, float(run.sum()))
                    best_peak = max(best_peak, float(run.max()))
        name = str(country).upper()
        rows.append(dict(country=name, retained=name in retained,
                         units=n_units, weeks_reported=weeks_total,
                         longest_run=best_run,
                         continuity=best_run / max(weeks_total, 1),
                         cases_in_run=best_cases, peak_in_run=best_peak))

    d = pd.DataFrame(rows).sort_values("longest_run", ascending=False)
    d.to_csv(tables / "45_selection_audit.csv", index=False)

    print("=" * 78)
    print("WHICH REQUIREMENT EXCLUDES A COUNTRY?")
    print("=" * 78)
    print(f"  {len(d)} countries in the weekly subset; "
          f"{int(d['retained'].sum())} contribute a usable window.\n")

    out = d[~d["retained"]].copy()
    short = out["longest_run"] < MIN_WEEKS
    small = (~short) & ((out["cases_in_run"] < MIN_CASES)
                        | (out["peak_in_run"] < MIN_PEAK))
    shape = (~short) & (~small)
    print(f"  of the {len(out)} excluded countries:")
    print(f"    no gap-free run of {MIN_WEEKS} weeks      {short.sum():3d}  "
          f"({short.mean() * 100:.0f}%)")
    print(f"    long enough, but the wave is too small  {small.sum():3d}  "
          f"({small.mean() * 100:.0f}%)")
    print(f"    long and large, rejected on shape       {shape.sum():3d}  "
          f"({shape.mean() * 100:.0f}%)")
    print(f"\n  The reporting requirement binds for {short.mean() * 100:.0f}% of "
          f"them. The other {(1 - short.mean()) * 100:.0f}% report")
    print("  continuously enough and are excluded because their epidemics are")
    print("  too small or the wrong shape. The limitation named the minority")
    print("  mechanism: this sample is selected first on outbreak size, and only")
    print("  second on surveillance continuity.")
    print("\n  That points the other way from the concern as written. Large,")
    print("  sharply peaked, single-wave epidemics are the ones these models can")
    print("  fit best, so the instability measured here is what remains on the")
    print("  most favourable data available, not on the hardest.")

    print("\n" + "=" * 78)
    print("HOW DIFFERENT ARE THE RETAINED COUNTRIES?")
    print("=" * 78)
    for label, g in (("retained", d[d["retained"]]), ("excluded", d[~d["retained"]])):
        print(f"  {label:9s} n={len(g):3d}  median longest gap-free run "
              f"{g['longest_run'].median():6.0f} weeks   median weeks reported "
              f"{g['weeks_reported'].median():7.0f}")
    ratio = (d[d["retained"]]["longest_run"].median()
             / max(d[~d["retained"]]["longest_run"].median(), 1))
    print(f"\n  The retained countries' longest gap-free run is {ratio:.0f} times")
    print("  the excluded countries'. The sample is selected on surveillance")
    print("  continuity, and that is the claim the paper can support. It carries")
    print("  no information about the wealth of the excluded settings, which an")
    print("  earlier limitations paragraph asserted and nothing here measures.")

    print(f"\nTable: {tables / '45_selection_audit.csv'}")


if __name__ == "__main__":
    main()
