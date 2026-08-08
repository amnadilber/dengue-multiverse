"""
Pipeline step 13 — how many usable outbreak windows exist worldwide?

The Pakistani analysis showed that three routine choices each reverse the
in-sample verdict on climate-driven transmission. On three windows that is an
anecdote. The question worth answering is how often it happens, and answering it
requires many outbreaks from many countries.

OpenDengue covers 100+ countries. This script asks how many of them yield a
window that can actually carry a fit, applying the same criteria the Pakistani
windows had to meet:

* an unbroken weekly series, because a gap displaces the alignment between cases
  and climate for every week after it;
* long enough to contain a whole wave rather than a fragment;
* a wave large enough that the counts constrain anything at all;
* a rise and fall, not a plateau or a monotone trend, since a model of one
  epidemic fitted to a series without one estimates nothing.

Nothing is fitted here. The output is the inventory that decides whether the
multi-outbreak study is feasible, and it is deliberately separate so that the
selection rule is visible and auditable rather than buried inside an analysis.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# `dengue_pk` must be imported before NumPy: see dengue_pk/_msvc_runtime.py.
from dengue_pk import load_config, resolve  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402



MIN_WEEKS = 30            # a wave plus quiet shoulders on both sides
MAX_WEEKS = 110           # longer than one season; two waves fitted as one
MIN_CASES = 300           # below this the counts carry little information
MIN_PEAK = 20             # a wave, not a trickle
PEAK_RATIO = 2.5          # a wave, not a plateau; see the sensitivity table printed below
SHOULDER_FRACTION = 0.25  # a wave ends where it falls to a quarter of its peak
SHOULDER_WEEKS = 8        # quiet weeks kept either side, to anchor start and end

cfg = load_config()
raw = resolve(cfg, "raw")
tables = resolve(cfg, "tables")
figures = resolve(cfg, "figures")

print("Loading OpenDengue temporal extract (about 500 MB)...")
df = pd.read_csv(raw / cfg["data"]["opendengue"]["csv_name"],
                 usecols=["adm_0_name", "adm_1_name", "adm_2_name",
                          "calendar_start_date", "calendar_end_date",
                          "dengue_total", "S_res"], low_memory=False)
df["start"] = pd.to_datetime(df["calendar_start_date"], errors="coerce")
df["end"] = pd.to_datetime(df["calendar_end_date"], errors="coerce")
df = df[(df["end"] - df["start"]).dt.days + 1 == 7]
df = df.dropna(subset=["start", "dengue_total"])
print(f"{len(df):,} weekly records across "
      f"{df['adm_0_name'].nunique()} countries\n")


def longest_unbroken(s: pd.Series):
    """Longest run of consecutive weekly observations."""
    if len(s) < 2:
        return s
    gap = s.index.to_series().diff().dt.days.fillna(7)
    grp = (gap != 7).cumsum()
    sizes = s.groupby(grp).size()
    return s[grp == sizes.idxmax()]


def waves_in(run: pd.Series):
    """Extract individual epidemic waves from an unbroken run.

    Cutting a long run at troughs is not enough. Where transmission never falls
    close to zero — Singapore, Colombia and Belize are endemic in this dataset —
    no trough qualifies and the whole record survives as one segment, giving a
    595-week "wave" that a single-epidemic model cannot describe and should never
    be asked to.

    Waves are therefore built outward from peaks rather than inward from troughs.
    Each sufficiently large local maximum is taken in turn, its window extended
    in both directions until the counts fall below a fraction of that peak or a
    higher peak is met, and the result kept only if it is long enough to contain
    the wave and short enough to be one.
    """
    n = len(run)
    if n < MIN_WEEKS:
        return []

    # Boundaries are found on a three-week average. On the raw counts a single
    # quiet week inside a wave terminates it, splitting one epidemic into two
    # fragments that neither model can fit.
    smooth = run.rolling(3, center=True, min_periods=1).mean().to_numpy(float)
    claimed = np.zeros(n, dtype=bool)
    out = []

    for pk in np.argsort(smooth)[::-1]:
        if claimed[pk] or smooth[pk] < MIN_PEAK:
            continue
        threshold = SHOULDER_FRACTION * smooth[pk]
        half = MAX_WEEKS // 2

        lo = pk
        while (lo > 0 and not claimed[lo - 1] and pk - lo < half
               and smooth[lo - 1] > threshold):
            lo -= 1
        hi = pk
        while (hi < n - 1 and not claimed[hi + 1] and hi - pk < half
               and smooth[hi + 1] > threshold):
            hi += 1

        # Include the quiet shoulders on either side: the model needs weeks
        # before take-off and after decline to constrain the initial condition
        # and the final size.
        lo2 = lo
        while lo2 > 0 and not claimed[lo2 - 1] and lo - lo2 < SHOULDER_WEEKS:
            lo2 -= 1
        hi2 = hi
        while hi2 < n - 1 and not claimed[hi2 + 1] and hi2 - hi < SHOULDER_WEEKS:
            hi2 += 1

        claimed[lo2:hi2 + 1] = True
        if MIN_WEEKS <= hi2 - lo2 + 1 <= MAX_WEEKS:
            out.append(run.iloc[lo2:hi2 + 1])

    return sorted(out, key=lambda s: s.index.min())


def usable(seg: pd.Series) -> tuple[bool, str]:
    total, peak = seg.sum(), seg.max()
    if total < MIN_CASES:
        return False, f"only {total:.0f} cases"
    if peak < MIN_PEAK:
        return False, f"peak of {peak:.0f}"
    median = np.median(seg.to_numpy(float))
    if median > 0 and peak / median < PEAK_RATIO:
        return False, f"peak only {peak / median:.1f}x the median week"
    # The peak must be interior: a wave still rising at the end of the series
    # has no decline to constrain the final size.
    pk = int(np.argmax(seg.to_numpy()))
    if pk < 3 or pk > len(seg) - 4:
        return False, "peak at the edge of the series"
    return True, ""


segments = []
for country, cdf in df.groupby("adm_0_name"):
    for level, mask in (("national", cdf["adm_1_name"].isna()),
                        ("admin1", cdf["adm_1_name"].notna()
                         & cdf["adm_2_name"].isna())):
        sub = cdf[mask]
        if sub.empty:
            continue
        units = ([(None, sub)] if level == "national"
                 else list(sub.groupby("adm_1_name")))
        for unit, udf in units:
            s = udf.groupby("start")["dengue_total"].sum().sort_index()
            for seg in waves_in(longest_unbroken(s)):
                segments.append((str(country), level,
                                 str(unit) if unit else str(country), seg))

# How many windows survive at each value of the peak-prominence criterion?
# Printed before the value is applied, so the choice is visible rather than
# tuned to whatever yields a convenient sample size. Trimming a window to its
# wave necessarily raises its median week, so a threshold calibrated on
# untrimmed segments is too strict here — this table is what shows by how much.
print("--- windows surviving at each peak-to-median threshold ---")
for ratio in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
    n = sum(1 for _, _, _, seg in segments
            if seg.sum() >= MIN_CASES and seg.max() >= MIN_PEAK
            and (np.median(seg.to_numpy(float)) == 0
                 or seg.max() / np.median(seg.to_numpy(float)) >= ratio))
    print(f"  peak >= {ratio:.1f}x median week: {n:3d} windows")
print(f"  (applying {PEAK_RATIO:.1f})\n")

rows, rejected = [], []
for country, cdf in df.groupby("adm_0_name"):
    for level, mask in (("national", cdf["adm_1_name"].isna()),
                        ("admin1", cdf["adm_1_name"].notna()
                         & cdf["adm_2_name"].isna())):
        sub = cdf[mask]
        if sub.empty:
            continue
        units = ([(None, sub)] if level == "national"
                 else list(sub.groupby("adm_1_name")))
        for unit, udf in units:
            s = udf.groupby("start")["dengue_total"].sum().sort_index()
            run = longest_unbroken(s)
            for seg in waves_in(run):
                ok, why = usable(seg)
                rec = dict(country=str(country), level=level,
                           unit=str(unit) if unit else str(country),
                           start=seg.index.min().date(),
                           end=seg.index.max().date(),
                           weeks=len(seg), cases=int(seg.sum()),
                           peak=int(seg.max()),
                           peak_week=int(np.argmax(seg.to_numpy())))
                (rows if ok else rejected).append({**rec, "reason": why})

inv = pd.DataFrame(rows).sort_values(["country", "start"])
inv.to_csv(tables / "12_global_windows.csv", index=False)

print(f"{len(inv)} usable windows from {inv['country'].nunique()} countries")
print(f"({len(rejected)} candidate segments rejected)\n")

print("--- by country ---")
by_country = (inv.groupby("country")
              .agg(windows=("weeks", "size"), total_cases=("cases", "sum"),
                   median_weeks=("weeks", "median"))
              .sort_values("windows", ascending=False))
print(by_country.head(25).to_string())

print(f"\n--- spread ---")
print(f"weeks per window:  median {inv['weeks'].median():.0f}, "
      f"range {inv['weeks'].min()}–{inv['weeks'].max()}")
print(f"cases per window:  median {inv['cases'].median():,.0f}, "
      f"range {inv['cases'].min():,}–{inv['cases'].max():,}")
print(f"national windows:  {(inv['level'] == 'national').sum()}")
print(f"admin-1 windows:   {(inv['level'] == 'admin1').sum()}")

if len(rejected):
    rej = pd.DataFrame(rejected)
    print(f"\n--- why segments were rejected ---")
    print(rej["reason"].str.replace(r"[\d.,]+", "N", regex=True)
          .value_counts().head(8).to_string())

fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
axes[0].hist(inv["weeks"], bins=30, color="#1f6f8b")
axes[0].set_xlabel("weeks in window")
axes[0].set_ylabel("windows")
axes[0].set_title(f"{len(inv)} usable windows")
axes[0].grid(alpha=0.3)

axes[1].barh(by_country.head(20).index[::-1],
             by_country.head(20)["windows"][::-1], color="#8c1c13")
axes[1].set_xlabel("usable windows")
axes[1].set_title("Countries contributing most")
axes[1].grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(figures / "11_global_windows.png", dpi=150)
print(f"\nFigure: {figures / '11_global_windows.png'}")
print(f"Table:  {tables / '12_global_windows.csv'}")
