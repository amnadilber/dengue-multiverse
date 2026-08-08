"""
Output integrity.

Results leave this project as CSV files, and a value that does not survive the
round trip is worse than a missing one: it is silently absent from every later
analysis while the file still looks complete.

This is not hypothetical. The constant-transmission model was originally labelled
``null``, which pandas parses back as ``NaN`` by default. Every one of its rows
vanished from the sensitivity table on re-import — 37 fits, several hours of
computation, present in the file and invisible to any reader who opened it the
obvious way. The label is now ``constant``; these tests keep it that way.
"""

import csv
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"

# Strings pandas converts to NaN unless explicitly told otherwise. Any of these
# used as a category label disappears on re-import.
PANDAS_NA_STRINGS = {
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
    "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a",
    "nan", "null", "", "NAN",
}


def table_files():
    return sorted(TABLES.glob("*.csv")) if TABLES.exists() else []


@pytest.mark.parametrize("path", table_files(), ids=lambda p: p.name)
def test_no_label_is_swallowed_by_the_csv_reader(path):
    """No cell in a result table may be a string pandas reads back as missing.

    Checked against the raw text rather than the parsed frame, because by the
    time pandas has parsed it the evidence is gone.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        pytest.skip("empty table")

    header, *body = rows
    offenders = {}
    for row in body:
        for column, cell in zip(header, row):
            if cell in PANDAS_NA_STRINGS and cell != "":
                offenders.setdefault(column, set()).add(cell)

    assert not offenders, (
        f"{path.name} contains labels that pandas reads as NaN: {offenders}. "
        f"Rename them — a category that disappears on re-import is worse than a "
        f"missing one, because the file still looks complete.")


@pytest.mark.parametrize("path", table_files(), ids=lambda p: p.name)
def test_tables_round_trip_without_losing_rows(path):
    """Reading a table back must preserve every row and every label column."""
    with open(path, newline="", encoding="utf-8") as fh:
        raw_rows = sum(1 for _ in csv.reader(fh)) - 1
    if raw_rows <= 0:
        pytest.skip("empty table")

    df = pd.read_csv(path)
    assert len(df) == raw_rows, (
        f"{path.name}: {raw_rows} data rows on disk, {len(df)} after parsing")

    if "model" in df.columns:
        assert df["model"].notna().all(), (
            f"{path.name}: the 'model' column has missing values after parsing, "
            f"which usually means a label collided with a pandas NA string")


def test_model_labels_are_the_expected_ones():
    """Guard the specific rename, so it cannot regress silently."""
    comparison = TABLES / "02_model_comparison.csv"
    if not comparison.exists():
        pytest.skip("run scripts/03_fit_classical.py first")
    models = set(pd.read_csv(comparison)["model"].dropna())
    assert models == {"climate", "constant"}, f"unexpected model labels: {models}"
