"""Well QC tab: the per-well QC table written by the pipeline's build_well_qc.py (well_qc.tsv).

The table is rendered as-is, so its columns follow whatever build_well_qc.py wrote (e.g. one
pct_low_ch<N> column per checked channel). The verdict counts also feed Tab 1.
"""

import pandas as pd

# Display order - the verdicts worth acting on first
VERDICT_ORDER = ["fail", "warn", "blacklisted", "pass"]

# Written as integers, but read back as float whenever a well has no value (e.g. not measured)
INT_COLUMNS = ["col", "n_cells"]


def load_well_qc(path):
    table = pd.read_csv(path, sep="\t", dtype={"plate": str, "well": str, "row": str})
    for col in INT_COLUMNS:
        if col in table.columns:
            table[col] = table[col].astype("Int64")
    for col in ("missing_cycles", "qc_flags_fail", "qc_flags_warn"):
        if col in table.columns:
            table[col] = table[col].fillna("")
    return table


def _format_cell(value, numeric):
    """(display text, sort key) for one cell. Blank values sort below every real number."""
    if pd.isna(value):
        return "", -1 if numeric else ""
    if numeric and isinstance(value, float):
        return f"{value:.2f}", value
    return str(value), value


def build_well_qc_tab(well_qc_path):
    table = load_well_qc(well_qc_path)

    rank = {v: i for i, v in enumerate(VERDICT_ORDER)}
    table = table.assign(_rank=table["qc_verdict"].map(rank).fillna(len(rank))) \
        .sort_values(["_rank", "plate", "row", "col"]).drop(columns="_rank")

    numeric = {col: pd.api.types.is_numeric_dtype(table[col]) for col in table.columns}
    rows = [
        {
            "verdict": record["qc_verdict"],
            "cells": [_format_cell(record[col], numeric[col]) for col in table.columns],
        }
        for record in table.to_dict("records")
    ]

    counts = table["qc_verdict"].value_counts()
    return {
        "available": True,
        "columns": [{"name": col, "numeric": numeric[col]} for col in table.columns],
        "rows": rows,
        "n_wells": len(table),
        "counts": {v: int(counts.get(v, 0)) for v in VERDICT_ORDER},
    }
