"""Tab 2 (registration) - qc'ed cell count, correlation distribution, and the
worst/best-aligning example images.

Registration correlation is measured per-cell, per registered channel pair, as
columns matching `qc_registration_pattern` (default "registration_corr") in
measure_intensity's object_features output. A cell counts as "qc'ed" only if every
matching column clears `qc_regcor`.

Example images are chosen by each field's own alignment rate - the percentage of its
cells clearing that same threshold - taking the `qc_n_sample_registration` worst and
best. run_registration.py already writes a PNG for every field (when rg_plot=true), so
this only picks which of them to embed; nothing is re-rendered.
"""

import glob
import logging
import os
import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from tglow.qc.assets import histogram_bar, image_to_data_uri, style_plot
from tglow.qc.registration import filter_registration_correlation, registration_correlation_columns

log = logging.getLogger(__name__)

# run_registration.py (see processes/registration.nf) writes these under
# <registration_images_dir>/<ref_plate>/<row letter>/<col>/, named
# "<field>_<qry_plate>_refch<N>_qrych<N>.png" - plate/row/col live in the
# directory structure rather than the filename, so the caption below is
# assembled from both.
FILENAME_RE = re.compile(r"^(?P<field>\d+)_(?P<qry_plate>.+)_refch(?P<ref_ch>\d+)_qrych(?P<qry_ch>\d+)\.png$")


def build_registration_stats(object_features, pattern, threshold):
    """Qc'ed cell count/percentage using the same filter as scaling-factor estimation."""
    corr_cols = registration_correlation_columns(object_features, pattern)
    if not corr_cols:
        return {"available": False, "corr_cols": []}

    qced = filter_registration_correlation(object_features, pattern, threshold)
    n_total = len(object_features)
    n_qced = len(qced)

    return {
        "available": True,
        "corr_cols": corr_cols,
        "n_total_cells": n_total,
        "n_qced_cells": n_qced,
        "pct_qced": round(100 * n_qced / n_total, 1) if n_total else 0,
    }


def build_correlation_density_plot(object_features, pattern, threshold):
    """Density histogram of each registration-correlation column, with a line at `threshold`.

    Binned here rather than in the browser (see assets.histogram_bar), which also
    means the overlaid traces need one shared set of bin edges - Plotly's own
    histogram traces negotiated that client-side.
    """
    corr_cols = registration_correlation_columns(object_features, pattern)
    if not corr_cols:
        return None

    values_by_col = {col: object_features[col].dropna().to_numpy() for col in corr_cols}
    all_values = np.concatenate([v for v in values_by_col.values() if v.size]) if any(v.size for v in values_by_col.values()) else np.array([])
    edges = np.histogram_bin_edges(all_values, bins=50) if all_values.size else 50

    fig = go.Figure()
    data_min, data_max = None, None
    for col, values in values_by_col.items():
        fig.add_trace(histogram_bar(values, bins=edges, density=True, name=col, opacity=0.6))
        if values.size:
            data_min = values.min() if data_min is None else min(data_min, values.min())
            data_max = values.max() if data_max is None else max(data_max, values.max())

    fig.add_vline(x=threshold, line_dash="dash", line_color="red",
                   annotation_text=f"qc_regcor={threshold}", annotation_position="top right")
    fig.update_layout(
        title="Registration correlation distribution",
        xaxis_title="Correlation",
        yaxis_title="Density",
        barmode="overlay",
    )
    if data_min is not None:
        # Correlation is mathematically bounded to [-1, 1] - clamp the axis range to
        # the real data span but never let it extend past that, even if the data
        # creeps slightly outside via floating point.
        fig.update_layout(xaxis=dict(range=[max(-1, data_min), min(1, data_max)]))

    # Lives in the sidebar (see the template) rather than tab-main, so it must not be
    # forced square/fixed-width - and a bit shorter than the default fits the sidebar
    # better than a plot sized for the main column. Tight margins so the plot fills
    # the sidebar's own narrow width instead of leaving Plotly's default ~80px
    # margins either side.
    style_plot(fig, square=False)
    fig.update_layout(height=300, margin=dict(l=45, r=15, t=32, b=35))
    return fig


def build_field_alignment(object_features, pattern, threshold, min_cells):
    """Per-field alignment rate: plate, well, field, n_cells, n_aligned, pct_aligned.

    "Aligned" is the same per-cell rule the rest of this tab uses -
    filter_registration_correlation, i.e. every matching correlation column >= threshold -
    so the percentage agrees with the qc_regcor shown beside it.

    Deliberately not image_features' n_cells_reg_corr_pass: nothing passes
    --registration-threshold to measure_intensity, so that column is computed at its
    hard-coded 0.4 default with a strict >, which would disagree with the threshold this
    tab reports. It also carries no per-field total to divide by.

    Fields with fewer than min_cells segmented cells are dropped - a field with 3 cells,
    one of which aligns, scores 33% and would otherwise crowd the worst-aligning list
    with sparse fields rather than genuine misalignment.
    """
    if not registration_correlation_columns(object_features, pattern):
        return pd.DataFrame(columns=["plate", "well", "field", "n_cells", "n_aligned", "pct_aligned"])

    keys = ["plate", "well", "field"]
    n_cells = object_features.groupby(keys).size().rename("n_cells")
    aligned = filter_registration_correlation(object_features, pattern, threshold).groupby(keys).size()

    df = n_cells.to_frame()
    # reindex, not a join - a field where nothing aligned has no row on the right and
    # must score 0 rather than drop out of the ranking entirely.
    df["n_aligned"] = aligned.reindex(df.index).fillna(0).astype(int)
    df["pct_aligned"] = 100 * df["n_aligned"] / df["n_cells"]

    return df[df["n_cells"] >= min_cells].reset_index()


def index_registration_images(registration_images_dir):
    """Map (plate, well, field) -> [png paths] over run_registration's output layout.

    run_registration.py writes <dir>/<plate>/<row letter>/<col>/<field>_<qry_plate>_refch<N>_qrych<N>.png,
    so plate/row/col come from the path and field/qry_plate from the filename. A field maps
    to a list because it gets one PNG per query plate (only one merged cycle is supported
    today, but the layout allows more).
    """
    if registration_images_dir is None or not os.path.isdir(registration_images_dir):
        return {}

    index = {}
    for path in sorted(glob.glob(os.path.join(registration_images_dir, "**", "*.png"), recursive=True)):
        match = FILENAME_RE.match(os.path.basename(path))
        if not match:
            log.warning(f"Registration image with unexpected name: {path}")
            continue

        plate, row_letter, col = os.path.normpath(path).split(os.sep)[-4:-1]
        well = f"{row_letter}{col.zfill(2)}"
        # object_features' field is an int; the filename yields a string
        index.setdefault((plate, well, int(match.group("field"))), []).append((path, match.group("qry_plate")))

    return index


def _samples_for(rows, image_index):
    """Turn ranked alignment rows into [{caption, data_uri}], skipping fields with no PNG."""
    images = []

    for row in rows.itertuples():
        for path, qry_plate in image_index.get((row.plate, row.well, int(row.field)), []):
            data_uri = image_to_data_uri(path)
            if data_uri is None:
                continue

            images.append({
                "caption": f"{row.plate} / {row.well} / field {row.field} (ref vs {qry_plate}) - "
                           f"{row.pct_aligned:.1f}% aligned ({row.n_cells} cells)",
                "data_uri": data_uri,
            })

    return images


def select_registration_samples(object_features, pattern, threshold, registration_images_dir, n_samples, min_cells):
    """The n_samples worst- and best-aligning fields, as {"worst": [...], "best": [...]}.

    Ranked by build_field_alignment's pct_aligned. Worst is listed worst-first and best
    best-first, so each gallery opens on its most extreme example. Random sampling (what
    this replaced) almost never surfaced a genuine registration failure; the extremes are
    the whole reason to look at these images.
    """
    alignment = build_field_alignment(object_features, pattern, threshold, min_cells)
    if alignment.empty:
        return {"worst": [], "best": []}

    image_index = index_registration_images(registration_images_dir)
    if not image_index:
        return {"worst": [], "best": []}

    ranked = alignment.sort_values(["pct_aligned", "plate", "well", "field"]).reset_index(drop=True)

    # With fewer than 2*n_samples eligible fields, head(n) and tail(n) would overlap and the
    # same field would appear in both galleries - split down the middle instead so the two
    # stay disjoint.
    split = n_samples if len(ranked) >= 2 * n_samples else len(ranked) // 2

    worst = ranked.head(split)
    best = ranked.iloc[split:][::-1].head(n_samples)

    return {"worst": _samples_for(worst, image_index), "best": _samples_for(best, image_index)}
