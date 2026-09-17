"""Tab 2 (registration) - qc'ed cell count, correlation distribution, sample images.

Registration correlation is measured per-cell, per registered channel pair, as
columns matching `qc_registration_pattern` (default "registration_corr") in
measure_intensity's object_features output. A cell counts as "qc'ed" only if every
matching column clears `qc_regcor`.
"""

import glob
import logging
import os
import random
import re

import plotly.graph_objects as go

from tglow.qc.assets import image_to_data_uri, style_plot
from tglow.qc.registration import filter_registration_correlation, registration_correlation_columns

log = logging.getLogger(__name__)

# Fixed seed so the sampled registration images are stable across re-runs with
# unchanged inputs (keeps -resume/report diffs meaningful instead of churning
# on every run).
SAMPLE_SEED = 42

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
    """Density histogram of each registration-correlation column, with a line at `threshold`."""
    corr_cols = registration_correlation_columns(object_features, pattern)
    if not corr_cols:
        return None

    fig = go.Figure()
    for col in corr_cols:
        fig.add_trace(go.Histogram(
            x=object_features[col].dropna(),
            histnorm="probability density",
            name=col,
            opacity=0.6,
        ))

    fig.add_vline(x=threshold, line_dash="dash", line_color="red",
                   annotation_text=f"qc_regcor={threshold}", annotation_position="top right")
    fig.update_layout(
        title="Registration correlation distribution",
        xaxis_title="Correlation",
        yaxis_title="Density",
        barmode="overlay",
    )
    # Lives in the sidebar (see the template) rather than tab-main, so it must not be
    # forced square/fixed-width - and a bit shorter than the default fits the sidebar
    # better than a plot sized for the main column.
    style_plot(fig, square=False)
    fig.update_layout(height=300)
    return fig


def sample_registration_images(registration_images_dir, n_samples):
    """Randomly sample up to n_samples registration PNGs (deterministic seed) as data URIs.

    Returns a list of {"caption": "<plate> / <well> / field <N> (ref vs <qry_plate>)",
    "data_uri": ...} dicts, in sampled order.
    """
    if registration_images_dir is None or not os.path.isdir(registration_images_dir):
        return []

    all_pngs = sorted(glob.glob(os.path.join(registration_images_dir, "**", "*.png"), recursive=True))
    if not all_pngs:
        return []

    sampled = random.Random(SAMPLE_SEED).sample(all_pngs, k=min(n_samples, len(all_pngs)))

    images = []
    for path in sampled:
        data_uri = image_to_data_uri(path)
        if data_uri is None:
            continue

        plate, row_letter, col = os.path.normpath(path).split(os.sep)[-4:-1]
        well = f"{row_letter}{col.zfill(2)}"

        match = FILENAME_RE.match(os.path.basename(path))
        if match:
            caption = f"{plate} / {well} / field {match.group('field')} (ref vs {match.group('qry_plate')})"
        else:
            log.warning(f"Registration image with unexpected name: {path}")
            caption = f"{plate} / {well} / {os.path.basename(path)}"

        images.append({"caption": caption, "data_uri": data_uri})

    return images
