"""Tab 1 (general QC) data - counts, cells-per-well heatmaps.

Loaded once from measure_intensity's per-plate-staged parquet output and shared
by tab_registration.py/tab_intensity.py so the object_features/image_features
frames aren't re-read from disk for every tab.
"""

import logging

import pandas as pd
import plotly.graph_objects as go

from tglow.qc.io import load_measurements
from tglow.qc.plate_layout import build_well_grid, style_heatmap_axes

log = logging.getLogger(__name__)


class MeasurementData:
    """Holds the concatenated object_features/image_features frames for the whole run."""

    def __init__(self, measurements_dir):
        self.object_features = load_measurements(measurements_dir, "object_features")
        self.image_features = load_measurements(measurements_dir, "image_features")


def load_blacklist(blacklist_path, plates=None):
    """Load the <plate>\\t<well> blacklist TSV (no header). Returns a DataFrame, or empty if None."""
    if blacklist_path is None:
        return pd.DataFrame(columns=["plate", "well"])

    blacklist = pd.read_csv(blacklist_path, sep="\t", header=None, names=["plate", "well"], dtype=str)
    if plates is not None:
        blacklist = blacklist[blacklist["plate"].isin(plates)]
    return blacklist


def count_cycles(registration_manifest_path):
    """Number of cycles per registration group: 1 (reference) + len(query_plates).

    Returns None if no registration manifest was provided. If groups disagree on
    cycle count, returns the max and logs a warning (mixed-cycle-count runs are
    unusual but not invalid).
    """
    if registration_manifest_path is None:
        return None

    manifest = pd.read_csv(registration_manifest_path, sep="\t", dtype=str)
    cycle_counts = manifest["query_plates"].str.split(",").apply(len) + 1

    if cycle_counts.nunique() > 1:
        log.warning(f"Registration manifest groups disagree on cycle count: {sorted(cycle_counts.unique())} - reporting the max")

    return int(cycle_counts.max())


def images_without_cells(image_features, object_features):
    """Count of (plate, well, field) image rows with zero matching object rows."""
    key_cols = ["plate", "well", "field"]
    cells_per_image = object_features.groupby(key_cols).size().rename("n_cells")
    merged = image_features[key_cols].merge(cells_per_image, on=key_cols, how="left")
    return int(merged["n_cells"].isna().sum())


def fields_per_well_summary(image_features):
    """Summarize fields-per-well as a single display string (e.g. "5" or "3-5 (mean 4.2)")."""
    counts = image_features.groupby(["plate", "well"])["field"].nunique()
    if counts.nunique() == 1:
        return str(int(counts.iloc[0]))
    return f"{int(counts.min())}-{int(counts.max())} (mean {counts.mean():.1f})"


def build_general_stats(measurements, blacklist_df, registration_manifest_path):
    """Build the Tab 1 headline stat dict."""
    image_features = measurements.image_features
    object_features = measurements.object_features

    n_wells = image_features.drop_duplicates(["plate", "well"]).shape[0]
    n_cells = len(object_features)

    return {
        "n_images": len(image_features),
        "n_cells": n_cells,
        "fields_per_well": fields_per_well_summary(image_features),
        "n_plates": image_features["plate"].nunique(),
        "n_wells": n_wells,
        "n_blacklisted_wells": len(blacklist_df),
        "n_images_without_cells": images_without_cells(image_features, object_features),
        "n_cycles": count_cycles(registration_manifest_path),
        "mean_cells_per_well": round(n_cells / n_wells, 1) if n_wells else 0,
    }


def build_cells_per_well_heatmaps(object_features, plate_formats):
    """One Plotly heatmap figure per plate: cell count per well.

    `plate_formats` is the dict plate -> (n_rows, n_cols) from
    plate_layout.infer_plate_formats - computed once, shared with every other
    heatmap in the report, so a given plate always renders at the same format.
    """
    figures = {}

    cells_per_well = (
        object_features.groupby(["plate", "row", "col", "well"])
        .size()
        .rename("n_cells")
        .reset_index()
    )

    for plate, plate_df in cells_per_well.groupby("plate"):
        grid, row_labels, col_labels = build_well_grid(
            plate_df, row_col="row", col_col="col", value_col="n_cells", agg="sum", plate_format=plate_formats[plate]
        )

        fig = go.Figure(data=go.Heatmap(
            z=grid,
            x=col_labels,
            y=row_labels,
            colorscale="Viridis",
            colorbar=dict(title="cells"),
            hovertemplate="row %{y} col %{x}<br>cells: %{z}<extra></extra>",
        ))
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(title=f"Cells per well - plate {plate}", xaxis_title="Column", yaxis_title="Row")
        style_heatmap_axes(fig)
        figures[plate] = fig

    return figures
