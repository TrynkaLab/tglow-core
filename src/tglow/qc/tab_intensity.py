"""Tab 5 ((unscaled) intensities) - per-plate/channel/feature well-mean heatmaps + distributions.

Computed on qc'ed cells only (registration correlation >= sc_registration_thresh, same filter as
Tab 2/tab_registration.py) from measure_intensity's unscaled object_features output.
Features are the min/q25/median/q75/mean/max per-channel stats measure_intensity writes as
ch<N>__<stat> columns. Median is the default selected feature (see DEFAULT_FEATURE) -
the template selects it explicitly rather than relying on dict order.
"""

import logging
import re

import plotly.graph_objects as go

from tglow.qc.assets import histogram_bar, style_plot
from tglow.qc.plate_layout import build_well_grid, style_heatmap_axes
from tglow.qc.registration import filter_registration_correlation

log = logging.getLogger(__name__)

# Display label -> measure_intensity stat suffix
FEATURES = {"min": "min", "q25": "q25", "median": "median", "q75": "q75", "mean": "mean", "max": "max"}
DEFAULT_FEATURE = "median"

CHANNEL_COLUMN_RE = re.compile(r"^ch(\d+)__(.+)$")


def qced_cells(object_features, pattern, threshold):
    """Filter object_features down to qc'ed cells for intensity stats."""
    return filter_registration_correlation(object_features, pattern, threshold)


def available_channels(object_features):
    """Channels (0-indexed, matching the ch<N>__ column convention) that have every FEATURES stat present."""
    channels = sorted({int(m.group(1)) for col in object_features.columns for m in [CHANNEL_COLUMN_RE.match(col)] if m})
    return [c for c in channels if all(f"ch{c}__{stat}" in object_features.columns for stat in FEATURES.values())]


def build_intensity_heatmaps(qced_df, channels, plate_formats):
    """dict[channel][feature_label][plate] -> Plotly heatmap of per-well mean(feature).

    `plate_formats` is the dict plate -> (n_rows, n_cols) from
    plate_layout.infer_plate_formats - computed once (from the full, unfiltered
    well population) and shared with every other heatmap in the report. Inferring
    independently from qced_df here would be wrong: registration filtering can
    remove every qc'ed cell from a plate's higher rows/columns, which would
    under-infer that plate's format compared to e.g. Tab 1's cells-per-well
    heatmap for the same plate.
    """
    heatmaps = {}

    for channel in channels:
        heatmaps[channel] = {}
        for label, stat in FEATURES.items():
            col = f"ch{channel}__{stat}"
            heatmaps[channel][label] = {}

            for plate, plate_df in qced_df.groupby("plate"):
                well_means = plate_df.groupby(["row", "col", "well"])[col].mean().reset_index()
                grid, row_labels, col_labels = build_well_grid(
                    well_means, row_col="row", col_col="col", value_col=col, agg="mean", plate_format=plate_formats[plate]
                )

                fig = go.Figure(data=go.Heatmap(
                    z=grid, x=col_labels, y=row_labels, colorscale="Viridis",
                    colorbar=dict(title=label),
                    hovertemplate="row %{y} col %{x}<br>" + label + ": %{z:.1f}<extra></extra>",
                ))
                fig.update_yaxes(autorange="reversed")
                fig.update_layout(title=f"Ch{channel} {label} intensity - plate {plate}", xaxis_title="Column", yaxis_title="Row")
                style_heatmap_axes(fig)
                style_plot(fig, square=False, white_bg=False)
                heatmaps[channel][label][plate] = fig

    return heatmaps


def build_intensity_distributions(qced_df, channels):
    """dict[channel][feature_label] -> Plotly histogram of the raw per-cell values (all plates/wells).

    Binned here rather than in the browser (see assets.histogram_bar) - one raw
    trace per channel x feature over every qc'ed cell is what makes the report
    enormous.
    """
    distributions = {}

    for channel in channels:
        distributions[channel] = {}
        for label, stat in FEATURES.items():
            col = f"ch{channel}__{stat}"
            fig = go.Figure(data=histogram_bar(qced_df[col], bins=50))
            fig.update_layout(title=f"Ch{channel} {label} intensity distribution (qc'ed cells)", xaxis_title=label, yaxis_title="Count")
            # Lives in the sidebar (see the template) rather than tab-main, so it must
            # not be forced square/fixed-width - and a bit shorter fits the sidebar
            # better than a plot sized for the main column. Tight margins so the plot
            # fills the sidebar's own narrow width instead of leaving Plotly's default
            # ~80px margins either side.
            style_plot(fig, square=False)
            fig.update_layout(height=300, margin=dict(l=45, r=15, t=32, b=35))
            distributions[channel][label] = fig

    return distributions
