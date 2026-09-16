"""Tab 7 (debris) - debris percentage vs threshold/mean ("mean/otsu") ratio per channel.

Debris statistics are merged directly into measure_intensity's image_features output
(one row per field) - see build_debris_statistics_wide in
measure_intensity_features_with_debris.py - so this tab reads straight off
MeasurementData.image_features, with no separate debris_statistics file to load.

An image "passes" for a given channel when its debris_percentage is at or below
qc_debris_max_pct AND its threshold_mean_ratio is at or above qc_debris_min_ratio -
both configurable, since what counts as an acceptable debris level/threshold
separation is experiment-specific.
"""

import logging
import re

import plotly.graph_objects as go

log = logging.getLogger(__name__)

CHANNEL_COLUMN_RE = re.compile(r"^ch(\d+)__debris_percentage$")


def available_channels(image_features):
    """Channels (1-indexed) that have both debris_percentage and threshold_mean_ratio."""
    channels = sorted(
        int(m.group(1)) for col in image_features.columns for m in [CHANNEL_COLUMN_RE.match(col)] if m
    )
    return [c for c in channels if f"ch{c}__threshold_mean_ratio" in image_features.columns]


def _channel_columns(image_features, channel):
    pct_col = f"ch{channel}__debris_percentage"
    ratio_col = f"ch{channel}__threshold_mean_ratio"
    return image_features[[pct_col, ratio_col]].dropna(), pct_col, ratio_col


def build_debris_scatter_plots(image_features, channels, debris_max_pct, ratio_min):
    """dict[channel] -> Plotly scatter of debris_percentage (x) vs threshold_mean_ratio (y), one point per field."""
    figures = {}

    for channel in channels:
        df, pct_col, ratio_col = _channel_columns(image_features, channel)

        fig = go.Figure(data=go.Scatter(
            x=df[pct_col], y=df[ratio_col], mode="markers",
            marker=dict(size=6, opacity=0.6),
        ))
        fig.add_vline(x=debris_max_pct, line_dash="dash", line_color="red",
                      annotation_text=f"max debris %={debris_max_pct}", annotation_position="top right")
        fig.add_hline(y=ratio_min, line_dash="dash", line_color="red",
                      annotation_text=f"min ratio={ratio_min}", annotation_position="bottom right")
        fig.update_layout(
            title=f"Channel {channel}: debris % vs threshold/mean ratio",
            xaxis_title="Debris percentage",
            yaxis_title="Threshold / mean ratio",
        )
        figures[channel] = fig

    return figures


def build_debris_pass_table(image_features, channels, debris_max_pct, ratio_min):
    """Per-channel rows: n_total, n_passed, pct_passed, mean debris % among passing images."""
    rows = []

    for channel in channels:
        df, pct_col, ratio_col = _channel_columns(image_features, channel)
        n_total = len(df)

        passed = df[(df[pct_col] <= debris_max_pct) & (df[ratio_col] >= ratio_min)]
        n_passed = len(passed)

        rows.append({
            "channel": channel,
            "n_total": n_total,
            "n_passed": n_passed,
            "pct_passed": round(100 * n_passed / n_total, 1) if n_total else 0,
            "mean_debris_pct_passed": round(passed[pct_col].mean(), 2) if n_passed else None,
        })

    return rows
