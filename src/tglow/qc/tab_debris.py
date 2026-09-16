"""Tab 7 (debris) - debris percentage vs threshold/mean ("mean/otsu") ratio per channel,
per-plate/channel summary table, and clickable highest-debris example images.

Debris statistics are merged directly into measure_intensity's image_features output
(one row per field) - see build_debris_statistics_wide in
measure_intensity_features_with_debris.py - so this tab reads straight off
MeasurementData.image_features, with no separate debris_statistics file to load.

An image counts as "with debris" for a given channel when its debris_percentage is
above debris_max_pct OR its threshold_mean_ratio is below ratio_min - both
configurable, since what counts as an acceptable debris level/threshold separation is
experiment-specific. "Without debris" is the complement (both criteria clear).
"""

import glob
import logging
import os
import re

import plotly.graph_objects as go

from tglow.qc.assets import image_to_data_uri

log = logging.getLogger(__name__)

CHANNEL_COLUMN_RE = re.compile(r"^ch(\d+)__debris_percentage$")

DEBRIS_SAMPLE_RE = re.compile(
    r"^(?P<plate>.+)_(?P<well>[A-Za-z]+\d+)_(?P<field>\d+)_ch(?P<channel>\d+)_pct(?P<pct>[\d.]+)_debris\.png$"
)


def available_channels(image_features):
    """Channels (1-indexed) that have both debris_percentage and threshold_mean_ratio."""
    channels = sorted(
        int(m.group(1)) for col in image_features.columns for m in [CHANNEL_COLUMN_RE.match(col)] if m
    )
    return [c for c in channels if f"ch{c}__threshold_mean_ratio" in image_features.columns]


def _channel_columns(image_features, channel):
    pct_col = f"ch{channel}__debris_percentage"
    ratio_col = f"ch{channel}__threshold_mean_ratio"
    return pct_col, ratio_col


def build_debris_scatter_plots(image_features, channels, debris_max_pct, ratio_min):
    """dict[channel] -> Plotly scatter of debris_percentage (x) vs threshold_mean_ratio (y), one point per field."""
    figures = {}

    for channel in channels:
        pct_col, ratio_col = _channel_columns(image_features, channel)
        df = image_features[[pct_col, ratio_col]].dropna()

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


def build_debris_summary_table(image_features, channels, debris_max_pct, ratio_min):
    """One row per (plate, channel): n_total, n_with_debris/n_without_debris (+ pct), mean threshold/background.

    "With debris" fails the pass criteria (debris_percentage > debris_max_pct OR
    threshold_mean_ratio < ratio_min); "without debris" is the complement - the
    inverse framing of the pass/fail check itself, not a different threshold.
    """
    rows = []

    for plate in sorted(image_features["plate"].unique()):
        plate_df = image_features[image_features["plate"] == plate]

        for channel in channels:
            pct_col, ratio_col = _channel_columns(image_features, channel)
            threshold_col = f"ch{channel}__threshold"
            background_col = f"ch{channel}__mean_intensity"

            df = plate_df[[pct_col, ratio_col]].dropna()
            n_total = len(df)

            without_debris = df[(df[pct_col] <= debris_max_pct) & (df[ratio_col] >= ratio_min)]
            n_without = len(without_debris)
            n_with = n_total - n_without

            row = {
                "plate": plate,
                "channel": channel,
                "n_total": n_total,
                "n_with_debris": n_with,
                "pct_with_debris": round(100 * n_with / n_total, 1) if n_total else 0,
                "n_without_debris": n_without,
                "pct_without_debris": round(100 * n_without / n_total, 1) if n_total else 0,
                "mean_threshold": round(plate_df[threshold_col].mean(), 2) if threshold_col in plate_df.columns else None,
                "mean_background": round(plate_df[background_col].mean(), 2) if background_col in plate_df.columns else None,
            }
            rows.append(row)

    return rows


def build_debris_param_summary(image_features, debris_max_pct, ratio_min):
    """Side-table params: debris method/cell-mask expansion (constant per run, read
    straight off image_features) plus the configured pass-criteria thresholds - same
    "just pass a dict through" pattern as tab_flatfield.build_flatfield_param_summary."""
    method = image_features["method"].iloc[0] if "method" in image_features.columns and len(image_features) else "n/a"
    expansion = image_features["cellmask_expansion"].iloc[0] if "cellmask_expansion" in image_features.columns and len(image_features) else "n/a"

    return {
        "Debris method": method,
        "Cell mask expansion": expansion,
        "Max debris %": debris_max_pct,
        "Min threshold/mean ratio": ratio_min,
    }


def sample_debris_images(debris_samples_dir):
    """Group the highest-debris overlay PNGs written by make_debris_overlay.py by channel.

    Returns dict channel(int) -> list of {"caption", "data_uri"}, sorted by the
    embedded debris percentage descending (worst first).
    """
    if debris_samples_dir is None or not os.path.isdir(debris_samples_dir):
        return {}

    by_channel = {}
    for path in sorted(glob.glob(os.path.join(debris_samples_dir, "*_debris.png"))):
        match = DEBRIS_SAMPLE_RE.match(os.path.basename(path))
        if not match:
            log.warning(f"Skipping debris sample image with unexpected name: {path}")
            continue

        data_uri = image_to_data_uri(path)
        if data_uri is None:
            continue

        channel = int(match.group("channel"))
        pct = float(match.group("pct"))
        by_channel.setdefault(channel, []).append({
            "caption": f"{match.group('plate')} / {match.group('well')} / field {match.group('field')} (debris%={pct})",
            "data_uri": data_uri,
            "_pct": pct,
        })

    for channel, images in by_channel.items():
        images.sort(key=lambda img: img["_pct"], reverse=True)
        for img in images:
            del img["_pct"]

    return by_channel
