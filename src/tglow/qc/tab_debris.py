"""Tab 7 (debris) - debris percentage vs threshold/mean ("mean/otsu") ratio per channel,
per-plate/channel summary table, and clickable example image viewers for the
"no_debris" (highest-ratio), "uncertain" (highest-debris) and "debris"
(highest-debris-that-still-passed) classes.

These three names are the single vocabulary for the classification: they key the
sample filenames select_debris_samples.py/make_debris_overlay.py write, the
summary_table fields, render.py's sample_images_* keys and the template's carousel
ids, and they are what the report displays (underscores shown as spaces). Renaming
one means renaming all of them.

Debris statistics are merged directly into measure_intensity's image_features output
(one row per field) - see build_debris_statistics_wide in
measure_intensity_features_with_debris.py - so this tab reads straight off
MeasurementData.image_features, with no separate debris_statistics file to load.

Every image is classified per channel into exactly one of three groups (both
ratio_min/debris_max_pct configurable, since what counts as an acceptable debris
level/threshold separation is experiment-specific):
- "no_debris" (threshold_mean_ratio < ratio_min): the threshold sits too close to
  the background mean to be a meaningful separation, which suggests there is no
  debris to find - so debris_percentage computed from it isn't meaningful either way.
- "debris" (threshold_mean_ratio >= ratio_min AND debris_percentage <
  debris_max_pct): a reliably measured, acceptable debris level - the normal,
  usable case.
- "uncertain" (threshold_mean_ratio >= ratio_min AND debris_percentage >=
  debris_max_pct): a reliable threshold but an unusually high debris fraction -
  often a sign the background itself is anomalously uniform/dark (making the
  debris mask spuriously large) rather than a genuine debris detection, so it's
  flagged as suspect rather than trusted outright.
"""

import glob
import logging
import os
import re

import plotly.graph_objects as go

from tglow.qc.assets import image_to_data_uri, style_plot

log = logging.getLogger(__name__)

CHANNEL_COLUMN_RE = re.compile(r"^ch(\d+)__debris_percentage$")

DEBRIS_SAMPLE_RE = re.compile(
    r"^(?P<plate>.+)_(?P<well>[A-Za-z]+\d+)_(?P<field>\d+)_ch(?P<channel>\d+)_(?P<sample_class>no_debris|uncertain|debris)_pct(?P<pct>[\d.]+)_ratio(?P<ratio>[\d.]+)_debris\.png$"
)

# "no_debris" samples are picked (and should be displayed) by highest ratio, not
# highest debris % - the whole point of that class is that its debris % isn't
# meaningful - so it needs its own sort key, unlike the other two classes.
SAMPLE_CLASS_SORT_KEY = {"no_debris": "_ratio"}
DEFAULT_SAMPLE_SORT_KEY = "_pct"


def available_channels(image_features):
    """Channels (0-indexed) that have both debris_percentage and threshold_mean_ratio."""
    channels = sorted(
        int(m.group(1)) for col in image_features.columns for m in [CHANNEL_COLUMN_RE.match(col)] if m
    )
    return [c for c in channels if f"ch{c}__threshold_mean_ratio" in image_features.columns]


def _channel_columns(image_features, channel):
    pct_col = f"ch{channel}__debris_percentage"
    ratio_col = f"ch{channel}__threshold_mean_ratio"
    return pct_col, ratio_col


def build_debris_scatter_plots(image_features, channels, debris_max_pct, ratio_min):
    """dict[channel] -> Plotly scatter of threshold_mean_ratio (x) vs debris_percentage (y), one point per field."""
    figures = {}

    for channel in channels:
        pct_col, ratio_col = _channel_columns(image_features, channel)
        df = image_features[[pct_col, ratio_col]].dropna()

        fig = go.Figure(data=go.Scatter(
            x=df[ratio_col], y=df[pct_col], mode="markers",
            marker=dict(size=6, opacity=0.6),
        ))
        fig.add_hline(y=debris_max_pct, line_dash="dash", line_color="red",
                      annotation_text=f"max debris %={debris_max_pct}", annotation_position="top right")
        fig.add_vline(x=ratio_min, line_dash="dash", line_color="red",
                      annotation_text=f"min ratio={ratio_min}", annotation_position="bottom right")
        fig.update_layout(
            title=f"Channel {channel}: threshold/mean ratio vs debris %",
            xaxis_title="Threshold / mean ratio",
            yaxis_title="Debris percentage",
        )
        style_plot(fig)
        figures[channel] = fig

    return figures


def build_debris_summary_table(image_features, channels, debris_max_pct, ratio_min):
    """One row per channel (across all plates): n_total, the 3-way classification counts (+ pct), mean debris %/threshold/background.

    See the module docstring for the no_debris/debris/uncertain
    definitions - every image with both stats present falls into exactly one.
    """
    rows = []

    for channel in channels:
        pct_col, ratio_col = _channel_columns(image_features, channel)
        threshold_col = f"ch{channel}__threshold"
        background_col = f"ch{channel}__mean_intensity"

        df = image_features[[pct_col, ratio_col]].dropna()
        n_total = len(df)

        n_no_debris = len(df[df[ratio_col] < ratio_min])
        n_debris = len(df[(df[ratio_col] >= ratio_min) & (df[pct_col] < debris_max_pct)])
        n_uncertain = len(df[(df[ratio_col] >= ratio_min) & (df[pct_col] >= debris_max_pct)])

        def _pct(n):
            return round(100 * n / n_total, 1) if n_total else 0

        row = {
            "channel": channel,
            "n_total": n_total,
            "n_no_debris": n_no_debris,
            "pct_no_debris": _pct(n_no_debris),
            "n_debris": n_debris,
            "pct_debris": _pct(n_debris),
            "n_uncertain": n_uncertain,
            "pct_uncertain": _pct(n_uncertain),
            "mean_debris_pct": round(df[pct_col].mean(), 2) if n_total else None,
            "mean_threshold": round(image_features[threshold_col].mean(), 2) if threshold_col in image_features.columns else None,
            "mean_background": round(image_features[background_col].mean(), 2) if background_col in image_features.columns else None,
        }
        rows.append(row)

    return rows


def build_debris_param_summary(image_features, debris_max_pct, ratio_min):
    """Side-table params: debris method/cell-mask expansion (constant per run, read
    straight off image_features - not exposed as pipeline params, so there's no
    exact param name to key them by) plus the configured pass-criteria thresholds,
    keyed by their exact param name (sc_debris_max_pct/sc_debris_min_ratio) - same
    "just pass a dict through" pattern as tab_flatfield.build_flatfield_param_summary."""
    method = image_features["method"].iloc[0] if "method" in image_features.columns and len(image_features) else "n/a"
    expansion = image_features["cellmask_expansion"].iloc[0] if "cellmask_expansion" in image_features.columns and len(image_features) else "n/a"

    return {
        "Debris method": method,
        "Cell mask expansion": expansion,
        "sc_debris_max_pct": debris_max_pct,
        "sc_debris_min_ratio": ratio_min,
    }


def sample_debris_images(debris_samples_dir):
    """Group the highest-debris overlay PNGs written by make_debris_overlay.py by (sample_class, channel).

    Returns dict sample_class(str) -> dict channel(int) -> list of {"caption", "data_uri"},
    each channel's list sorted worst-of-class first: by debris percentage descending,
    except "no_debris" which sorts by ratio descending (see SAMPLE_CLASS_SORT_KEY).
    """
    if debris_samples_dir is None or not os.path.isdir(debris_samples_dir):
        return {}

    by_class = {}
    for path in sorted(glob.glob(os.path.join(debris_samples_dir, "*_debris.png"))):
        match = DEBRIS_SAMPLE_RE.match(os.path.basename(path))
        if not match:
            log.warning(f"Skipping debris sample image with unexpected name: {path}")
            continue

        data_uri = image_to_data_uri(path)
        if data_uri is None:
            continue

        sample_class = match.group("sample_class")
        channel = int(match.group("channel"))
        pct = float(match.group("pct"))
        ratio = float(match.group("ratio"))
        by_class.setdefault(sample_class, {}).setdefault(channel, []).append({
            "caption": f"{match.group('plate')} / {match.group('well')} / field {match.group('field')} (debris%={pct:.2f}, ratio={ratio:.2f})",
            "data_uri": data_uri,
            "_pct": pct,
            "_ratio": ratio,
        })

    for sample_class, channels in by_class.items():
        sort_key = SAMPLE_CLASS_SORT_KEY.get(sample_class, DEFAULT_SAMPLE_SORT_KEY)
        for images in channels.values():
            images.sort(key=lambda img: img[sort_key], reverse=True)
            for img in images:
                del img["_pct"]
                del img["_ratio"]

    return by_class
