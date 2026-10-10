"""Tab 6 (scaling factors) - scale factor barplot + sigmoid curves, from scaling_index.tsv.

Only produced when sc_autoscale is enabled (calculate_scaling_factors.py is the
only thing that writes scaling_index.tsv - manual scaling never does).
"""

import logging

import os

import numpy as np
import pandas as pd
import plotly.colors
import plotly.graph_objects as go
from scipy.stats import gaussian_kde

from tglow.qc.assets import style_plot
from tglow.utils.tglow_utils import sigmoid

log = logging.getLogger(__name__)

SIGMOID_INPUT_COLS = ["plate", "channel", "well", "field", "lower", "upper"]


def load_scaling_index(scaling_index_path):
    return pd.read_csv(scaling_index_path, sep="\t", index_col=0)


def load_sigmoid_inputs(sigmoid_inputs_path):
    """sigmoid_inputs.tsv from calculate_scaling_factors, or an empty frame when absent/empty (e.g. a stub run's touched file)."""
    if sigmoid_inputs_path is None or not os.path.exists(sigmoid_inputs_path) or os.path.getsize(sigmoid_inputs_path) == 0:
        return pd.DataFrame(columns=SIGMOID_INPUT_COLS)

    df = pd.read_csv(sigmoid_inputs_path, sep="\t", dtype={"plate": str})
    return df


def build_scale_factor_barplot(scaling_index):
    """Grouped barplot: x=channel, y=scale_factor, grouped/colored by plate."""
    fig = go.Figure()
    for plate, plate_df in scaling_index.groupby("ref_plate"):
        plate_df = plate_df.sort_values("channel")
        fig.add_trace(go.Bar(
            x=["ch" + str(c) for c in plate_df["channel"]],
            y=plate_df["scale_factor"],
            name=str(plate),
        ))
    fig.update_layout(
        title="Scaling factor per channel",
        xaxis_title="Channel",
        yaxis_title="Scale factor",
        barmode="group",
    )
    style_plot(fig, square=False)
    return fig


def _is_borrowed_sigmoid(row):
    """True when a plate's slope/bias were filled from the channel mean rather than fitted.

    calculate_scaling_factors.py's fill_missing_sigmoid fills two cases: an inverted fit
    (x2 <= x1 - x1/x2 are still written, but slope/bias are not) and a plate with no fit at
    all (no control samples, so x1/x2 are NaN). Either way x1/x2 don't describe the curve.
    """
    x1, x2 = row.get("sigmoid_x1"), row.get("sigmoid_x2")
    return pd.isna(x1) or pd.isna(x2) or x2 <= x1


def _peak_scaled_density(values, grid, marker_x=()):
    """Gaussian KDE of values, divided by its maximum so the peak is 1, or None when it can't be fit.

    Returns (kde, peak, x, y). x is grid plus the in-range data values and marker_x: a
    distribution much narrower than the grid spacing (a tight background, often with many
    tied integer values) peaks between grid points, so on the grid alone the peak is
    underestimated - markers then land above 1 - and a marker at an arbitrary x falls
    between vertices, off the drawn line. Evaluating at the data and markers resolves both.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 2 or np.ptp(values) == 0:
        return None

    try:
        kde = gaussian_kde(values)
    except np.linalg.LinAlgError:
        return None

    in_range = values[(values >= grid[0]) & (values <= grid[-1])]
    marker_x = np.asarray(marker_x, dtype=float)
    x = np.union1d(np.union1d(grid, in_range), marker_x[np.isfinite(marker_x)])
    y = kde(x)
    peak = y.max()
    if not np.isfinite(peak) or peak <= 0:
        return None

    return kde, peak, x, y / peak


def _rgba(hex_color, alpha):
    r, g, b = plotly.colors.hex_to_rgb(hex_color)
    return f"rgba({r},{g},{b},{alpha})"


def _marker_label(row, plate, which):
    """Hover text for an x1/x2 marker, naming the quantile and feature when scaling_index carries them."""
    value = row.get(f"sigmoid_{which}")
    side = "lower" if which == "x1" else "upper"
    quantile, feature = row.get(f"sigmoid_{side}_quantile"), row.get(f"sigmoid_{side}_feature")
    label = f"{plate} {which} = {value:.1f}"
    if not pd.isna(quantile) and not pd.isna(feature):
        label += f" (q{quantile:g} of {feature})"
    return label


def build_sigmoid_plots(scaling_index, sigmoid_inputs=None, n_points=200, default_tol=1e-3):
    """dict[channel] -> Plotly figure with one sigmoid curve per plate.

    The x-range comes from the curve itself - [0, x where the sigmoid reaches 1 - tol] - not
    from sigmoid_x2: for a fitted plate the two coincide, but a plate whose slope/bias were
    borrowed from the channel mean keeps its own (inverted, or missing) x2, which can sit far
    below the borrowed midpoint and plot nothing but ~0. Borrowed curves are drawn dashed.

    sigmoid_inputs (calculate_scaling_factors' sigmoid_inputs.tsv) adds, per plate, the
    control-image distributions x1/x2 were taken from: the sigmoid_lower_feature values
    and sigmoid_upper_feature values, as shaded KDEs (no outline) peak-scaled to 1 so they
    share the sigmoid's 0-1 axis, with markers where x1/x2 land. The x-range is then
    extended to cover them. Channels without inputs are drawn exactly as before.
    """
    figures = {}
    if sigmoid_inputs is None:
        sigmoid_inputs = pd.DataFrame(columns=SIGMOID_INPUT_COLS)

    # One colour per plate, shared by its curve, densities and markers, and stable across channels
    palette = plotly.colors.qualitative.Plotly
    plates = sorted(scaling_index["ref_plate"].astype(str).unique())
    plate_color = {plate: palette[i % len(palette)] for i, plate in enumerate(plates)}

    for channel, channel_df in scaling_index.groupby("channel"):
        # fill_missing_sigmoid copies slope/bias but not tol, so a borrowed row falls back to
        # the channel's fitted tol
        channel_tol = channel_df["sigmoid_tol"].dropna() if "sigmoid_tol" in channel_df else pd.Series(dtype=float)
        channel_tol = channel_tol.median() if not channel_tol.empty else default_tol

        channel_inputs = sigmoid_inputs[sigmoid_inputs["channel"] == channel]
        has_inputs = not channel_inputs.empty

        # Curve parameters first, so the shared x-range is known before anything is drawn
        curves = []
        for _, row in channel_df.iterrows():
            slope, bias = row.get("sigmoid_slope"), row.get("sigmoid_bias")
            if pd.isna(slope) or pd.isna(bias) or slope <= 0:
                continue

            tol = row.get("sigmoid_tol")
            tol = channel_tol if pd.isna(tol) else tol
            x_max = bias + np.log(1 / tol - 1) / slope
            if x_max <= 0:
                continue
            curves.append((row, slope, bias, x_max))

        x_end = max([c[3] for c in curves], default=0)
        if has_inputs:
            # Cover the densities too, otherwise everything above a median x2 is cut off
            x_end = max(x_end, np.nanquantile(channel_inputs[["lower", "upper"]].to_numpy(dtype=float), 0.99))
        grid = np.linspace(0, x_end, n_points) if x_end > 0 else None

        density_traces, curve_traces, marker_traces = [], [], []
        plates_with_curve = {str(c[0]["ref_plate"]) for c in curves}

        for row, slope, bias, x_max in curves:
            plate = str(row["ref_plate"])
            borrowed = _is_borrowed_sigmoid(row)
            name = plate + (" (channel mean)" if borrowed else "")
            # Without inputs each curve keeps its own range, exactly as before
            x = grid if has_inputs else np.linspace(0, x_max, n_points)
            y = sigmoid(x, slope, bias)
            curve_traces.append(go.Scatter(x=x, y=y, mode="lines", name=name, legendgroup=plate,
                                           line=dict(dash="dash" if borrowed else "solid", color=plate_color[plate])))

        if has_inputs and grid is not None:
            for _, row in channel_df.iterrows():
                plate = str(row["ref_plate"])
                plate_inputs = channel_inputs[channel_inputs["plate"].astype(str) == plate]
                if plate_inputs.empty:
                    continue

                color = plate_color[plate]
                # A plate with inputs but no drawable curve still needs a legend entry
                show_legend = plate not in plates_with_curve

                for which, col in [("x1", "lower"), ("x2", "upper")]:
                    value = row.get(f"sigmoid_{which}")
                    density = _peak_scaled_density(plate_inputs[col], grid, [] if pd.isna(value) else [value])
                    if density is not None:
                        kde, peak, x, y = density
                        density_traces.append(go.Scatter(
                            x=x, y=y, mode="lines", fill="tozeroy", fillcolor=_rgba(color, 0.15),
                            line=dict(width=0), name=plate, legendgroup=plate,
                            showlegend=show_legend, hoverinfo="skip"))
                        show_legend = False

                    if pd.isna(value):
                        continue
                    # value is one of the density's x points, so the marker sits on a vertex of the line
                    marker_y = kde(value)[0] / peak if density is not None else 1.0
                    marker_traces.append(go.Scatter(
                        x=[value], y=[marker_y], mode="markers", legendgroup=plate, showlegend=False,
                        marker=dict(color=color, size=10),
                        name=_marker_label(row, plate, which), hovertemplate="%{fullData.name}<extra></extra>"))

        # Densities underneath, then the curves, then the markers on top
        fig = go.Figure(data=density_traces + curve_traces + marker_traces)

        fig.update_layout(
            title=f"Channel {channel} sigmoid soft-threshold curve",
            xaxis_title="Intensity",
            yaxis_title="Sigmoid weight / scaled density" if has_inputs else "Sigmoid weight",
        )
        # square=False drops style_plot's fixed 480x480 so the curve spans the full column
        # width: it is read along the intensity axis, and several plates' curves sit close
        # together in x, so horizontal room is what makes them separable. Width is left to the
        # container rather than pinned, which would overflow the sidebar layout on a narrower
        # window; only the height is fixed.
        style_plot(fig, square=False)
        fig.update_layout(height=460)
        figures[channel] = fig

    return figures
