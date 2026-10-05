"""Tab 6 (scaling factors) - scale factor barplot + sigmoid curves, from scaling_index.tsv.

Only produced when sc_autoscale is enabled (calculate_scaling_factors.py is the
only thing that writes scaling_index.tsv - manual scaling never does).
"""

import logging

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from tglow.qc.assets import style_plot
from tglow.utils.tglow_utils import sigmoid

log = logging.getLogger(__name__)


def load_scaling_index(scaling_index_path):
    return pd.read_csv(scaling_index_path, sep="\t", index_col=0)


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


def build_sigmoid_plots(scaling_index, n_points=200, default_tol=1e-3):
    """dict[channel] -> Plotly figure with one sigmoid curve per plate.

    The x-range comes from the curve itself - [0, x where the sigmoid reaches 1 - tol] - not
    from sigmoid_x2: for a fitted plate the two coincide, but a plate whose slope/bias were
    borrowed from the channel mean keeps its own (inverted, or missing) x2, which can sit far
    below the borrowed midpoint and plot nothing but ~0. Borrowed curves are drawn dashed.
    """
    figures = {}

    for channel, channel_df in scaling_index.groupby("channel"):
        # fill_missing_sigmoid copies slope/bias but not tol, so a borrowed row falls back to
        # the channel's fitted tol
        channel_tol = channel_df["sigmoid_tol"].dropna() if "sigmoid_tol" in channel_df else pd.Series(dtype=float)
        channel_tol = channel_tol.median() if not channel_tol.empty else default_tol

        fig = go.Figure()
        for _, row in channel_df.iterrows():
            slope, bias = row.get("sigmoid_slope"), row.get("sigmoid_bias")
            if pd.isna(slope) or pd.isna(bias) or slope <= 0:
                continue

            tol = row.get("sigmoid_tol")
            tol = channel_tol if pd.isna(tol) else tol
            x_max = bias + np.log(1 / tol - 1) / slope
            if x_max <= 0:
                continue

            borrowed = _is_borrowed_sigmoid(row)
            name = str(row["ref_plate"]) + (" (channel mean)" if borrowed else "")
            x = np.linspace(0, x_max, n_points)
            y = sigmoid(x, slope, bias)
            fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=name,
                                     line=dict(dash="dash" if borrowed else "solid")))

        fig.update_layout(
            title=f"Channel {channel} sigmoid soft-threshold curve",
            xaxis_title="Intensity",
            yaxis_title="Sigmoid weight",
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
