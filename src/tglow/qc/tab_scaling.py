"""Tab 6 (scaling factors) - scale factor barplot + sigmoid curves, from scaling_index.tsv.

Only produced when sc_autoscale is enabled (calculate_scaling_factors.py is the
only thing that writes scaling_index.tsv - manual scaling never does).
"""

import logging

import numpy as np
import pandas as pd
import plotly.graph_objects as go

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
    return fig


def build_sigmoid_plots(scaling_index, n_points=200):
    """dict[channel] -> Plotly figure with one sigmoid curve per plate, x in [0, sigmoid_x2]."""
    figures = {}

    for channel, channel_df in scaling_index.groupby("channel"):
        fig = go.Figure()
        for _, row in channel_df.iterrows():
            x2 = row.get("sigmoid_x2")
            if pd.isna(x2) or x2 <= 0:
                continue
            x = np.linspace(0, x2, n_points)
            y = sigmoid(x, row["sigmoid_slope"], row["sigmoid_bias"])
            fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=str(row["ref_plate"])))

        fig.update_layout(
            title=f"Channel {channel} sigmoid soft-threshold curve",
            xaxis_title="Intensity",
            yaxis_title="Sigmoid weight",
        )
        figures[channel] = fig

    return figures
