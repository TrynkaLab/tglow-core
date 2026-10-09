"""Tests for the Tab 6 sigmoid plots, with and without the sigmoid_inputs.tsv densities."""

import numpy as np
import pandas as pd

from tglow.qc.tab_scaling import build_sigmoid_plots, load_sigmoid_inputs
from tglow.utils.tglow_utils import sigmoid_params


def _scaling_index(rows):
    df = pd.DataFrame(rows)
    df.index = df["ref_plate"] + "-" + df["channel"].astype(str)
    return df


def _row(plate, channel, x1, x2, tol=1e-3):
    """A scaling_index row as calculate_scaling_factors writes it, slope/bias fitted from x1/x2."""
    par = sigmoid_params(x1, x2, tol=tol)
    return {"ref_plate": plate, "channel": channel, "sigmoid_x1": x1, "sigmoid_x2": x2,
            "sigmoid_slope": par["slope"], "sigmoid_bias": par["bias"], "sigmoid_tol": tol,
            "sigmoid_lower_quantile": 0.95, "sigmoid_upper_quantile": 0.5,
            "sigmoid_lower_feature": f"ch{channel}__background_q75", "sigmoid_upper_feature": f"ch{channel}__otsu_log"}


def _inputs(rng, plate, channel, n, lower_mean=200, upper_mean=1500):
    return pd.DataFrame({"plate": plate, "channel": channel, "well": "A01", "field": 1,
                         "lower": rng.normal(lower_mean, 30, n), "upper": rng.normal(upper_mean, 300, n)})


def _by_kind(fig):
    density = [t for t in fig.data if t.fill == "tozeroy"]
    markers = [t for t in fig.data if t.mode == "markers"]
    curves = [t for t in fig.data if t.mode == "lines" and t.fill != "tozeroy"]
    return curves, density, markers


def test_without_inputs_draws_one_curve_per_plate():
    si = _scaling_index([_row("P1", 0, 250, 1500), _row("P2", 0, 300, 1400)])
    fig = build_sigmoid_plots(si)[0]
    curves, density, markers = _by_kind(fig)
    assert len(curves) == 2 and not density and not markers
    assert fig.layout.yaxis.title.text == "Sigmoid weight"


def test_densities_are_peak_scaled_and_markers_sit_on_x1_x2():
    rng = np.random.default_rng(0)
    si = _scaling_index([_row("P1", 0, 250, 1500), _row("P2", 0, 300, 1400)])
    inputs = pd.concat([_inputs(rng, "P1", 0, 40), _inputs(rng, "P2", 0, 40)])
    fig = build_sigmoid_plots(si, inputs)[0]
    curves, density, markers = _by_kind(fig)

    assert len(curves) == 2 and len(density) == 4 and len(markers) == 4
    for trace in density:
        assert np.isclose(max(trace.y), 1.0)
    assert sorted(m.x[0] for m in markers) == sorted([250, 1500, 300, 1400])
    # Densities drawn first so the curves and markers stay on top
    assert fig.data[0].fill == "tozeroy" and fig.data[-1].mode == "markers"
    # x-range extended to cover the upper density, beyond where the curves saturate
    assert max(fig.data[0].x) >= np.quantile(inputs[["lower", "upper"]].to_numpy(), 0.99) - 1e-9
    assert "q0.95 of ch0__background_q75" in markers[0].name


def test_plate_shares_one_colour_and_legend_group():
    rng = np.random.default_rng(1)
    si = _scaling_index([_row("P1", 0, 250, 1500)])
    fig = build_sigmoid_plots(si, _inputs(rng, "P1", 0, 30))[0]
    assert {t.legendgroup for t in fig.data} == {"P1"}
    assert sum(1 for t in fig.data if t.showlegend is not False) == 1


def test_single_image_plate_gets_marker_only():
    rng = np.random.default_rng(2)
    si = _scaling_index([_row("P1", 0, 250, 1500)])
    fig = build_sigmoid_plots(si, _inputs(rng, "P1", 0, 1))[0]
    curves, density, markers = _by_kind(fig)
    assert len(curves) == 1 and not density and len(markers) == 2
    assert all(m.y[0] == 1.0 for m in markers)


def test_inverted_fit_still_shows_densities_and_dashed_curve():
    rng = np.random.default_rng(3)
    # Borrowed channel-mean curve, but this plate's own x2 <= x1
    row = _row("P1", 0, 250, 1500)
    row.update({"sigmoid_x1": 1600, "sigmoid_x2": 1400, "sigmoid_tol": np.nan})
    si = _scaling_index([row, _row("P2", 0, 250, 1500)])
    inputs = pd.concat([_inputs(rng, "P1", 0, 30, lower_mean=1600, upper_mean=1400), _inputs(rng, "P2", 0, 30)])
    fig = build_sigmoid_plots(si, inputs)[0]
    curves, density, markers = _by_kind(fig)
    p1_curve = [c for c in curves if c.legendgroup == "P1"][0]
    assert p1_curve.line.dash == "dash"
    assert len([d for d in density if d.legendgroup == "P1"]) == 2


def test_channel_without_inputs_is_unchanged():
    rng = np.random.default_rng(4)
    si = _scaling_index([_row("P1", 0, 250, 1500), _row("P1", 1, 400, 2000)])
    figs = build_sigmoid_plots(si, _inputs(rng, "P1", 0, 30))
    curves, density, markers = _by_kind(figs[1])
    assert len(curves) == 1 and not density and not markers


def test_load_sigmoid_inputs_handles_missing_and_empty(tmp_path):
    assert load_sigmoid_inputs(None).empty
    assert load_sigmoid_inputs(str(tmp_path / "NO_SIGMOID_INPUTS")).empty
    empty = tmp_path / "sigmoid_inputs.tsv"
    empty.write_text("")
    assert load_sigmoid_inputs(str(empty)).empty
    header_only = tmp_path / "header.tsv"
    header_only.write_text("plate\tchannel\twell\tfield\tlower\tupper\n")
    assert load_sigmoid_inputs(str(header_only)).empty
