"""Tab 4 (deconvolution) - PSF slices + before/after sample images.

The before/after PNGs themselves are produced by `make_decon_before_after.py`
(run once per sampled well by the `qc_decon_samples` Nextflow process, see
processes/qc.nf) - this module only builds the PSF figures and locates/loads the
before/after PNGs that script already wrote, named
``<plate>_<well>_ch<channel>_before_after.png``.
"""

import glob
import logging
import os
import re

import numpy as np
import plotly.graph_objects as go
import tifffile

from tglow.qc.assets import image_to_data_uri

log = logging.getLogger(__name__)

BEFORE_AFTER_RE = re.compile(r"^(?P<plate>.+)_(?P<well>[A-Za-z]+\d+)_ch(?P<channel>\d+)_before_after\.png$")


def build_psf_figures(psf_paths):
    """psf_paths: dict of 0-indexed channel -> PSF tiff path. Returns dict channel -> {"zx": fig, "xy": fig}."""
    figures = {}

    for channel, path in psf_paths.items():
        if not os.path.isfile(path):
            log.warning(f"PSF file for channel {channel} not found: {path}")
            continue

        psf = tifffile.imread(path)

        if psf.ndim != 3:
            log.warning(f"PSF for channel {channel} is not a 3D (Z,Y,X) stack (shape {psf.shape}) - skipping ZX slice")
            xy = psf if psf.ndim == 2 else np.max(psf, axis=0)
            zx = None
        else:
            y_mid = psf.shape[1] // 2
            zx = psf[:, y_mid, :]
            xy = np.max(psf, axis=0)

        entry = {}
        if xy is not None:
            entry["xy"] = go.Figure(data=go.Heatmap(z=xy, colorscale="gray", showscale=False))
            entry["xy"].update_layout(title=f"Channel {channel} PSF - XY (max projected)", yaxis=dict(autorange="reversed", scaleanchor="x"))
        if zx is not None:
            entry["zx"] = go.Figure(data=go.Heatmap(z=zx, colorscale="gray", showscale=False))
            entry["zx"].update_layout(title=f"Channel {channel} PSF - ZX slice", yaxis=dict(autorange="reversed", scaleanchor="x"))

        figures[channel] = entry

    return figures


def build_before_after_images(decon_samples_dir):
    """Group the before/after PNGs written by make_decon_before_after.py by channel.

    Returns dict channel(int) -> list of {"plate": ..., "well": ..., "data_uri": ...}.
    """
    if decon_samples_dir is None or not os.path.isdir(decon_samples_dir):
        return {}

    by_channel = {}
    for path in sorted(glob.glob(os.path.join(decon_samples_dir, "*_before_after.png"))):
        match = BEFORE_AFTER_RE.match(os.path.basename(path))
        if not match:
            log.warning(f"Skipping decon sample image with unexpected name: {path}")
            continue

        data_uri = image_to_data_uri(path)
        if data_uri is None:
            continue

        channel = int(match.group("channel"))
        by_channel.setdefault(channel, []).append({
            "plate": match.group("plate"),
            "well": match.group("well"),
            "data_uri": data_uri,
        })

    return by_channel
