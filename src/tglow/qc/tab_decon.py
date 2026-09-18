"""Tab 4 (deconvolution) - before/after sample images.

The before/after PNGs themselves are produced by `make_decon_before_after.py`
(run once per sampled field by the `qc_decon_samples` Nextflow process, see
processes/qc.nf) - this module locates/loads the before/after PNGs that script
already wrote, named ``<plate>_<well>_<field>_ch<channel>_cells<n_cells>_before_after.png``.

The sampled fields are the ones with the most segmented cells (see
select_decon_samples.py), so each channel's examples are shown densest-first -
decon's effect is easiest to judge where there is actually signal to resolve.
"""

import glob
import logging
import os
import re

from tglow.qc.assets import image_to_data_uri

log = logging.getLogger(__name__)

BEFORE_AFTER_RE = re.compile(
    r"^(?P<plate>.+)_(?P<well>[A-Za-z]+\d+)_(?P<field>\d+)_ch(?P<channel>\d+)_cells(?P<n_cells>\d+)_before_after\.png$"
)


def build_before_after_images(decon_samples_dir):
    """Group the before/after PNGs written by make_decon_before_after.py by channel.

    Each channel's samples are sorted by cell count descending (the densest field
    first), rather than by filename - filename order would sort by plate name.

    Returns dict channel(int) -> list of {"caption": ..., "data_uri": ...}.
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
        n_cells = int(match.group("n_cells"))
        by_channel.setdefault(channel, []).append({
            "caption": f"{match.group('plate')} / {match.group('well')} / field {match.group('field')} ({n_cells} cells)",
            "n_cells": n_cells,
            "data_uri": data_uri,
        })

    for samples in by_channel.values():
        samples.sort(key=lambda sample: -sample["n_cells"])

    return by_channel
