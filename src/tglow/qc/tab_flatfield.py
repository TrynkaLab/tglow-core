"""Tab 3 (flatfields) - estimation parameters + per-channel flatfield/evaluation PNGs.

run_flatfield_estimation.py writes its outputs under
``<rn_publish_dir>/flatfields/<plate>/<plate>_ch<channel>/`` - this is true whether
ff_global_flatfield is set or not, since stage_global_flatfield copies the single
shared global model into every plate's own folder under that same naming
convention. So when global, every plate's copy is identical - we only need to
show one representative plate per channel instead of repeating it.

No metadata file is written alongside the PNGs, so the "what parameters made this
flatfield" summary comes directly from the pipeline's own ff_* params, passed in
here rather than read from disk.
"""

import glob
import logging
import os

log = logging.getLogger(__name__)

FLATFIELD_PNG = "flat_and_darkfield.png"
# model_evaluation_nimg_<nimg_validate>_nbin_<bins>.png - nimg_validate/bins vary by
# run, so glob for it rather than hardcoding the exact filename.
EVALUATION_PNG_GLOBS = ["model_evaluation_nimg_*.png", "all_imgs_max_proj_pre_post.png"]


def build_flatfield_param_summary(ff_params):
    """ff_params is a dict of the raw ff_* pipeline params - just pass through for display."""
    return ff_params


def _find_first(directory, patterns):
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.join(directory, pattern)))
        if matches:
            return matches[0]
    return None


def build_flatfield_images(flatfields_dir, plate_channels, ff_global_flatfield):
    """Locate flatfield/evaluation PNGs per (plate, channel).

    `plate_channels` is a dict of plate -> list of 0-indexed channel ints (as used
    in the `<plate>_ch<channel>` folder naming). Returns a dict keyed by channel,
    each holding a list of per-plate entries: {"plate": ..., "flatfield_png": ...,
    "evaluation_png": ...} (png fields are data-URI-ready file paths, or None if
    missing). When ff_global_flatfield is set, only the first plate per channel is
    included (every plate's copy is identical).
    """
    if flatfields_dir is None or not os.path.isdir(flatfields_dir):
        return {}

    by_channel = {}
    for plate in sorted(plate_channels):
        for channel in plate_channels[plate]:
            by_channel.setdefault(channel, [])

            if ff_global_flatfield and by_channel[channel]:
                # Already have a representative plate for this channel - global
                # flatfields are byte-identical copies, so skip the rest.
                continue

            model_dir = os.path.join(flatfields_dir, plate, f"{plate}_ch{channel}")
            flatfield_png = os.path.join(model_dir, FLATFIELD_PNG)
            if not os.path.isfile(flatfield_png):
                flatfield_png = None

            evaluation_png = _find_first(model_dir, EVALUATION_PNG_GLOBS) if os.path.isdir(model_dir) else None

            by_channel[channel].append({
                "plate": plate,
                "label": "Global (shared across all plates)" if ff_global_flatfield else plate,
                "flatfield_png": flatfield_png,
                "evaluation_png": evaluation_png,
            })

    return by_channel
