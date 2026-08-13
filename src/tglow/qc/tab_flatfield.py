"""Tab 3 (flatfields) - estimation parameters + per-channel flatfield/evaluation PNGs.

run_flatfield_estimation.py writes its outputs under
``<rn_publish_dir>/flatfields/<plate>/<plate>_ch<channel>/`` - this is true whether
ff_global_flatfield is set or not, since stage_global_flatfield copies the shared
global model into every plate's own folder under that same naming convention.

When ff_global_flatfield is set, "global" is *per registration cycle*, not
necessarily across the whole run: with a registration manifest,
flatfield_estimation.nf groups plates by cycle index (all reference plates share
one model, all first-query plates share another, etc.), so there can be more than
one distinct global model in play (e.g. one per cycle). We can't tell from the
pipeline params alone how many distinct models exist or which plates share one -
so entries are deduplicated by content hash of their flatfield PNG instead of by
plate: plates whose flatfield is byte-identical are shown once (labeled with every
plate that shares it); plates with a genuinely different model each get their own
entry. Without registration (single cycle == whole run), this still collapses to
one entry when ff_global_flatfield is set, same as before.

No metadata file is written alongside the PNGs, so the "what parameters made this
flatfield" summary comes directly from the pipeline's own ff_* params, passed in
here rather than read from disk.
"""

import glob
import hashlib
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


def _file_hash(path):
    if path is None or not os.path.isfile(path):
        return None
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def build_flatfield_images(flatfields_dir, plate_channels, ff_global_flatfield):
    """Locate flatfield/evaluation PNGs per (plate, channel).

    `plate_channels` is a dict of plate -> list of 0-indexed channel ints (as used
    in the `<plate>_ch<channel>` folder naming). Returns a dict keyed by channel,
    each holding a list of entries: {"plates": [...], "label": ..., "flatfield_png":
    ..., "evaluation_png": ...} (png fields are data-URI-ready file paths, or None
    if missing) - one entry per distinct model (see module docstring for why
    dedup is by content hash rather than by the ff_global_flatfield flag alone).
    """
    if flatfields_dir is None or not os.path.isdir(flatfields_dir):
        return {}

    by_channel = {}
    for plate in sorted(plate_channels):
        for channel in plate_channels[plate]:
            model_dir = os.path.join(flatfields_dir, plate, f"{plate}_ch{channel}")
            flatfield_png = os.path.join(model_dir, FLATFIELD_PNG)
            if not os.path.isfile(flatfield_png):
                flatfield_png = None

            evaluation_png = _find_first(model_dir, EVALUATION_PNG_GLOBS) if os.path.isdir(model_dir) else None

            entries = by_channel.setdefault(channel, [])
            content_hash = _file_hash(flatfield_png) if ff_global_flatfield else None

            existing = next((e for e in entries if content_hash is not None and e["_hash"] == content_hash), None)
            if existing is not None:
                existing["plates"].append(plate)
            else:
                entries.append({
                    "plates": [plate],
                    "flatfield_png": flatfield_png,
                    "evaluation_png": evaluation_png,
                    "_hash": content_hash,
                })

    for entries in by_channel.values():
        for entry in entries:
            entry["label"] = ", ".join(entry["plates"])
            del entry["_hash"]

    return by_channel
