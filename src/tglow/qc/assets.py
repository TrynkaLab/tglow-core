"""Base64-inlining helpers so the QC report can be a single, offline-viewable HTML file."""

import base64
import logging
import mimetypes
import os

import cv2
import numpy as np

log = logging.getLogger(__name__)

# Every embedded image gets downscaled to fit within this on either side before
# base64-inlining - otherwise a run with many high-res flatfield/registration/decon
# PNGs makes for a multi-hundred-MB report. Large enough to stay readable at the
# report's own display size (images render at max-width: 100% of a ~320-580px tile).
MAX_IMAGE_DIMENSION = 900


def image_to_data_uri(path, max_dimension=MAX_IMAGE_DIMENSION):
    """Read an image file, downscale it if needed, and return it as a data: URI (or None if missing)."""
    if path is None or not os.path.isfile(path):
        return None

    image = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if image is None:
        # Not decodable as an image (or an unsupported format) - fall back to inlining as-is.
        log.warning(f"Could not decode {path} as an image - embedding original bytes without downscaling")
        mime, _ = mimetypes.guess_type(path)
        with open(path, "rb") as fh:
            encoded = base64.b64encode(fh.read()).decode("ascii")
        return f"data:{mime or 'application/octet-stream'};base64,{encoded}"

    height, width = image.shape[:2]
    longest_side = max(height, width)
    if longest_side > max_dimension:
        scale = max_dimension / longest_side
        image = cv2.resize(image, (max(1, round(width * scale)), max(1, round(height * scale))), interpolation=cv2.INTER_AREA)

    ok, buffer = cv2.imencode(".png", image)
    if not ok:
        raise RuntimeError(f"Failed to re-encode {path} as PNG")

    encoded = base64.b64encode(np.asarray(buffer)).decode("ascii")
    return f"data:image/png;base64,{encoded}"
