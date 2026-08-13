"""Base64-inlining helpers so the QC report can be a single, offline-viewable HTML file."""

import base64
import mimetypes
import os


def image_to_data_uri(path):
    """Read an image file and return it as a data: URI, or None if it doesn't exist."""
    if path is None or not os.path.isfile(path):
        return None

    mime, _ = mimetypes.guess_type(path)
    mime = mime or "image/png"

    with open(path, "rb") as fh:
        encoded = base64.b64encode(fh.read()).decode("ascii")

    return f"data:{mime};base64,{encoded}"
