"""Base64-inlining and shared plot helpers, so the QC report can be a single,
offline-viewable HTML file that doesn't grow without bound."""

import base64
import logging
import mimetypes
import os

import cv2
import numpy as np
import plotly.graph_objects as go

log = logging.getLogger(__name__)

# Every embedded image gets downscaled to fit within this on either side before
# base64-inlining - otherwise a run with many high-res flatfield/registration/decon
# PNGs makes for a multi-hundred-MB report. The carousels display at up to 900 CSS px
# (.carousel's max-width in the template), so 1200 leaves ~1.33x pixel density: crisp
# on a HiDPI screen and with a little room to zoom. Flatfield tiles are smaller again
# (~320px grid items) and are oversampled either way. Trades off directly against
# IMAGE_ENCODE_PARAMS below - more pixels at lower quality generally reads sharper than
# the reverse, which is why that quality dropped to 70 as this went up.
MAX_IMAGE_DIMENSION = 1200

# Embedded images are re-encoded as WebP rather than PNG: they're all 8-bit display
# renderings (matplotlib figures, overlay crops), where WebP measured 4-5x smaller than
# the equivalent PNG - which matters doubly once base64 inflates whatever we produce by
# a further 33%. PNG stays the fallback for anything that isn't 8-bit, where a lossy
# re-encode would need a normalization choice this helper has no business making, and
# for OpenCV builds without WebP support.
#
# Quality 70 rather than 80 is the other half of raising MAX_IMAGE_DIMENSION to 1200:
# q70 at 1200px costs about 1.4x what q80 at 900px did (q80 at 1200px would be 2.1x)
# while carrying a third more pixels, and on a downscaled grayscale frame the extra
# resolution reads sharper than the extra quality. The thing to re-check if this is
# lowered further is the debris overlay's thin red contour - fine linework is what
# lossy compression damages first.
IMAGE_FORMAT = ".webp"
IMAGE_ENCODE_PARAMS = [cv2.IMWRITE_WEBP_QUALITY, 70]
IMAGE_MIME = "image/webp"

_webp_available = None


def _webp_supported():
    """One-off probe - WebP support depends on how the installed OpenCV was built.

    Cached (and warned about once) rather than per image, so a build without it
    degrades to PNG quietly instead of logging for every sample in the report.
    """
    global _webp_available

    if _webp_available is None:
        try:
            ok, _ = cv2.imencode(IMAGE_FORMAT, np.zeros((2, 2, 3), np.uint8), IMAGE_ENCODE_PARAMS)
            _webp_available = bool(ok)
        except cv2.error:
            _webp_available = False

        if not _webp_available:
            log.warning("This OpenCV build cannot encode WebP - embedding images as PNG instead (larger report)")

    return _webp_available


def _encode_image(image):
    """Encode an image array for inlining, returning (buffer, mime).

    Drops a fully-opaque alpha channel first: matplotlib writes RGBA PNGs by
    default, and an all-255 alpha plane is a quarter of the pixel data carrying
    no information. A non-opaque alpha is kept, since flattening it would turn
    transparent regions black.
    """
    if image.ndim == 3 and image.shape[2] == 4 and (image[:, :, 3] == 255).all():
        image = image[:, :, :3]

    if image.dtype == np.uint8 and _webp_supported():
        ok, buffer = cv2.imencode(IMAGE_FORMAT, image, IMAGE_ENCODE_PARAMS)
        if ok:
            return buffer, IMAGE_MIME
        log.warning(f"WebP encoding produced no data for a {image.shape} {image.dtype} image - falling back to PNG")

    ok, buffer = cv2.imencode(".png", image)
    return (buffer, "image/png") if ok else (None, None)


def histogram_bar(values, bins=50, density=False, **bar_kwargs):
    """A pre-binned histogram as a go.Bar trace.

    go.Histogram bins client-side, which means the trace embeds every raw value in
    the HTML. With one trace per channel/feature over every qc'ed cell that becomes
    the single largest thing in the report - hundreds of MB on a large run, versus
    the ~8 KB this produces regardless of cell count. The rendered chart is
    equivalent; what's lost is Plotly re-binning as you zoom.

    bins accepts anything np.histogram does, so callers overlaying several
    distributions can pass a shared set of edges to keep them comparable.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return go.Bar(x=[], y=[], **bar_kwargs)

    counts, edges = np.histogram(values, bins=bins, density=density)
    return go.Bar(x=(edges[:-1] + edges[1:]) / 2, y=counts, width=np.diff(edges), **bar_kwargs)


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

    buffer, mime = _encode_image(image)
    if buffer is None:
        raise RuntimeError(f"Failed to re-encode {path} for inlining")

    encoded = base64.b64encode(np.asarray(buffer)).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def style_plot(fig, square=True, white_bg=True):
    """Shared QC report plot styling - a title centered and tight to the plot area,
    plus (disable for heatmaps, which keep their default light-grey plot area so
    empty wells stay visible) a roughly square aspect ratio and a white background."""
    fig.update_layout(title=dict(x=0.5, xanchor="center", font=dict(size=13)), margin=dict(t=32))
    if square:
        fig.update_layout(autosize=False, width=480, height=480)
    if white_bg:
        fig.update_layout(plot_bgcolor="white", paper_bgcolor="white")
    return fig
