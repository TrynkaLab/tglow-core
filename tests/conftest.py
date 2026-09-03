"""Shared fixtures for the tglow_io characterization tests.

The fixture writes a tiny synthetic plate through `AICSImageWriter` itself
(rather than calling the underlying imaging library directly), so the same
fixture and tests run unmodified whether `tglow_io.py` is backed by
aicsimageio or bioio.
"""

import numpy as np
import pytest

from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader, AICSImageWriter

PLATE = "PLATE1"
CHANNEL_NAMES = ["C0", "C1"]
STACK_SHAPE = (2, 3, 8, 8)  # CZYX


@pytest.fixture
def plate_dir(tmp_path):
    """Write a single CZYX field into a /plate/row/col/field.ome.tiff tree."""
    writer = AICSImageWriter(
        str(tmp_path), channel_names=CHANNEL_NAMES, skip_imagestats=True
    )

    rng = np.random.default_rng(0)
    stack = rng.integers(0, 65535, size=STACK_SHAPE, dtype=np.uint16)

    query = ImageQuery(PLATE, 1, 1, "001")
    writer.write_stack(stack, query)

    return {
        "path": str(tmp_path),
        "stack": stack,
        "query": query,
    }
