"""Characterization tests for AICSImageReader / AICSImageWriter in tglow_io.py.

Written against the aicsimageio-backed implementation to lock in current
behavior, then re-run unchanged against the bioio-backed implementation to
prove the migration is behavior-preserving. Covers exactly the attributes
(`.dims`, `.channel_names`, `.dtype`) that `compound_image_provider.py` and
`processed_image_provider.py` rely on downstream.
"""

import os

import numpy as np

from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader, AICSImageWriter

# Kept in sync with the fixture in conftest.py by value, not by import, since
# test modules aren't part of a package under pytest's default rootless mode.
PLATE = "PLATE1"
CHANNEL_NAMES = ["C0", "C1"]
STACK_SHAPE = (2, 3, 8, 8)  # CZYX


def test_get_img_exposes_dims_channel_names_dtype(plate_dir):
    reader = AICSImageReader(plate_dir["path"])
    img = reader.get_img(plate_dir["query"])

    assert img.dims["C"][0] == STACK_SHAPE[0]
    assert "C" in img.dims.order
    assert list(img.channel_names) == CHANNEL_NAMES
    assert img.dtype == np.uint16


def test_read_image_full_stack(plate_dir):
    reader = AICSImageReader(plate_dir["path"])
    query = ImageQuery(PLATE, 1, 1, "001")

    result = reader.read_image(query)

    assert result.shape == STACK_SHAPE
    np.testing.assert_array_equal(result, plate_dir["stack"])


def test_read_image_channel_only(plate_dir):
    reader = AICSImageReader(plate_dir["path"])
    query = ImageQuery(PLATE, 1, 1, "001", channel="0")

    result = reader.read_image(query)

    assert result.shape == STACK_SHAPE[1:]
    np.testing.assert_array_equal(result, plate_dir["stack"][0])


def test_read_image_plane_only(plate_dir):
    reader = AICSImageReader(plate_dir["path"])
    query = ImageQuery(PLATE, 1, 1, "001", plane="1")

    result = reader.read_image(query)

    assert result.shape == (STACK_SHAPE[0], STACK_SHAPE[2], STACK_SHAPE[3])
    np.testing.assert_array_equal(result, plate_dir["stack"][:, 1, :, :])


def test_read_image_channel_and_plane(plate_dir):
    reader = AICSImageReader(plate_dir["path"])
    query = ImageQuery(PLATE, 1, 1, "001", channel="0", plane="1")

    result = reader.read_image(query)

    assert result.shape == STACK_SHAPE[2:]
    np.testing.assert_array_equal(result, plate_dir["stack"][0, 1, :, :])


def test_read_stack_resets_channel_and_plane(plate_dir):
    reader = AICSImageReader(plate_dir["path"])
    query = ImageQuery(PLATE, 1, 1, "001", channel="0", plane="1")

    result = reader.read_stack(query)

    assert query.channel is None
    assert query.plane is None
    assert result.shape == STACK_SHAPE
    np.testing.assert_array_equal(result, plate_dir["stack"])


def test_write_image_stats(tmp_path):
    writer = AICSImageWriter(str(tmp_path), channel_names=CHANNEL_NAMES)
    rng = np.random.default_rng(1)
    stack = rng.integers(0, 65535, size=STACK_SHAPE, dtype=np.uint16)
    query = ImageQuery(PLATE, 1, 1, "001")

    writer.write_stack(stack, query)
    writer.write_image_stats(query)

    stats_path = os.path.join(str(tmp_path), PLATE, "A", "1", "intensity_stats.tsv")
    assert os.path.exists(stats_path)

    with open(stats_path) as fh:
        lines = fh.readlines()

    # header + one row per channel
    assert len(lines) == 1 + STACK_SHAPE[0]
