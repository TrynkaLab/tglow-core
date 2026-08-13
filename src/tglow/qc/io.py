"""Shared parquet-loading helpers for measure_intensity output.

Both ``calculate_scaling_factors.py`` and the QC report read the same
per-plate directory layout produced by ``measure_intensity`` + ``stage_as_plate``
(one subdirectory per plate, each holding that plate's well-level
object_features.parquet/image_features.parquet files).
"""

import glob
import logging
import os

import pandas as pd

log = logging.getLogger(__name__)


def load_measurements(input_dir, name_fragment):
    """Concatenate <input_dir>/<plate>/<well_subdir>/*<name_fragment>*.parquet across every plate.

    input_dir is the directory of per-plate subdirectories produced by stage_as_plate.
    Every well contributes an identically-named object_features.parquet/
    image_features.parquet, so stage_as_plate stages each well's pair into its own
    numbered subdirectory under the plate dir (a flat layout would collide on
    filenames) - it doesn't matter which well a file came from here, since every
    parquet file already carries its own plate/row/col/field/well columns.
    """
    paths = sorted(glob.glob(os.path.join(input_dir, "*", "*", f"*{name_fragment}*.parquet")))

    if not paths:
        raise RuntimeError(f"No {name_fragment} parquet files found under {input_dir}")

    return pd.concat([pd.read_parquet(p) for p in paths], ignore_index=True)
