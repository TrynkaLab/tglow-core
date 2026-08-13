"""Plate-format-aware layout helpers for QC report heatmaps.

No plate-format (24/96/384/1536 well) concept exists elsewhere in the pipeline -
wells are always addressed directly by row letter/column number. For a per-plate
heatmap we still need to know the plate's full grid shape (so e.g. a 96-well plate
renders as an 8x12 grid with empty wells left blank, not just a bounding box around
the wells actually present).
"""

import string

import numpy as np

# rows x cols for each supported plate format
PLATE_FORMATS = {
    24: (4, 6),
    96: (8, 12),
    384: (16, 24),
    1536: (32, 48),
}


def infer_plate_format(max_row, max_col, override="auto"):
    """Return the (n_rows, n_cols) grid for a plate given the max row/col index seen.

    `override` is either "auto" (infer from max_row/max_col, rounding up to the
    smallest standard format that fits) or one of the PLATE_FORMATS keys (as an
    int or string) to force a specific layout.
    """
    if override not in (None, "auto"):
        fmt = int(override)
        if fmt not in PLATE_FORMATS:
            raise ValueError(f"qc_plate_format must be one of {sorted(PLATE_FORMATS)} or 'auto', got: {override}")
        return PLATE_FORMATS[fmt]

    for fmt in sorted(PLATE_FORMATS):
        n_rows, n_cols = PLATE_FORMATS[fmt]
        if max_row <= n_rows and max_col <= n_cols:
            return (n_rows, n_cols)

    # Larger than any known format - fall back to a tight bounding box
    return (max_row, max_col)


def row_letter_to_index(row_letter):
    """'A' -> 0, 'B' -> 1, ..., 'AA' -> 26, ... (0-indexed)."""
    idx = 0
    for ch in row_letter.strip().upper():
        idx = idx * 26 + (string.ascii_uppercase.index(ch) + 1)
    return idx - 1


def _row_col_indices(rows, cols):
    """0-indexed (row, col) Series pair - row labels may be letters ('A', 'B', ...)
    or already-numeric row indices (as produced by ImageQuery, 1-indexed)."""
    rows = rows.astype(str)
    cols = cols.astype(int)

    if rows.str.isnumeric().all():
        row_idx = rows.astype(int) - 1
    else:
        row_idx = rows.map(row_letter_to_index)

    return row_idx, cols - 1


def infer_plate_formats(df, plate_col="plate", row_col="row", col_col="col", override="auto"):
    """Infer each plate's (n_rows, n_cols) format ONCE, from one authoritative/complete
    well population (e.g. image_features, which lists every imaged field regardless
    of any downstream cell-level filtering).

    Every heatmap for a given plate should reuse this same dict rather than each
    independently inferring from whatever (possibly filtered/partial) subset of rows
    it happens to see - otherwise the same physical plate can infer to different
    formats in different tabs (e.g. a qc-filtered subset that happens to lack any
    surviving cells in the plate's higher rows/columns would under-infer the format).

    Returns dict plate -> (n_rows, n_cols).
    """
    formats = {}
    for plate, plate_df in df.groupby(plate_col):
        row_idx, col_idx = _row_col_indices(plate_df[row_col], plate_df[col_col])
        max_row = int(row_idx.max()) + 1
        max_col = int(col_idx.max()) + 1
        formats[plate] = infer_plate_format(max_row, max_col, override=override)
    return formats


# Fixed figure size for every plate heatmap, regardless of that plate's own
# (n_rows, n_cols) - without an explicit size, Plotly auto-sizes each figure to
# its own content, so different plate formats (or even the same format loaded in
# a different browser layout pass) render at visibly different widget
# dimensions. All 4 standard plate formats share the same 2:3 row:col ratio, so
# one fixed width/height keeps cells square-ish across every format too.
HEATMAP_WIDTH = 520
HEATMAP_HEIGHT = 350


def style_heatmap_axes(fig):
    """Box border around the plot area, no internal gridlines, and a size fixed
    across every plate/format (see HEATMAP_WIDTH/HEIGHT above)."""
    axis_style = dict(
        showgrid=False, zeroline=False, showline=True, linewidth=1,
        linecolor="rgba(136, 136, 136, 0.5)", mirror=True, ticks="",
    )
    fig.update_xaxes(**axis_style)
    fig.update_yaxes(**axis_style)
    fig.update_layout(autosize=False, width=HEATMAP_WIDTH, height=HEATMAP_HEIGHT)
    return fig


def build_well_grid(df, row_col="row", col_col="col", value_col=None, agg="mean", plate_format="auto"):
    """Build a (n_rows, n_cols) grid of `value_col` (aggregated by `agg`) for one plate's wells.

    `df` is expected to hold one row per well (or per-cell/per-image, in which case
    it's aggregated to one value per well first). `plate_format` is either an
    "auto"/None/int-string override (inferred from `df`'s own max row/col - see
    infer_plate_format) or a precomputed (n_rows, n_cols) tuple (see
    infer_plate_formats - preferred, so every heatmap for one plate agrees).
    Returns (grid, row_labels, col_labels) - grid cells with no data are NaN.
    """
    row_idx, col_idx = _row_col_indices(df[row_col], df[col_col])

    if isinstance(plate_format, tuple):
        n_rows, n_cols = plate_format
    else:
        max_row = int(row_idx.max()) + 1
        max_col = int(col_idx.max()) + 1
        n_rows, n_cols = infer_plate_format(max_row, max_col, override=plate_format)

    grid = np.full((n_rows, n_cols), np.nan)

    if value_col is not None:
        well_values = df.assign(_row=row_idx, _col=col_idx).groupby(["_row", "_col"])[value_col].agg(agg)
        for (r, c), val in well_values.items():
            if r < n_rows and c < n_cols:
                grid[r, c] = val
    else:
        # No value column - just mark well presence (1) vs absence (NaN)
        for r, c in zip(row_idx, col_idx):
            if r < n_rows and c < n_cols:
                grid[r, c] = 1

    row_labels = [string.ascii_uppercase[i % 26] if i < 26 else f"{string.ascii_uppercase[i // 26 - 1]}{string.ascii_uppercase[i % 26]}" for i in range(n_rows)]
    col_labels = [str(i + 1) for i in range(n_cols)]

    return grid, row_labels, col_labels
