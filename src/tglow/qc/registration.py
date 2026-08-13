"""Registration-correlation filtering, shared by scaling-factor estimation and the QC report.

Both consumers of measure_intensity's object_features parquet need the same "keep
only cells whose per-channel-pair registration correlation clears a threshold"
filter - calculate_scaling_factors.py uses it (sc_registration_thresh/pattern) to
decide which cells to compute scaling factors from, and the QC report uses it
(qc_regcor/qc_registration_pattern) to report qc'ed cell counts/distributions.
"""

import logging

log = logging.getLogger(__name__)


def filter_registration_correlation(cell_df, pattern, threshold):
    """Drop objects whose registration correlation falls below threshold in any matching channel column.

    Columns matching `pattern` (substring match, e.g. "registration_corr") report the
    correlation between channels used for registration. An object is kept only if every
    matching column is >= threshold; a missing/NaN value counts as a failure.
    """
    corr_cols = [c for c in cell_df.columns if pattern in c]

    if not corr_cols:
        log.warning(f"No columns matching registration feature pattern '{pattern}' found - skipping registration-based filtering")
        return cell_df

    passes = cell_df[corr_cols].ge(threshold).all(axis=1)
    n_removed = (~passes).sum()
    log.info(
        f"Registration filtering ({corr_cols}, threshold={threshold}): "
        f"removing {n_removed}/{len(cell_df)} objects"
    )

    return cell_df.loc[passes].reset_index(drop=True)


def registration_correlation_columns(cell_df, pattern):
    """Return the list of columns in cell_df matching pattern (substring match)."""
    return [c for c in cell_df.columns if pattern in c]
