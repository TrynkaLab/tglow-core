"""Collect script warnings into a TSV so they survive past .command.err.

Warnings emitted with log.warning alone end up in the task's work directory and are
effectively invisible. WarningLog logs exactly as before and additionally records each
warning, so the script can write a scaling_warnings.tsv next to its other outputs - which
the QC report then reads back and renders (see tglow.qc.render.build_scaling_tab).

Named warning_log rather than warnings to avoid any confusion with the stdlib module.
"""

import logging

import pandas as pd

log = logging.getLogger(__name__)

COLUMNS = ["source", "category", "channel", "message"]


class WarningLog:
    """Accumulates warnings for one script run, then writes them as a TSV.

    `source` identifies the emitting script and is stamped onto every row, so several
    scripts' files can be concatenated without losing track of where each came from.
    """

    def __init__(self, source):
        self.source = source
        self.entries = []

    def warn(self, message, category="general", channel=None):
        """Record a warning and log it, exactly as a bare log.warning would have."""
        log.warning(message)
        self.entries.append({
            "source": self.source,
            "category": category,
            "channel": "" if channel is None else channel,
            "message": message,
        })

    def write(self, path):
        """Write the collected warnings as a TSV.

        Always writes, even with nothing collected - a present-but-empty file is how the
        QC report tells "this script ran and had nothing to report" apart from "this
        script never ran", and it keeps the Nextflow output declaration unconditional.
        """
        df = pd.DataFrame(self.entries, columns=COLUMNS)
        df.to_csv(path, sep="\t", index=False)
        log.info(f"Wrote {len(self.entries)} warning(s) to {path}")
        return df


def load_warnings(paths):
    """Concatenate zero or more scaling_warnings.tsv files into one DataFrame.

    Missing/unreadable paths are skipped rather than raised on - the report should still
    render if a warnings file went astray.
    """
    frames = []

    for path in paths or []:
        try:
            df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        except Exception as exc:
            log.warning(f"Could not read warnings file {path}: {exc}")
            continue
        if not df.empty:
            frames.append(df)

    if not frames:
        return pd.DataFrame(columns=COLUMNS)

    return pd.concat(frames, ignore_index=True)
