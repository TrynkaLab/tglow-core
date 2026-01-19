"""Tglow core package.

This package provides I/O helpers and utilities for indexing and reading
high-content imaging datasets (plate / well / field layouts) and parsers for
PerkinElmer export formats.

Public subpackages:
- ``tglow.io``: image readers, writers and parsers
- ``tglow.utils``: helper utilities for numeric conversions and registration

Keep the top-level module lightweight; import specific components from
submodules (e.g. ``from tglow.io import AICSImageReader``).
"""

__all__ = [
	"io",
	"utils",
]
