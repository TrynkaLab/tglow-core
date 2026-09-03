# Migrate tglow-core (and tglow-pipeline/gamma) from aicsimageio to bioio

> Work for this plan happens on the `bioio-migration` branch of `tglow-core` (created off
> `main`), not directly on `main`. Note the package-root directory inside this repo is also
> named `main/` — that's an unrelated, pre-existing naming coincidence with the git branch.

## Context

`tglow-core` wraps `aicsimageio` in two classes — `AICSImageReader` / `AICSImageWriter` in
`src/tglow/io/tglow_io.py` — used to read/write the `/plate/row/col/field.ome.tiff` CZYX
layout that the whole tglow pipeline is built on. `aicsimageio` is now in maintenance mode and
its dependency chain (with `hyperactive`/`gradient-free-optimizers`/`pandas` pins for BaSiCPy) is
already flagged as a pain point in `README.md`. `bioio` is the actively maintained successor from
the same authors, with a near drop-in API (`AICSImage` → `BioImage`) but a split plugin
architecture (`bioio` + format-specific packages like `bioio-ome-tiff`).

Actual usage is narrower than "the whole package": only `tglow_io.py` touches `aicsimageio`
directly, only OME-TIFF is ever read/written (no `.czi`/`.lif`/`.nd2`/`.zarr` anywhere), and the
two downstream consumers (`compound_image_provider.py`, `processed_image_provider.py`) only touch
`.dims`, `.dtype`, `.channel_names` on the object `AICSImageReader.get_img()` returns — they never
import `aicsimageio` themselves. There is currently **no test suite** in `tglow-core`.

Per user decision, scope is:
- `tglow-core` (this repo) — full migration.
- `tglow-pipeline` — **`gamma/` only**. `gamma/bin/convert_pe_raw.py`, `gamma/bin/utils/downsample.py`,
  and `gamma/bin/rescale_images.py` import `aicsimageio` directly (not just via tglow-core) and
  must be updated so `gamma/` has zero residual `aicsimageio` dependency. `beta/` is left as-is
  (it will keep working as long as `aicsimageio` stays importable in that env — not our concern here).
- A minimal pytest safety net is added for `tglow_io.py` since none exists today.

## Step 1 — Add a characterization test harness (before touching the implementation)

No tests exist, so write them first against the *current* `aicsimageio`-backed code to lock in
behavior, then reuse unchanged after the swap to prove parity.

- New `tests/` package (currently excluded from packaging in `pyproject.toml` — that's fine,
  it just means tests aren't shipped in the wheel).
- `tests/conftest.py`: build a tiny synthetic `/plate/row/col/field.ome.tiff` fixture tree
  (e.g. 2 channels × 3 planes × 8×8) using whatever writer is live at test-run time, laid out to
  match `AICSImageReader.__build_index__`'s expected `plate/row/col/field.ome.tiff` structure.
- `tests/test_tglow_io.py`: cover `AICSImageReader.get_img/read_image/read_stack` (all four
  `query.channel`/`query.plane` combinations at `tglow_io.py:406-424`, including the
  channel=None/plane=not-None branch at line 419) and `AICSImageWriter.write_stack/write_image_stats`
  round-tripping a stack. Assert on `.dims['C'][0]`, `.dims.order`, `.channel_names`, `.dtype` — the
  exact attributes `compound_image_provider.py:47,50` and `processed_image_provider.py:185,189,239`
  rely on — so any behavioral drift in the new library is caught here, not downstream.
- Run once against current `aicsimageio` to confirm they pass (establishes the baseline).
- Add `pytest` to `lib/requirements.txt`/`pyproject.toml` as a dev dependency (or
  `[project.optional-dependencies].test`).

## Step 2 — Fix the latent bug while touching that code

`tglow_io.py:419`: `img.get_image_data("CYX", T=0, Z=int(query.plane)).compute()` — calls
`.compute()` on what `get_image_data` returns as a plain `numpy.ndarray` (only
`get_image_dask_data` returns a dask array with `.compute()`). This is a pre-existing bug (would
raise `AttributeError` if `query.channel is None and query.plane is not None` is ever hit) that
has no test coverage today. Since Step 1 adds a test for exactly this branch, fix it in the same
pass: drop the `.compute()` call.

## Step 3 — Swap tglow-core's dependency and imports

- `pyproject.toml:26` and `lib/requirements.txt:1`: replace `"aicsimageio>=4.14.0"` with
  `"bioio>=1.0"` and `"bioio-ome-tiff>=1.0"` (only reader/writer plugin actually needed — no
  `bioio-czi`/`bioio-lif`/etc since only OME-TIFF is used anywhere in this codebase).
- `pyproject.toml:10`: bump `requires-python = ">=3.9, <4"` → `">=3.10, <4"` —
  `bioio-ome-tiff` requires Python ≥3.10. Also bump the `3.10` classifier line if a floor changes.
- `src/tglow/io/tglow_io.py:21-23`:
  ```python
  from bioio import BioImage
  from bioio_ome_tiff.writers import OmeTiffWriter
  ```
  Drop the `OmeTiffReader` import entirely (`aicsimageio.readers.ome_tiff_reader.OmeTiffReader`
  at line 23) — it was already dead (only referenced in a commented-out line at 401).
- `tglow_io.py:389,402`: `AICSImage(...)` → `BioImage(...)`.
- `tglow_io.py:408,414,419,424`: `get_image_data`/`get_image_dask_data` calls are unchanged —
  bioio's `Reader` base class preserves these method names/signatures from aicsimageio v4.
- `tglow_io.py:614-620`: `OmeTiffWriter.save(...)` call is unchanged — bioio-ome-tiff's
  `OmeTiffWriter.save()` keeps the same signature (`dim_order`, `channel_names`,
  `physical_pixel_sizes`, `image_names`).
- Leave `AICSImageReader`/`AICSImageWriter` class names as-is — do **not** rename to
  `BioImageReader`/`BioImageWriter`. They're consumed by ~20 files across
  `tglow-pipeline/{beta,gamma}/bin/**`; renaming is pure churn with no functional benefit and
  would break `beta/` for nothing.
- Rerun Step 1's tests unchanged against bioio to confirm parity (this is the actual proof the
  migration is behavior-preserving).

## Step 4 — Update tglow-core docs

- `README.md`: remove the "Notes and migration to BioIO" section (lines 59-60) and the
  aicsimageio-specific sentence in "Known issues" (line 63) now that the migration is done; update
  line 10-11 ("wrappers around `aicsimageio`" → "wrappers around `bioio`"); keep both project links
  in "References" (lines 69-70) since bioio still descends from aicsimageio conceptually, or drop
  the aicsimageio link if preferred.
- Bump `pyproject.toml:7` version (e.g. `0.1.3` → `0.2.0`) — public API is unchanged but the
  runtime dependency and minimum Python version both changed, which is release-note-worthy for a
  published PyPI package.

## Step 5 — Fix direct aicsimageio imports in tglow-pipeline/gamma

Three files import `aicsimageio` directly (not just via tglow-core's wrapper), so they'll break
once tglow-core drops the dependency:

- `gamma/bin/convert_pe_raw.py:8`: `from aicsimageio.types import PhysicalPixelSizes` →
  `from bioio_base.types import PhysicalPixelSizes`.
- `gamma/bin/utils/downsample.py:9`: same change.
- `gamma/bin/rescale_images.py:11`: `from aicsimageio.writers import OmeTiffWriter` →
  `from bioio_ome_tiff.writers import OmeTiffWriter`.
- Check whatever environment/requirements file installs `gamma/`'s Python dependencies (no
  `requirements.txt`/Dockerfile was found under `tglow-pipeline` referencing `aicsimageio` — it's
  picked up transitively from `tglow-core`) and confirm `bioio-ome-tiff`/`bioio-base` end up
  installed there once `tglow-core` is bumped to the new version.

## Verification

1. `pip install -e ".[test]"` (or equivalent) in a Python ≥3.10 env, on the `bioio-migration`
   branch, then `pytest tests/` — confirm the Step 1 tests pass against the new bioio-backed
   implementation.
2. `grep -rn "aicsimageio" .` (this repo) and `grep -rn "aicsimageio" ../../tglow-pipeline/gamma/`
   — both should return zero hits when done.
3. If a real sample plate is available, run `AICSImageReader(...).read_stack(...)` and
   `AICSImageWriter(...).write_stack(...)` against it end-to-end once manually (not just the
   synthetic fixture) to sanity-check real OME-TIFF metadata (channel names, physical pixel sizes)
   round-trips correctly — this is the one thing a tiny synthetic fixture can't fully prove.
4. Smoke-run one `gamma/bin` script that hits each of the three directly-fixed files
   (`convert_pe_raw.py`, `downsample.py`, `rescale_images.py`) against a small input to confirm
   the import swap didn't break argument handling downstream of the import.
