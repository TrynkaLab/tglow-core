

# 0.2.0

- Migrated `AICSImageReader` / `AICSImageWriter` from `aicsimageio` to `bioio` + `bioio-ome-tiff`. Public API (class names, method signatures) is unchanged.
- Bumped minimum supported Python to 3.10 (required by `bioio-ome-tiff`).
- Fixed a latent bug in `AICSImageReader.read_image` where the channel=None/plane=not-None code path called `.compute()` on a plain numpy array.
- Added a pytest suite covering `AICSImageReader`/`AICSImageWriter` (previously untested).
- Added `bioio-tifffile` as a dependency - `bioio-ome-tiff` hard-requires valid OME-XML metadata, so it rejects plain TIFFs with none at all (e.g. Cellpose's `*_cp_masks.tiff` mask output), which broke `mask_reader` in `ProcessedImageProvider` after the migration above. `BioImage`'s own plugin auto-detection already falls back to any other registered plugin when one rejects a file, so this fixes mask reading with no code change - `bioio-ome-tiff` is still preferred (and unaffected) for real OME-TIFF images.
- Added rendering for QC report in the tglow pipeline

# 0.1.4
Should be backwards compatible with 0.1.3

- Added util functions for sigmoid and rescale_stack, rescale_stack_inplace
- ProccessedImageProvider now uses rescale_stack_inplace
- Removed dependency for basicpy, and replaced with a reader that loads the npz file direct


# Prior to 0.1.3

Did not keep a good record 