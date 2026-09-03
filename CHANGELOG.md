

# 0.2.0

- Migrated `AICSImageReader` / `AICSImageWriter` from `aicsimageio` to `bioio` + `bioio-ome-tiff`. Public API (class names, method signatures) is unchanged.
- Bumped minimum supported Python to 3.10 (required by `bioio-ome-tiff`).
- Fixed a latent bug in `AICSImageReader.read_image` where the channel=None/plane=not-None code path called `.compute()` on a plain numpy array.
- Added a pytest suite covering `AICSImageReader`/`AICSImageWriter` (previously untested).


# 0.1.4
Should be backwards compatible with 0.1.3

- Added util functions for sigmoid and rescale_stack, rescale_stack_inplace
- ProccessedImageProvider now uses rescale_stack_inplace
- Removed dependency for basicpy, and replaced with a reader that loads the npz file direct


# Prior to 0.1.3

Did not keep a good record 