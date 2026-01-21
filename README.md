# Tglow: Core Python component of the tglow imaging pipeline

`tglow-core` is the Python core component of the Tglow high-content imaging (HCI)
analysis pipeline. It provides utilities to index and read multi-well plate images
and parsers for PerkinElmer (Opera Phenix / Operetta) exports. The package is
used by the `tglow-pipeline` workflows to load, preprocess and write OME-TIFF
images arranged in the common `/plate/row/col/field.ome.tiff` (CYZX) layout.

Key features
- Read and write CYZX / ZYX / YX image arrays via `AICSImageReader` /
	`AICSImageWriter` (wrappers around `aicsimageio`)
- Parse Revity/PerkinElmer `Index.xml` exports (`PerkinElmerParser`) and convert to a simple, Python-friendly index
- Convert large Revity/PerkinElmer exports to a much lower number of `/plate/row/col/field.ome.tiff` files
- Index and query plate/row/col/field image layouts using an `ImageQuery` object
- Utilities for registration, flatfield correction and numeric conversions designed to work with `tglow-pipeline`


# Installation
I recommend installing the published PyPI release where possible:

```bash
pip install tglow-core
```

To install the latest development version from the repository (editable install):

```bash
git clone https://github.com/TrynkaLab/tglow-core
cd tglow-core
pip install -e .
```

# Basic usage

Build an index from a PerkinElmer export and read a single image:

```python
from tglow.io.tglow_io import PerkinElmerRawReader
from tglow.io.image_query import ImageQuery

reader = PerkinElmerRawReader('path/to/Index.xml', '/data/exports')
iq = ImageQuery.from_plate_well('plate1', 'A01')
image = reader.read_image(iq)  # returns a numpy array
```

Read and write an OME-TIFF stack organized by plate/row/col/field:

```python
from tglow.io.tglow_io import AICSImageReader, AICSImageWriter
from tglow.io.image_query import ImageQuery

reader = AICSImageReader('/data/plates')
writer = AICSImageWriter('/output/plates')
iq = ImageQuery('plate1', 1, 1, 'field001')
stack = reader.read_stack(iq)
writer.write_stack(stack, iq)
```

# Notes and migration to BioIO
- This package currently wraps `aicsimageio`. As that project has been superseded by newer tooling, consider migrating to `bioio` or equivalent in future releases.

# Known issues
There is a known issue with BaSiCpy (https://github.com/peng-lab/BaSiCPy/issues/162). This requires using specific, older versions of `hyperactive` and `gradient-free-optimizers`, which can in turn require an older `pandas` version. The same goes for `aicsimageio`. The dependency chain can be inconvenient; I'll update this project and migrate to BioIO libraries as newer releases become available.

# Acknowledgements
- Martin Prete: initial XML parsing code adapted for this project

# References
- https://github.com/AllenCellModeling/aicsimageio
- https://github.com/bioio-devs/bioio
