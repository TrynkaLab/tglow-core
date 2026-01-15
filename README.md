# Tglow: Core python component of tglow imaging pipeline.

`tglow-core` is the Python core component of the Tglow high-content imaging (HCI)
analysis pipeline. It provides utilities to index and read multi-well plate images
and parsers for PerkinElmer (Opera Phenix / Operetta) exports. The package is
used by the `tglow-pipeline` workflows to load, preprocess and write OME-TIFF
images arranged in the common `/plate/row/col/field.ome.tiff` layout.

Key features
- Index plate/row/col/field image layouts and return an `ImageQuery` object
- Read and write CZYX / ZYX / YX image arrays via `AICSImageReader` /
	`AICSImageWriter` (wrappers around `aicsimageio`)
- Parse PerkinElmer `Index.xml` exports (`PerkinElmerParser`) and convert to a
	simple Python-friendly index
- Utilities for registration, flatfield correction and numeric conversions



# Installation


```bash
conda create -n tglow python==3.10 
conda activate tglow
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
- This package currently wraps `aicsimageio`. As that project has been
	superseded by newer tooling, consider migrating to `bioio` or equivalent in
	future releases.

# Acknowledgements
- Martin Prete: initial XML parsing code adapted for this project

# References
- https://github.com/AllenCellModeling/aicsimageio
- https://github.com/bioio-devs/bioio
