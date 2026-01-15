# Tglow: Core python component of tglow imaging pipeline.

Core python package for the python component of the tglow HCI analysis pipeline. Mainly contains wrappers arround AICSImageIO to index folders in /plate/row/col/field.ome.tiff and contains code to read and write Revity Opera Phenix and Operetta XML files. It is a core dependency of the tglow-pipeline. 

NOTE: This uses the aicsimageio which has been EOL since December 2025. We plan to migrate to its successor bioIO in the near future.


# Usage

By itself it can be used to interface with raw Opera Phenix or Opera Operetta exports with the class tglow.io.PerkinElmerRawReader which takes an exports Index.xml or Index.xml.idx. /plate/row/col/field.ome.tiff can be read and written using tglow.io.AICSImageReader and and tglow.io.AICSImageWriter. The class tglow.io.ImageQuery is a bean used to store plate/row/col/field/channel/plane information and can be passed to AICSImageReader.read_stack(), AICSImageReader.read_image() and AICSImageWriter.write_stack().


# Install instructions

Create a new conda enviroment

```
conda create -n tglow python==3.10
conda activate tglow
```

Clone the repo into a suitable install dir, then:

```
git clone https://github.com/TrynkaLab/tglow-core
cd tglow-core
pip install -e .
```

# Acknowledgements
We thank Martin Prete who supplied an inital version of the XML parsing script on which we iterated.


# References
- https://github.com/AllenCellModeling/aicsimageio
- https://github.com/bioio-devs/bioio
