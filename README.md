# Tglow: Core python component of tglow imaging pipeline.

Core python package for the python component of the tglow HCI analysis pipeline. Mainly contains wrappers arround AICSImageIO to index folders in /plate/row/col/field.ome.tiff and contains code to read and write Revity Opera Phenix and Operetta XML files. By itself its not very usefull but its a dependency for the tglow pipeline

# Install instructions

Create a new conda enviroment

```
conda create -n tglow python==3.10
conda activate tglow
```

Clone the repo into a suitable install dir, then:

```
git clone https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-core
cd tglow-core
pip install -e .
```
