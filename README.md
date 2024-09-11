# Tglow: Core python component of tglow imaging pipeline.

Core repo for the python component of the tglow HCI analysis pipeline.

<br>
<br>

# Install instructions

Discalimer, this is not written to be very portable at the moment as it is using conda as package manager, and uses some configurable paths for custom scripts. Will containerize this at some point. Also the GPU installations might require some manual configuring within conda to get them to wrok.

## High level dependencies

- Nextflow
- Conda

On farm22 can load through module system 

Due to conflicting versions we need two enviroments, one for the core tglow pipeline and one for cellprofiler
Check core/lib/requirements.txt for the dependencies for the tglow enviroment (this is installed by default)
Check core/lib/requirements_cellprofiler.txt for the cellprofiler enviroment 

### tglow enviroment

Dependencies (check core/lib/requirements.txt)
- python 3.10
- cellpose 3
- cuda 12
- pytorch (gpu enabled)
- basicpy
- pystackreg
- AICSImageio
- numpy
- skimage
- matplotlib
- scipy

Create a new conda enviroment

```
conda create -n tglow python==3.10
conda activate tglow
```

Clone the repo into a suitable install dir, then:

```
git clone https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-core
cd tglow-core/core
pip install .
```

To enable GPU install the GPU version of pytorch, following their instructions. I had more luck with the pip install
then the conda version.

IF you want to be able to edit the python package without the need to re-install, install with:

```
pip install -e .
```
Keep track of the install path of the conda enviroment, as this will need to be provided to nextflow


### Cellprofiler env

Dependencies (check requirements_cellprofiler.txt)
- python 3.8
- cellprofiler
- cellpose 2 (used as a cp plugin for brightfield images, at this time cellpose 3 doesn't work)
- tiffifle
- pystackreg 
- numpy
- skimage
- matplotlib


Deactivate the previous conda environments (twice to make sure)
```
conda deactivate
conda deactivate
```
Create the cellprofiler environment:
```
conda create -n cellprofiler python==3.8
```
open core/setup.py and edit the following line:

```
requirement_path = f"{lib_folder}/lib/requirements.txt"
```
to

```
requirement_path = f"{lib_folder}/lib/requirements_cellprofiler.txt"
```
then run:
```
pip install ./
```
Keep track of the install paths of the conda enviroments, as this will need to be provided to nextflow




