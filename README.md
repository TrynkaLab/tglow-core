# Tglow: Core python component of tglow imaging pipeline.

Core package for the python component of the tglow HCI analysis pipeline.

# Install instructions

Discalimer, this is not written to be very portable at the moment as it is using conda as package manager, and uses some configurable paths for custom scripts. Will try to containerize this at some point. Also the GPU installations might require some manual configuring within conda to get them to wrok.

## High level dependencies

- Nextflow
- Conda

On farm22 can load through module system 

Due to conflicting versions we need two enviroments, one for the core tglow pipeline and one for cellprofiler
Check core/lib/requirements.txt for the dependencies for the tglow enviroment (this is installed by default)
Check core/lib/requirements_cellprofiler.txt for the cellprofiler enviroment 

### Environment 1: tglow enviroment

High level dependencies (check core/lib/requirements.txt for python dependencies installed during package install)
- python 3.10
- cuda 12
- pytorch (gpu enabled)
- opencl (should be installed with RedLionFish and clij2-fft, but your milage may vary)
- pyopencl (should be installed with RedLionFish and clij2-fft, but your milage may vary)

These are in requirements.txt, but might need some fiddling to get to work on GPU
- clij2-fft (https://pypi.org/project/clij2-fft/, https://github.com/clij/clij2-fft)
- RedLionFish (https://pypi.org/project/RedLionfish/)

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

To enable GPU install the GPU version of pytorch (https://pytorch.org/get-started/locally/), following their instructions. 
I had more luck with the pip install then the conda version. This is required if you want to run deconvolution using CLIJ2
RedLionFish or use Cellpose with GPU. 

Keep track of the install path of the conda enviroment, as this will need to be provided to nextflow.

### Environment 2: Cellprofiler env

High level dependencies (check requirements_cellprofiler.txt for python dependencies installed during package install)
- python 3.9
- cellprofiler

Deactivate the previous conda environments (twice to make sure)
```
conda deactivate
conda deactivate
```
Create the cellprofiler environment:
```
conda create -n cellprofiler python==3.9
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
pip install -e ./
```

Then install cellprofiler using the pip install provided instructions:
https://github.com/CellProfiler/CellProfiler/wiki/Source-installation-(Linux)


Keep track of the install paths of the conda enviroments, as this will need to be provided to nextflow.





