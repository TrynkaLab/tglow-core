# Tglow-core

Main repo for the TeamTrynka imaging pipelines. Still very much a work in progress, more documentation to follow when things mature.

# Repo structure

- cellprofiler: cellprofiler pipelines for specific instances (projects/datasets)
- core: core component of the pipeline
    * runners: python / bash executor scripts used by the bash and nextflow pipelines
    * tglow: core python package for IO tasks
- nextflow: nextflow pipeline
- r-core: R functions for analyzing the output features
- templates: template scripts implementing the old bash based pipeline (2d)

# Workflows (nextflow)

1. stage:

This workflow takes a perkin elmer (currently only works for phenix) index xml and stages the files into a plate/row/col/field.ome.tiff file structure. IO is handled through AICSImageio, channel names and pixel sizes are extracted from the index.xml and set as metadata items in the ome tiff. Some additional files are staged for reproducabillity and ease of reading channel orders etc. This output is also intended to be backed up to iRODS

Proccess:
1. Create a manifest with the wells to run to re-use later as a nextflow channel
2. Read the files from /nfs and save directly into the above format


2. run_pipeline:

This takes as input the images produced during staging and then runs the following processes.

Processes:
1. Basicpy [optional]: parallelized on the plate + channel level, flatfields are saved, no images are stored
2. Registration (pystackreg) [optional]: parallelized on the well level (all fields). Only registration matrices are saved, no images are stored
3. Cellpose: parallelized on the well level (all fields), runs on GPU. If registering, currently only runs on the reference plate. Must run a cell segmentation, can optionally provide nucleus channel as well.
4. TBD deconvonvelution: Will spin out another copy of the data, runs on GPU
5. cellprofiler: 
    a. Stage the files into a cellprofiler compatable format and apply flatfields and registration if applicable
    b. run feature extraction


# Install instructions (nextflow)

This is not written to be very portable at the moment, especially the old bash pipeline is very sanger farm specific.
However the Nextflow pipeline should be more portable.

## High level dependencies

- Nextflow
- Conda

On farm22 can load through module system 

Due to conflicting versions we need two enviroments, one for the core tglow pipeline and one for cellprofiler
Check core/lib/requirements.txt for the dependencies for the tglow enviroment (this is installed by default)
Check core/lib/requirements_cellprofiler.txt for the cellprofiler enviroment 

### tglow enviroment

Dependencies (check requirements.txt)
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

conda create -n tglow python==3.10
conda activate tglow


Clone the repo into a suitable install dir, then:

git clone <this repo>
cd <this repo>/core
pip install .

To enable GPU install the GPU version of pytorch, following their instructions. I had more luck with the pip install
then the conda version

If you want to be able to edit the python package, install with:

pip install -e .

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

conda deactivate
conda deactivate

Create the cellprofiler environment:

conda create -n cellprofiler python==3.8

open core/setup.py and edit the following line:

requirement_path = f"{lib_folder}/lib/requirements.txt"
to
requirement_path = f"{lib_folder}/lib/requirements_cellprofiler.txt"

then run:

pip install ./

Keep track of the install path of the conda enviroment, as this will need to be provided to nextflow


# Running the pipeline (bash)

Make a copy of the templates to the location where you want to run, and edit each of the scripts so they are configured properly to your data.

# Running the pipline (nextflow)

The nextflow pipline runs in two stages. Some of it is not done through fully "proper" nextflow, as nextflow is very storage heavy and can easily store redundant copies of data which is not great for imaging. 

Hence most processes rely on the storeDir option which serves as a permanent cache between runs. Keep in mind that as long as nextflow finds the files in this storeDir, processes are NOT re-run, so you will have to remove them manually if you want to re-do some steps! This is the case for the stageing of the raw images, the deconvolution, registration and basicpy. This is not "proper" nextflow, but it made the most sense in this case. Many of the processes rely on eachother yet don't have a natural  a > b > c structure but represent a complicated tree. 

Parallelization is done on the per well level to not overload the system with 20 second jobs wasting a lot of time on overheads (loading conda, initializing python etc). This makes some of the processes implicitly assume all the fields are supposed to be run. Given the same well always has the same number of fields this should work fine. But nextflow itself is not aware of the fields! This pipeline doesn't currently support seperate fields between cycles, but this should ideally never happen anyway. If it does happen, first stage the data, then remove or rename the fields so they match manually in the storeDir, then it should work fine.

The first workflow involves staging the data from an Harmony export on NFS. This is done with the -entry stage

nextflow <path/to/tglow.nf> -entry stage <params>

There are a bunch of configurations to be set, have a look at the nextflow.config for detaills. In general, you make a manifest in the following form

<plate> <index.xml> <channels> <basicpy_channels> <cellpose_nucleus> <cellpose_cell>

that tells the pipeline where the data lives, and which channels are which and how to apply them. The file can run multiple plates at the same time by adding more lines.


# TODO items:

- stage image folders from iRODS bacup instead of PE export
- 

