# Tglow-core

Main repo for the TeamTrynka imaging pipelines. Still very much a work in progress, more documentation to follow when things mature.

Likely will re-organize this in the future into tglow-core and tglow-feature, with tglow-core being dataset invariant and tglow-feature being dataset specific scripts for feature extraction. It all depends a bit how flexible we make the codebase and how much effort we put into proper instancing over "hacking" a new instance. 


# Repo structure

- cellprofiler: cellprofiler pipelines for specific instances (projects/datasets)
- core: core component of the pipeline
    * runners: python / bash executor scripts used by the bash and nextflow pipelines
    * tglow: core python package for IO tasks
- nextflow: nextflow pipeline
- r-core: R functions for analyzing the output features
- templates: template scripts implementing the old bash based pipeline (2d)


# Install instructions (nextflow)

This is not written to be very portable at the moment, especially the old bash pipeline is very sanger farm specific.
However the Nextflow pipeline should be more portable.

## High level dependencies

- Nextflow
- Conda

On farm22 can load through module system 

## Setup conda enviroments

Due to conflicting versions we need two enviroments, one for the core tglow pipeline and one for cellprofiler

### tglow enviroment
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

Once those are installed install the tglow core component

git clone <this repo>
cd <this repo>/core
pip install .

If you want to be able to edit the python package, install with:

pip install -e .

### Cellprofiler env

- python 3.8
- cellprofiler
- cellpose 2 (used as a cp plugin for brightfield images, at this time cellpose 3 doesn't work)
- tiffifle
- pystackreg 
- numpy
- skimage
- matplotlib


# Running the pipeline (bash)

Make a copy of the templates to the location where you want to run, and edit each of the scripts so they are configured properly to your data.

# Running the pipline (nextflow)

The nextflow pipline runs in two stages. Some of it is not done through fully "proper" nextflow, as nextflow is very storage heavy and can easily store redundant copies of data which is not great for imaging. 

Hence most processes rely on the storeDir option which serves as a permanent cache between runs. Keep in mind that as long as nextflow finds the files in this storeDir, processes are NOT re-run, so you will have to remove them manually if you want to re-do some steps! This is the case for the stageing of the raw images, the deconvolution, registration and basicpy. This is not proper nextflow, but it made the most sense in this case. Many of the processes rely on eachother yet don't have a natural  a > b > c structure but represent a complicated tree. 

Parallelization is done on the per well level to not overload the system with 20 second jobs wasting a lot of time on overheads (loading conda, initializing python etc). This makes some of the processes implicitly assume all the fields are supposed to be run. Given the same well always has the same number of fields this should work fine. But nextflow itself is not aware of the fields! This pipeline doesn't currently support seperate fields between cycles, but this should ideally never happen anyway. If it does happen, first stage the data, then remove or rename the fields so they match manually in the storeDir, then it should work fine.

The first workflow involves staging the data from an Harmony export on NFS. This is done with the -entry stage

nextflow <path/to/tglow.nf> -entry stage <params>

There are a bunch of configurations to be set, have a look at the nextflow.config for detaills. In general, you make a manifest in the following form

<plate> <index.xml> <channels> <basicpy_channels> <cellpose_nucleus> <cellpose_cell>

that tells the pipeline where the data lives, and which channels are which and how to apply them. The file can run multiple plates at the same time by adding more lines.
