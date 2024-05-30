# Tglow-core

Main repo for the TeamTrynka imaging pipelines. Still very much a work in progress, more documentation to follow when things mature.

<br>
<br>

# Repo structure

- cellprofiler: cellprofiler pipelines for specific instances (projects/datasets)
- core: core component of the pipeline
    * runners: python / bash executor scripts used by the bash and nextflow pipelines
    * tglow: core python package for IO tasks
- nextflow: nextflow pipeline
- r-core: R functions for analyzing the output features
- templates: template scripts implementing the old bash based pipeline (2d)

<br>
<br>

# Workflows (nextflow)

## 1. stage:

This workflow takes a perkin elmer (currently only works for phenix) index xml and stages the files into a plate/row/col/field.ome.tiff file structure. IO is handled through AICSImageio, channel names and pixel sizes are extracted from the index.xml and set as metadata items in the ome tiff. Some additional files are staged for reproducabillity and ease of reading channel orders etc. This output is also intended to be backed up to iRODS

Proccess:
1. Create a manifest with the wells to run to re-use later as a nextflow channel
2. Read the files from /nfs and save directly into the above format

## 2. run_pipeline:

This takes as input the images produced during staging and then runs the following processes.

Processes:
1. Basicpy [optional]: parallelized on the plate + channel level, flatfields are saved, no images are stored
2. Registration (pystackreg) [optional]: parallelized on the well level (all fields). Only registration matrices are saved, no images are stored
3. Cellpose: parallelized on the well level (all fields), runs on GPU. If registering, currently only runs on the reference plate. Must run a cell segmentation, can optionally provide nucleus channel as well.
4. TBD deconvonvelution: Will spin out another copy of the data, runs on GPU
5. feature extraction / cellprofiler: 
    1. Stage the files into a cellprofiler compatable format, apply flatfields, bring together imaging cylces and apply registration if applicable
    2. run feature extraction (cellprofiler) or other custom script

<br>
<br>


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
git clone <this repo>
cd <this repo>/core
pip install .
```

To enable GPU install the GPU version of pytorch, following their instructions. I had more luck with the pip install
then the conda version

If you want to be able to edit the python package, install with:

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
Keep track of the install path of the conda enviroment, as this will need to be provided to nextflow


# Running the pipeline (bash)

Make a copy of the templates to the location where you want to run, and edit each of the scripts so they are configured properly to your data.

# Running the pipline (nextflow)

## Quick explanation

Install as indicated above

Adjust the required configurations either in a config file or on the nextflow command (see nettflow docs for more detaills). See nextflow/nextflow.config for available parameters.

prepare a manifest.tsv with one line per plate:
```
plate	index_xml	channels	bp_channels	cp_nucl_channel cp_cell_channel
ref_plate   index.xml   1,2,4,5 3,5 5   2  
qry_plate1  index.xml  2,3,4 none   none    none
qry_plate2  index.xml  3,4 3,4 3   4  

```
- plate: The exact plate name in the export
- index_xml: path the perkin elmer (or Martin's version) index xml file
- channels: list of available channels
- bp_channels: channels to run basicpy for
- cp_nucl_channel:optional channel for nuclei segmentation
- cp_cell_channel: cellpose channel for cell segmentation


First stage the data so its accesible on the lustre

```
OUTPUT=./output

nextflow run tglow.nf \
-profile lsf \
-w ${OUTPUT}/workdir \
-resume \
-entry stage \
--rn_manifest manifest.tsv \
--rn_publish_dir ${OUTPUT}/results \
--rn_image_dir ${OUTPUT}/results/images
```

After the data is staged into rn_image_dir run the pipeline
```
OUTPUT=./output

nextflow run tglow.nf \
-profile lsf \
-w ${OUTPUT}/workdir \
-resume \
-entry run_pipeline \
--rn_manifest manifest.tsv \
--rn_publish_dir ${OUTPUT}/results \
--rn_image_dir ${OUTPUT}/results/images
```

## Registering multiple cycles of imaging

To register imaging cycles additionally provide a registration manifest which links plates together. This should have the form:

```
reference_plate	reference_channel	query_plates	query_channels
ref_plate   5   qry_plate1,qry_plate2   3,4
```

- reference_plate: Plate name of reference plate
- reference_channel: 1 based index of channel to use in registering (nucleus)
- query_plates: comma seperated list of plate names to register against reference
- query_channels: comma seperated list of 1 based channel indices of channel to register (nucleus)

If this is provided, in downstream tasks (cellprofiler/feature extraction) data is treated as one plate and the query plates are treated as extra channels, with their indices increasing seqeuntially in the order specified in the manifest. 

Query plates must be provided in the manifest.tsv!
Currently cellpose is only run on reference plates, even if the channels are provided in the registration manifest. If needed, will update cellpose to run on the registered data so other cycle channels can be used in segmenting.

## Long explanation
The nextflow pipline runs in two stages. Some of it is not done through fully "proper" nextflow, as nextflow is very storage heavy and can easily store redundant copies of data which is not great for imaging. 

Hence most processes rely on the storeDir option which serves as a permanent cache between runs. Keep in mind that as long as nextflow finds the files in this storeDir, processes are NOT re-run, so you will have to remove them manually if you want to re-do some steps! This is the case for the stageing of the raw images, the deconvolution, registration and basicpy. This is not "proper" nextflow, but it made the most sense in this case. Many of the processes rely on eachother yet don't have a natural  a > b > c structure but represent a complicated tree. 

Parallelization is done on the per well level to not overload the system with 20 second jobs wasting a lot of time on overheads (loading conda, initializing python etc). This makes some of the processes implicitly assume all the fields are supposed to be run. Given the same well always has the same number of fields this should work fine. But nextflow itself is not aware of the fields! This pipeline doesn't currently support seperate fields between cycles, but this should ideally never happen anyway. If it does happen, first stage the data, then remove or rename the fields so they match manually in the storeDir, then it should work fine.

The first workflow involves staging the data from an Harmony export on NFS. This is done with the -entry stage

```
nextflow <path/to/tglow.nf> -entry stage <params>
```

There are a bunch of configurations to be set, have a look at the nextflow.config for detaills. In general, you make a manifest in the following form

```
<plate> <index.xml> <channels> <basicpy_channels> <cellpose_nucleus> <cellpose_cell>
```
that tells the pipeline where the data lives, and which channels are which and how to apply them. The file can run multiple plates at the same time by adding more lines.


# TODO items:

- stage image folders from iRODS bacup instead of PE export
- Deconvolution
- Detailled description of parameters

