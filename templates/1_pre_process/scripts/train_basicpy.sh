#!/bin/bash
# Usage
# train_basicpy.sh <output_root> <channel>

# Configure bash to work with conda
source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/conda/basicpy

output=$1
channel=$2

# Setup variables
PIPELINE_DIR="/software/teamtrynka/tglow-core"
IMG_RAW="../output/raw"

# Run script
python ${PIPELINE_DIR}/core/train_basicpy.py \
--input ${IMG_RAW} \
--output ${output} \
--output_prefix ch${channel}/whole_plate_nimg200_mergen100 \
--channel ${channel} \
--nimg 200 \
--merge_n 100