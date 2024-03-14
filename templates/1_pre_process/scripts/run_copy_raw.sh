#!/bin/bash
# Usage
# copy_raw.py authored by Martin Prete
# Also edited the job names to make them easier to parse
####################################################################
# Script invocation:
# python actions/copy_raw.py \
#       --input_file="/path/to/raw/export/Images/Index.xml" \
#       --output_path="/lustre/scratch/raw" 
####################################################################

# Configure bash to work with conda
source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/basicpy

# Setup variables
PIPELINE_DIR="/software/teamtrynka/tglow-core"

python ${PIPELINE_DIR}/core/copy_raw.py $0
