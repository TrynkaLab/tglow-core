#!/bin/bash
# Usage
# run_companion_ome.sh <batch_file> <stitched input dir>
# batch file should be a list of plate names, one per line

# Configure bash to work with conda
source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh

# Important so the correct python version is loaded
conda deactivate

# Activate conda env
conda activate /software/teamtrynka/cellprofiler

PIPELINE_DIR="/software/teamtrynka/tglow-core"
BATCH="batch_files/ki67_test_plates.files"
INPUT_DIR="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/KI67_TEST/1_pre_process/output/stitched"

while read batch;
do

path="${INPUT_DIR}/${batch}/index.xml"

CMD="python ${PIPELINE_DIR}/core/companion.py \
--input_file ${path}"

echo $CMD
eval $CMD

done < $BATCH
