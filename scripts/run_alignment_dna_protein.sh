#!/bin/bash

source /mnt/beegfs/pipelines/single-cell_dna/tapestri_pipeline/v2/etc/profile.d/conda.sh
conda activate /mnt/beegfs/pipelines/single-cell_dna/tapestri_pipeline/v2
tapestri dna+protein run --output-folder $1 --config $2 --overwrite
conda deactivate
