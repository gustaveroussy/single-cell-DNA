import os
import functools
from snakemake.utils import min_version
from pathlib import Path
import yaml
import glob
from typing import Any, Dict, List
min_version("7.5")
import sys

### Importing Configuration Files --- ###

PIPELINE_FOLDER=Path(snakemake.workflow.srcdir(""))
PIPELINE_COMPASS_PREPROCESS=str(PIPELINE_FOLDER)+"/envs/conda/preprocess.yaml"
#SING_IMG="/mnt/beegfs/pipelines/single-cell_dna/1.0"+"/envs/singularity/SCDNA_compass.simg"
SING_IMG_DSB="/mnt/beegfs/pipelines/single-cell_dna/1.1"+"/envs/singularity/single_cell_dna_dsb.simg"
CONDA_MOSAIC_ENV=str(PIPELINE_FOLDER)+"/envs/conda/environment_mosaic_droplet.yml"
CONDA_DSB_ENV=str(PIPELINE_FOLDER)+"/envs/conda/environement_dsb.yaml"
CONFIG_FILE_PATH=sys.argv[6]
SING_IMG=str(PIPELINE_FOLDER)+"/envs/singularity/tree_building_container.simg"

lib_to_import=str(PIPELINE_FOLDER)+"/common/" 
sys.path.append(lib_to_import)
from utils import *

GLOBAL_TMP = config['tmp'] if 'tmp' in config else "/tmp"
if os.path.normpath(GLOBAL_TMP) != "/tmp" :
    if os.path.exists(GLOBAL_TMP) :
        sys.stderr.write("Temporary directory is set to: " + GLOBAL_TMP + "\n")
    else :
        sys.stderr.write(GLOBAL_TMP + " doesn't exist! Temporary directory is set to /tmp \n")
        GLOBAL_TMP = "/tmp"

sample_list=[sample for sample in config["sample"] if sample in os.listdir(config["input_sample_path"])]

i_steps=config["steps"]

include: "rules/rule_all.smk"

rule all:
    input:
        **get_targets(i_steps)
    message:
        "pipeline goes all way long through your step(s)"

if "Alignment" in i_steps:
    include: "rules/alignment.smk"

if "filtering" in i_steps:
    include: "rules/preprocessing.smk"

if "SNV_CNV" in i_steps:
    include: "rules/SNV_CNV.smk"
    
if "PROTEIN" in i_steps:
    include: "rules/PROTEIN.smk"
    
if "ALL" in i_steps:
    include: "rules/ALL.smk"
    
if "phylogeny" in i_steps:

    include: "rules/phylogeny.smk"