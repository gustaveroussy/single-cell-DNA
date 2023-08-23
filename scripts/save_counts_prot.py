import os
import sys
import missionbio.mosaic as ms
import subprocess
import missionbio.mosaic.io as mio
from pathlib import Path
from yaml.loader import SafeLoader
import yaml
import pandas as pd
import numpy as np
import re
from pathlib import Path
import os
import ast

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

lib_to_import=re.sub(r"/\w*.py","/",sys.argv[0])
sys.path.append(lib_to_import)
from function_utils import *

config=sys.argv[-1]
config=ast.literal_eval(config)

sample_path = sys.argv[1]
sample = ms.load(sample_path, raw=False, apply_filter=True, single=True)

try:
        method_list=config["prot_norm_dimred"]["normalization"]
except KeyError:
        method_list=["CLR"]

if "DSB" in method_list:
    sample.protein.get_attribute('read_counts',constraint='row+col').to_csv(sys.argv[2], sep="\t")