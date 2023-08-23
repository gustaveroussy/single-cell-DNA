import os
import sys
import missionbio.mosaic as ms
import subprocess
import missionbio.mosaic.utils as mutils
from missionbio.mosaic.constants import COLORS
import missionbio.mosaic.io as mio
import plotly
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import plotly.io as pio
from pathlib import Path
from yaml.loader import SafeLoader
import yaml
import seaborn as sns
import pandas as pd
import numpy as np
import re
from pathlib import Path
import os
import ast
import argparse


parser = argparse.ArgumentParser(description='single-cell DNA-seq protein normalization (CLR,asinh,NSP)')

parser.add_argument("--input_h5", help="input h5 path file",required=True,action='store')
parser.add_argument("--output_h5", help="output h5 path file",required=True,action='store')
parser.add_argument("--sample_name", help="wilcard sample",required=True,action='store')
parser.add_argument("--config_file", help="config file",required=True,action='store')


os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

lib_to_import=sys.argv[0].replace("/scripts/prot_normalization.py", "/")
sys.path.append(lib_to_import+"/common")
from utils import *

args = parser.parse_args()

config=args.config_file
config=ast.literal_eval(config)


try:
    norm_meth=config["prot_norm_dimred"]["normalization"]
except KeyError:
    norm_meth=["CLR"]

sample_path = args.input_h5
sample = ms.load(sample_path, raw=False, apply_filter=True, single=True)

for i_method_norm in norm_meth:
        if i_method_norm in ["CLR", "NSP", "asinh"]:
                sample.protein.normalize_reads(i_method_norm)
                sample.protein.layers['normalized_counts_'+str(i_method_norm)]=sample.protein.layers['normalized_counts']

del sample.protein.layers['normalized_counts']

mio.save(sample=sample,path=args.output_h5,raw=False)

path = Path(config["output_sample_path"]+"/"+args.sample_name+"/prot/diagnostic/")
path.mkdir(parents=True, exist_ok=True)

for i_method_norm in norm_meth:
        fig=sample.protein.ridgeplot(attribute='normalized_counts_'+i_method_norm,features=sample.protein.ids())
        fig.write_html(config["output_sample_path"]+"/"+args.sample_name+"/prot/diagnostic/"+"ridgeplot_norm_"+i_method_norm+".html")