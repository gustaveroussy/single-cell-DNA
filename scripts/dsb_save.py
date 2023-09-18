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
from scipy import stats

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

lib_to_import=sys.argv[0].replace("/scripts/dsb_save.py", "/")
sys.path.append(lib_to_import+"common")
from utils import *

config=sys.argv[-1]
config=ast.literal_eval(config)

try:
        method_list=config["prot_norm_dimred"]["normalization"]
except KeyError:
        method_list=["CLR"]

if "DSB" in method_list:
    sample_path = sys.argv[1]
    sample = ms.load(sample_path, raw=False, apply_filter=True, single=True)
    
    df = pd.read_csv(sys.argv[2], sep=',')
    df.rename(columns = {'Unnamed: 0':'id'}, inplace = True)
    df.index=df['id'].values
    df.drop('id', inplace=True, axis=1)
    df=df.T
    
    if "Fc..RI.." in df.columns:
        
        tmp_ids=[i_string.replace('Fc..RI..','FcεRIα') for i_string in np.array(df.columns)]
        id_prot=[i_string.replace('.','-') for i_string in tmp_ids]
        sample.protein = sample.protein[np.array(df.index),id_prot]
        sample.protein.col_attrs['id']=np.array(id_prot)
        
    else:
        
        sample.protein = sample.protein[np.array(df.index),sample.protein.ids()]
        sample.protein.layers['normalized_counts_DSB']=df.values
        df_z = df.select_dtypes(include='number').apply(stats.zscore)
        sample.protein.layers['normalized_counts_DSB_znorm']=df_z.values
    

mio.save(sample,path=sys.argv[3],raw=False)