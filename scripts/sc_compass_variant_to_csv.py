import os
import sys
import missionbio.mosaic as ms
from yaml.loader import SafeLoader
import yaml
import os
from pathlib import Path
import re
import ast

lib_to_import=sys.argv[0].replace("/scripts/sc_compass_variant_to_csv.py", "/")
sys.path.append(lib_to_import+"common")
from utils import *

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

config=sys.argv[4]
config=ast.literal_eval(config)

sample_path = sys.argv[1]

sample = ms.load(sample_path, raw=False, apply_filter=True, single=True)

annotation_var=pd.read_csv(sys.argv[2])

sample.dna= sample.dna[sample.dna.barcodes(),annotation_var["Variant ID"].values]

sample_name= [sample.dna.row_attrs['sample_name'][0]]*len(sample.dna.col_attrs['CHROM'])

data = {'sample ID': sample_name,
       'chr':sample.dna.col_attrs['CHROM'],
       'start':sample.dna.col_attrs['POS'],
       'ref allele':sample.dna.col_attrs['REF'],
       'alt allele':sample.dna.col_attrs['ALT']}

df = pd.DataFrame(data)

df.to_csv(sys.argv[3],index=False)