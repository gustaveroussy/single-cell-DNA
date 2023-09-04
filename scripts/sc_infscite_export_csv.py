import os
import sys
import missionbio.mosaic as ms
from yaml.loader import SafeLoader
import yaml
import os
from pathlib import Path
import re
import ast
import argparse


parser = argparse.ArgumentParser(description='infSCITE preprocessing')

parser.add_argument("--input_h5", help="input h5 path file",required=True,action='store')
parser.add_argument("--ouput_path", help="outputpath file",required=True,action='store')
parser.add_argument("--sample_name", help="sample_name",required=True,action='store')

args = parser.parse_args()

lib_to_import=sys.argv[0].replace("/scripts/sc_infscite_export_csv.py", "/")
sys.path.append(lib_to_import+"/common")
from utils import *

sample = ms.load(args.input_h5, raw=False, apply_filter=True, single=True)

#try:
#    mutation_matrix_type=config["infSCITE"]["mutation_matrix_type"]
#except KeyError:
#    mutation_matrix_type=["default"]

Path(args.ouput_path).mkdir(parents=True, exist_ok=True)

df=pd.DataFrame(sample.dna.layers["NGT"]).T

Path(args.ouput_path).mkdir(parents=True, exist_ok=True)

df.to_csv(args.ouput_path+"/mutation_matrix.csv",index=False,header=False,sep=" ")

parameter_df=pd.DataFrame({"ado":[sample.dna.metadata["ado_rate"]],"fdr":[0.001],"doublet_rate":[0.08],"n_sample":[len(sample.dna.barcodes())],"n_variants":[len(sample.dna.ids())]})

parameter_df.to_csv(args.ouput_path+"//infscite_parameters.csv",index=False,header=True)


gene_list_name=[i.replace(":", "_") for i in sample.dna.ids()]
gene_list=pd.DataFrame(gene_list_name)
gene_list.to_csv(args.ouput_path+"/gene_list.geneNames",index=False,header=False)