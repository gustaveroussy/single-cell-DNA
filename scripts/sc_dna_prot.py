import os
import sys
import missionbio.mosaic as ms
from missionbio.mosaic.constants import COLORS
import missionbio.mosaic.io as mio
from yaml.loader import SafeLoader
import yaml
import re
import os
import ast
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description='Protein PCA & UMAP & Clustering')

parser.add_argument("--input_h5", help="input h5 path file",required=True,action='store')
parser.add_argument("--output_h5", help="output h5 path file",required=True,action='store')
parser.add_argument("--sample_name", help="wilcard sample",required=True,action='store')
parser.add_argument("--config_file", help="config file",required=True,action='store')

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

lib_to_import=sys.argv[0].replace("/scripts/sc_dna_prot.py", "/")
sys.path.append(lib_to_import+"common")
from utils import *

args = parser.parse_args()

config=args.config_file
config=ast.literal_eval(config)

try:
    normalisation_list=config["prot_norm_dimred"]["normalization"]
except KeyError:
    normalisation_list=["CLR"]
        
try:
    prot_clustering_method_list=config["prot_norm_dimred"]["clustering_method"]
except KeyError:
    prot_clustering_method_list=["dbscan"]


sample = ms.load(args.input_h5, raw=False, apply_filter=True, single=True)

max_components=len(sample.protein.col_attrs["id"])-1


for i_method_norm in normalisation_list:
        for i_pca_components in range(2,max_components):
                sample.protein.run_pca(components=i_pca_components, 
                        attribute='normalized_counts_'+str(i_method_norm),
                        show_plot=False,
                        output_label="pca_norm_"+str(i_method_norm)+"_ndim_"+str(i_pca_components))

sample=make_umap(sample,'protein')


for i_method in prot_clustering_method_list:
    sample=make_clustering(assay=sample,arg_attribute_assay="protein",arg_method_clustering=i_method,arg_max_components=max_components)
                            
                     
for i_method in normalisation_list:
    for i_method_clustering in prot_clustering_method_list:
        str_path_result_path=config["output_sample_path"]+"/"+args.sample_name+"/prot/clustering/"+i_method+"/"+i_method_clustering+"/"
        directory_result=Path(str_path_result_path).mkdir(parents=True, exist_ok=True)
        
        plot_clustering(assay=sample,
        arg_attribute_assay="protein",
        arg_method_clustering=i_method_clustering,
        arg_max_components=max_components,
        args_directory_result=str_path_result_path,
        args_normalization=i_method)
        
        fig=sample.protein.ridgeplot(attribute='normalized_counts_'+i_method,features=sample.protein.ids())
        
        Path(config["output_sample_path"]+"/"+args.sample_name+"/prot/ridgleplot/").mkdir(parents=True, exist_ok=True)
        fig.write_html(config["output_sample_path"]+"/"+args.sample_name+"/prot/ridgleplot/ridgleplot_+"+i_method_norm+".html")

mio.save(sample=sample,path=args.output_h5,raw=False)