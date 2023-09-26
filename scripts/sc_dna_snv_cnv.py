import sys
import missionbio.mosaic as ms
import missionbio.mosaic.io as mio
import pandas as pd
import re
import ast
from factor_analyzer.factor_analyzer import FactorAnalyzer
from sklearn.impute import KNNImputer
import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description='Filtering low variants quality & plotting diagnostic plots')

parser.add_argument("--input_h5", help="input h5 path file",required=True,action='store')
parser.add_argument("--output_h5", help="output h5 path file",required=True,action='store')
parser.add_argument("--sample_name", help="wilcard sample",required=True,action='store')
parser.add_argument("--config_file", help="config file",required=True,action='store')

args = parser.parse_args()

lib_to_import=sys.argv[0].replace("/scripts/sc_dna_snv_cnv.py", "/")
sys.path.append(lib_to_import+"common")
from utils import *

config=args.config_file
config=ast.literal_eval(config)

try:
    method=config["snv_norm_dimred"]["method_dimred"]
except KeyError:
    method="pca"
try:
    snv_component_range=config["snv_norm_dimred"]["max_dims"]
except KeyError:
    snv_component_range=6
try:
    clustering_method_list=config["snv_norm_dimred"]["clustering_method"]
except KeyError:
    clustering_method_list=["dbscan"]
try:
    bool_predict_missing_value=config["filtering"]["predict_missing_value"]
except KeyError:
    bool_predict_missing_value=False
try:
    cnv_component_range=config["cnv_norm_dimred"]["max_dims"]
except KeyError:
    cnv_component_range=6
try:
    cnv_clustering_method_list=config["cnv_norm_dimred"]["clustering_method"]
except KeyError:
    cnv_clustering_method_list=["dbscan"]

        
sample_path = args.input_h5
sample = ms.load(sample_path, raw=False, apply_filter=False, single=True)

print("\nSTART SNV analysis")
if bool_predict_missing_value == True:
# -- predicting missing value --
    X=sample.dna.layers["AF_MISSING"]
    # replace -50 value with NaN
    X_nan = np.where(X==-50, np.nan, X)
    # K-nearest neighbors imputer
    imputer = KNNImputer(n_neighbors=20)
    imput_matrix=imputer.fit_transform(X_nan)
    #replace with new predict matrix
    sample.dna.layers["AF_MISSING"]=imput_matrix

print("\nReduction Dimension")
# select reduction dimension
if method == "pca":
    for i_component in range(2,snv_component_range+1):
            sample.dna.run_pca(components=i_component, 
                    attribute="AF_MISSING",
                    show_plot=False,
                    output_label="pca_ndim_"+str(i_component))
if method == "fa":
    for i_component in range(2,snv_component_range+1):
        fa = FactorAnalyzer(n_factors=i_component,rotation="varimax",method="minres")
        sample.dna.row_attrs["fa_ndim_"+str(i_component)]=fa.fit_transform(pd.DataFrame(sample.dna.layers["AF_MISSING"]))

print("\nUMAP")
# projection of each pca components into umap
sample=make_umap(sample,'dna')
print("\nClustering in progress")
# clustering using different methods
for i_method in clustering_method_list:
    sample=make_clustering(assay=sample,
                            arg_attribute_assay="dna",
                            arg_method_clustering=i_method,
                            arg_max_components=snv_component_range)
print("\nClustering done")
        
print("\nSaving figures")
for i_method in clustering_method_list:
    str_path_result_path=config["output_sample_path"]+"/"+args.sample_name+"/dna/clustering/"+method+"/"+i_method+"/"
    directory_result=Path(str_path_result_path).mkdir(parents=True, exist_ok=True)
    plot_clustering(assay=sample,arg_attribute_assay="dna",arg_method_clustering=i_method,arg_max_components=snv_component_range,args_directory_result=str_path_result_path,args_normalization=None)
print("\nEND SNV analysis")   

print("\nSTART CNV analysis")
reads = sample.cnv.get_attribute('read_counts', constraint='row+col')
working_amplicons = (reads.median() > 0).values
sample.cnv = sample.cnv[:, working_amplicons]
sample.cnv = sample.cnv[sample.dna.barcodes(),:]
sample.cnv.normalize_reads()
sample.cnv.compute_ploidy(diploid_cells=sample.dna.barcodes())

print("\nReduction Dimension")
# select reduction dimension
for i_component in range(2,cnv_component_range+1):
    sample.cnv.run_pca(components=i_component, 
        attribute="normalized_counts",
        show_plot=False,
        output_label="pca_ndim_"+str(i_component))

print("\nUMAP")
# projection of each pca components into umap
sample=make_umap(sample,'cnv')
print("\nClustering in progress")
for i_method in cnv_clustering_method_list:
    sample=make_clustering(assay=sample,
                            arg_attribute_assay="cnv",
                            arg_method_clustering=i_method,
                            arg_max_components=cnv_component_range)
print("\nClustering done")

for i_method in cnv_clustering_method_list:
    str_path_result_path=config["output_sample_path"]+"/"+args.sample_name+"/cnv/clustering/"+i_method+"/"
    directory_result=Path(str_path_result_path).mkdir(parents=True, exist_ok=True)
    plot_clustering(assay=sample,arg_attribute_assay="cnv",arg_method_clustering=i_method,arg_max_components=cnv_component_range,args_directory_result=str_path_result_path,args_normalization=None)


print("\nEND CNV analysis")
mio.save(sample=sample,path=args.output_h5,raw=False)