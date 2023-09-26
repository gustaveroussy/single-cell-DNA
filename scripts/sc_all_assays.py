import os
import sys
import missionbio.mosaic as ms
import missionbio.mosaic.io as mio
from yaml.loader import SafeLoader
import yaml
import os
from pathlib import Path
import re
import ast
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity
from factor_analyzer.factor_analyzer import FactorAnalyzer
from sklearn.cluster import OPTICS
from sklearn.impute import KNNImputer
import argparse

parser = argparse.ArgumentParser(description='Filtering low variants quality & plotting diagnostic plots')

parser.add_argument("--input_h5", help="input h5 path file",required=True,action='store')
parser.add_argument("--output_h5", help="output h5 path file",required=True,action='store')
parser.add_argument("--sample_name", help="wilcard sample",required=True,action='store')
parser.add_argument("--config_file", help="config file",required=True,action='store')

args = parser.parse_args()

#os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

lib_to_import=sys.argv[0].replace("/scripts/sc_all_assays.py", "/")
sys.path.append(lib_to_import+"common")
from utils import *

config=args.config_file
config=ast.literal_eval(config)

try:
    max_snv_dim=config["snv_norm_dimred"]["max_dims"]
except KeyError:
    max_snv_dim=6


if config["type_analysis"] == "dna+protein":
    sample_dna=ms.load(config["output_sample_path"]+"/"+args.sample_name+"/h5/dna/dimred_QC_max_n_dim_"+str(max_snv_dim)+".h5", raw=False, apply_filter=True, single=True)
    sample_prot=ms.load(args.input_h5,raw=False, apply_filter=True, single=True)
    sample=ms.load(config["output_sample_path"]+"/"+args.sample_name+"/"+args.sample_name+".dna+protein.h5", raw=False, apply_filter=True, single=True)

else:
    sample_dna=ms.load(args.input_5, raw=False, apply_filter=True, single=True)
    sample=ms.load(config["output_sample_path"]+"/"+args.sample_name+"/results/"+args.sample_name+".dna.h5", raw=False, apply_filter=True, single=True)
    
try:
    snv_method_dimred=config["all_norm_dimred"]["snv"]["method_dimred"]
except KeyError:
   snv_method_dimred="pca" 
try:
    snv_dims_components=config["all_norm_dimred"]["snv"]["dims"]
except KeyError:
    snv_dims_components=6
try:
    snv_method_clustering=config["all_norm_dimred"]["snv"]["clustering_method"]
except KeyError:
    snv_method_clustering="dbscan"
try:
    snv_method_clustering_res=config["all_norm_dimred"]["snv"]["res"]
except KeyError:
    snv_method_clustering_res=0.5
try:
    cnv_dims_pca_components=config["all_norm_dimred"]["cnv"]["dims"]
except KeyError:
    cnv_dims_pca_components=6
try:
    cnv_method_clustering=config["all_norm_dimred"]["cnv"]["clustering_method"]
except KeyError:
    cnv_method_clustering="dbscan"
try:
    cnv_method_clustering_res=config["all_norm_dimred"]["cnv"]["res"]
except KeyError:
    cnv_method_clustering_res=0.5
try:
    prot_norm=config["all_norm_dimred"]["prot"]["normalization"]
except KeyError:
    prot_norm="CLR"
try:
    prot_dims_pca_components=config["all_norm_dimred"]["prot"]["dims"]
except KeyError:
    prot_dims_pca_components=6
try:
    prot_method_clustering=config["all_norm_dimred"]["prot"]["clustering_method"]
except KeyError:
    prot_method_clustering="dbscan"
try:
    prot_method_clustering_res=config["all_norm_dimred"]["prot"]["res"]
except KeyError:
    prot_method_clustering_res=0.5
try:
    bool_predict_missing_value=config["all_norm_dimred"]["snv"]["predict_missing_value"]
except KeyError:
    bool_predict_missing_value=False

print("SNV")

sample.dna=sample.dna[sample_dna.dna.barcodes(), sample_dna.dna.ids()]

if bool_predict_missing_value == True:
    X=sample.dna.layers["AF_MISSING"]
    X_nan = np.where(X==-50, np.nan, X)
    imputer = KNNImputer(n_neighbors=5)
    imput_matrix=imputer.fit_transform(X_nan)
    sample.dna.layers["AF_MISSING"]=imput_matrix

if snv_method_dimred == "pca":
    sample.dna.run_pca(components=snv_dims_components, 
                    attribute="AF_MISSING",
                    show_plot=False,
                    output_label="pca")
                        
if snv_method_dimred == "fa":
    fa = FactorAnalyzer(n_factors=snv_dims_components,rotation="varimax",method="minres")
    sample.dna.row_attrs["fa"]=fa.fit_transform(pd.DataFrame(sample.dna.layers["AF_MISSING"]))

sample.dna.run_umap(attribute=snv_method_dimred,
                        min_dist=0.0, 
                        n_neighbors=50,
                        random_state=52,
                        output_label="umap")
    
if snv_method_clustering in ["dbscan","hdbscan","graph-community"]:
    if snv_method_clustering == "dbscan":
        snv_res_clustering_method={'eps':snv_method_clustering_res}
    if snv_method_clustering == "hdbscan":
        snv_res_clustering_method={'min_cluster_size':snv_method_clustering_res}
    if snv_method_clustering == "graph-community":
        snv_res_clustering_method={'k':snv_method_clustering_res}
    
    sample.dna.cluster(attribute="umap",method=snv_method_clustering,**snv_res_clustering_method)
                            
    if snv_method_clustering == "OPTICS":
        model = OPTICS(eps=0.5, min_samples=snv_method_clustering_res)
        X=pd.DataFrame(sample.dna.row_attrs["umap"])
        y_label = model.fit_predict(X)
        y_label=np.array([str(i) for i in y_label])
        y_label = np.array([len(np.unique(y_label)) if item == '-1' else item for item in y_label])
        sample.dna.row_attrs["label"]=y_label
        
    if snv_method_clustering == "leiden" and snv_method_dimred == "pca":
        
        clustering_leiden(assay=sample,attribute=snv_method_clustering_res,resolution=snv_method_clustering_res,k=200,show_plot=False)
        
    if snv_method_clustering == "leiden" and snv_method_dimred == "fa": 
        
        print("Leiden is only intented to work with PCA as reduction dimension")
                        
print("CNV")
reads = sample.cnv.get_attribute('read_counts', constraint='row+col')
working_amplicons = (reads.median() > 0).values

sample.cnv = sample.cnv[:, working_amplicons]


cnv_select_barcodes_and_normalization(assay=sample,
                                                barcodes_list=sample.dna.barcodes())
    
    
sample.cnv.run_pca(components=cnv_dims_pca_components, 
        attribute='normalized_counts',
        show_plot=False)


sample.cnv.run_umap(attribute="pca",
                        min_dist=0.0, 
                        n_neighbors=50,
                        random_state=52,
                        output_label="umap")

if cnv_method_clustering == "dbscan":
    cnv_res_clustering_method={'eps':cnv_method_clustering_res}
if cnv_method_clustering == "hdbscan":
    cnv_res_clustering_method={'min_cluster_size':cnv_method_clustering_res}
if cnv_method_clustering == "graph-community":
    cnv_res_clustering_method={'k':cnv_method_clustering_res}
    
sample.cnv.cluster(attribute="pca",
                    method=cnv_method_clustering,
                    **cnv_res_clustering_method)
                        
print("PROTEIN")
if config["type_analysis"] == "dna+protein":
    #sample_prot.protein.col_attrs['id']= np.array([i_string.replace('.','-') for i_string in sample_prot.protein.col_attrs['id']])
    sample.protein=sample.protein[sample_prot.protein.barcodes(),sample_prot.protein.ids()]

    if prot_norm in ["NSP", "CLR","asinh"]:
        sample.protein.normalize_reads(prot_norm)
    if prot_norm == "DSB":
        sample.protein.layers["normalized_counts"]=sample_prot.protein.layers["normalized_counts_DSB"]


    sample.protein.run_pca(components=prot_dims_pca_components, 
                attribute="normalized_counts",
                show_plot=False,
                output_label="pca")

    sample.protein.run_umap(attribute="pca",
                        min_dist=0.0, 
                        n_neighbors=50,
                        random_state=52,
                        output_label="umap")

    if prot_method_clustering == "dbscan":
        prot_res_clustering_method={'eps':prot_method_clustering_res}
    if prot_method_clustering == "hdbscan":
        prot_res_clustering_method={'min_cluster_size':prot_method_clustering_res}
    if prot_method_clustering == "graph-community":
        prot_res_clustering_method={'k':prot_method_clustering_res}

    sample.protein.cluster(attribute="umap",
                            method=prot_method_clustering,
                            **prot_res_clustering_method)
        

print("Making figures")
annotation_df= pd.read_csv(config["output_sample_path"]+"/"+args.sample_name+"/dna/annotation/QC_annotation.csv",sep=",")

annotation_df["Variant ID"].values

if config["type_analysis"] == "dna+protein":
    
    intersect_barcodes = np.intersect1d(sample.dna.barcodes(), sample.protein.barcodes())
    sample.protein = sample.protein[intersect_barcodes,sample.protein.ids()]
    sample.dna= sample.dna[intersect_barcodes,annotation_df["Variant ID"].values]
    
    intersect_barcodes2 = np.intersect1d(sample.dna.barcodes(), sample.cnv.barcodes())
    sample.cnv = sample.cnv[intersect_barcodes2,sample.cnv.ids()]
    sample.dna= sample.dna[intersect_barcodes2,annotation_df["Variant ID"].values]

if config["type_analysis"] == "dna":
    
    intersect_barcodes = np.intersect1d(sample.dna.barcodes(), sample.cnv.barcodes())
    sample.cnv = sample.cnv[intersect_barcodes,sample.cnv.ids()]
    sample.dna= sample.dna[intersect_barcodes,annotation_df["Variant ID"].values]
    
mio.save(sample=sample,path=args.output_h5,raw=False)
    
annotation = sample.dna.get_annotations()

#Store the annotations in the dna assay as a new column attribute
for col, content in annotation.items():
    sample.dna.add_col_attr(col, content.values)
        
sample.dna.set_ids_from_cols(["Gene", "CHROM", "POS", "REF", "ALT"])
    
str_path_result_path=config["output_sample_path"]+"/"+args.sample_name+"/all/"
directory_result=Path(str_path_result_path).mkdir(parents=True, exist_ok=True)
    
pval_cnv, tstat_cnv = sample.cnv.test_signature("normalized_counts")
pval_cnv
pval_cnv = pval_cnv + 10 ** -50
pvals_cnv = -np.log10(pval_cnv) * (tstat_cnv > 0)
    
if len(np.unique(sample.cnv.row_attrs["label"])) > 1:
    pval_index_list=list(enumerate(pval_cnv.mean(axis=0).values))
        
    amplicons_to_keep=np.array([sample.cnv.ids()[amp] for amp in [pval_index_list[i][0] for i in range(len(pval_index_list)) if pval_index_list[i][1] < 0.05]])
        
    sample.cnv = sample.cnv[:,amplicons_to_keep]
    
    
fig=sample.heatmap(clusterby='dna', sortby='dna', flatten=False)
fig.data[2].colorscale = 'magma'
fig.write_html(str_path_result_path+"combined_sortby_dna.html")

fig = sample.heatmap(clusterby='cnv', sortby='cnv', flatten=False)
fig.data[2].colorscale = 'magma'
fig.write_html(str_path_result_path+"combined_sortby_cnv.html")
    
if config["type_analysis"] == "dna+protein":    
    fig = sample.heatmap(clusterby='protein', sortby='protein', flatten=False)
    fig.data[2].colorscale = 'magma'
    fig.write_html(str_path_result_path+"combined_sortby_protein.html")

    protein_umap = sample.protein.row_attrs["umap"]
    feats = sample.dna.ids()[::1] # plot the first 3 variants
    fig=sample.dna.scatterplot(attribute=protein_umap, colorby="AF", features=feats)
    fig.write_html(str_path_result_path+"UMAP_prot_colorby_allelic_fraction.html")
    
    snv_umap = sample.dna.row_attrs["umap"]
    feats = sample.protein.ids()[::1] # plot the first 3 variants
    fig=sample.protein.scatterplot(attribute=snv_umap, colorby="normalized_counts", features=feats)
    fig.write_html(str_path_result_path+"UMAP_snv_colorby_protein_expression.html")
    
    sample.protein.row_attrs["dna_label"]=sample.dna.row_attrs["label"]
    fig=sample.protein.ridgeplot(attribute='normalized_counts',features=sample.protein.ids(),splitby="dna_label")
    fig.write_html(str_path_result_path+"ridgeplot_splitby_snv_clustering.html")
    
    fig=sample.protein.scatterplot(attribute="umap",colorby="label")
    fig.write_html(str_path_result_path+"UMAP_prot_colorby_label.html")

# Then, plot the umap colored by NGT_FILTERED for certain variants
cnv_umap = sample.cnv.row_attrs["umap"]
feats = sample.dna.ids()[::1] # plot the first 3 variants
fig=sample.dna.scatterplot(attribute=cnv_umap, colorby="NGT", features=feats)
fig.write_html(str_path_result_path+"UMAP_cnv_colorby_genotype.html")

feats = sample.dna.ids()[::1] # plot the first 3 variants
fig=sample.dna.scatterplot(attribute="umap", colorby="AF_MISSING", features=feats)
fig.write_html(str_path_result_path+"UMAP_snv_colorby_allelic_fraction.html")
    
fig=sample.dna.scatterplot(attribute="umap", colorby="NGT", features=feats)
fig.write_html(str_path_result_path+"UMAP_snv_colorby_genotype.html")

#sample.dna.layers[""] # plot the first 3 variants
#fig = sample.dna.scatterplot(attribute="umap", colorby="AF_MISSING", features=feats)
#fig.write_html(str_path_result_path+"UMAP_snv_colorby_allelic_fraction.html")

df_dna=pd.DataFrame(sample.dna.layers["AF_MISSING"],columns=sample.dna.col_attrs["id"],index=sample.dna.row_attrs["barcode"])
df_dna["label"]=sample.dna.row_attrs["label"]
    
fig,axes=plt.subplots(len(sample.dna.col_attrs["id"]),1,figsize=(1*6,len(sample.dna.col_attrs["id"])*6))

for i in range(len(sample.dna.col_attrs["id"])):
    sns.boxplot(ax=axes[i],data=df_dna, x="label", y=sample.dna.col_attrs["id"][i])
    axes[i].set_title (sample.dna.col_attrs["id"][i],fontdict={'fontsize': 10})
    
fig.savefig(str_path_result_path+"boxplot_allelic_fraction_clusterby_label.pdf")
plt.close(fig)
        
fig=sample.dna.scatterplot(attribute="umap",colorby="label")
fig.write_html(str_path_result_path+"UMAP_snv_colorby_label.html")
    
fig=sample.cnv.scatterplot(attribute="umap",colorby="label")
fig.write_html(str_path_result_path+"UMAP_cnv_colorby_label.html")
  
fig = sample.heatmap(clusterby='dna', sortby='dna', flatten=True,quantify_dna_mut=True,flatten_samples=True)
fig.data[2].colorscale = 'magma'

#fig.layout.xaxis3.ticktext = sample.cnv.col_attrs['gene_name'].copy()
fig.write_html(str_path_result_path+"combined_sortby_dna_flatten.html")

try:
    list_variants_of_interest=config["all_norm_dimred"]["variants_of_interest"]
except KeyError:
    list_variants_of_interest=None
try:
    chr_of_interest=config["all_norm_dimred"]["chr_of_interest"]
except KeyError:
    chr_of_interest=None

if list_variants_of_interest != None:
    make_label_for_each_variants(sample,list_variants_of_interest)
    
    cartesian_product_of_variants_of_interest=make_cartesian_product_of_variants(sample)
    
    sample.dna.row_attrs["label_genotype"]=make_label(sample,cartesian_product_of_variants_of_interest)
    
    new_palette={np.unique(sample.dna.row_attrs["label_genotype"])[i]:ms.COLORS[i] for i in range(0,len(np.unique(sample.dna.row_attrs["label_genotype"])))}
    
    sample.dna.set_palette(new_palette)
    
    sample.cnv.set_palette(new_palette)
    
    sample.cnv.compute_ploidy(diploid_cells=sample.dna.barcodes())
        
    fig=sample.dna.scatterplot(attribute="umap",colorby="label_genotype")
    fig.write_html(str_path_result_path+"UMAP_snv_colorby_voi.html")
    
    fig=sample.dna.heatmap(attribute="AF_MISSING",splitby="label_genotype")
    fig.write_html(str_path_result_path+"heatmap_snv_colorby_label_genotype.html")
    
    sample.cnv.row_attrs["label_genotype"]=sample.dna.row_attrs["label_genotype"]
    
    fig=sample.cnv.heatmap('ploidy',features=chr_of_interest,convolve=100,splitby="label_genotype")
    fig.write_html(str_path_result_path+"heatmap_ploidy_cnv_colorby_voi_convolve.html")
    
    fig=sample.cnv.heatmap('ploidy', features=chr_of_interest,splitby="label_genotype")
    fig.write_html(str_path_result_path+"heatmap_ploidy_cnv_colorby_voi.html")
    
    str_path_result_path_genotype_label=config["output_sample_path"]+"/"+args.sample_name+"/all/label_genotype/"
    directory_result=Path(str_path_result_path).mkdir(parents=True, exist_ok=True)
    
    #fig=sample.heatmap(clusterby='label_genotype', sortby='dna', flatten=False)
    #fig.data[2].colorscale = 'magma'
    #fig.write_html(str_path_result_path+"combined_sortby_label_genotype.html")

#for i in regex_search_pattern_in_list(r"^\w*_label$",sample.dna.row_attrs.keys()):
#    fig=sample.dna.heatmap(attribute="AF_MISSING",splitby=i)
#    fig.write_html(str_path_result_path_genotype_label+"heatmap_snv_splitby_"+i+"_allelic_fraction.html")
    
#    fig=sample.dna.scatterplot(attribute="umap",colorby=i)
#    fig.write_html(str_path_result_path_genotype_label+"UMAP_snv_colorby_"+i+".html")
    

#for i in np.unique(sample.cnv.row_attrs["label_genotype"]):
#    fig=sample.cnv.plot_ploidy(i)
#    fig.write_html(str_path_result_path_genotype_label+"cnv_ploidy_plot_"+i+".html")

    if config["type_analysis"] == "dna+protein":
        sample.protein.row_attrs["label_genotype"]=sample.dna.row_attrs["label_genotype"]
        sample.protein.set_palette(new_palette)
        table_distribution = {'clusters':sample.protein.row_attrs["label"],'annotation':sample.protein.row_attrs["label_genotype"]}
        table_distribution_df = pd.DataFrame(table_distribution,index=sample.protein.barcodes())
        table_distribution_df=pd.crosstab(index=table_distribution_df['annotation'], columns=table_distribution_df['clusters'])
        
        table_distribution_df=pd.DataFrame(table_distribution_df.T.values,columns=np.array(table_distribution_df.T.columns),index=np.array(table_distribution_df.T.index))
        
        new=pd.DataFrame()
        columns_to_compute=table_distribution_df.columns
        for i in columns_to_compute:
            new[i] = (table_distribution_df[i] / table_distribution_df[i].sum()) * 100
        
        
        new.T.plot(kind='bar',stacked=True,color=sample.protein.get_palette()).legend(loc='upper left')
        plt.legend(loc=(1.05, 0.5))
        plt.gcf().set_size_inches(10, 10)
        plt.savefig(str_path_result_path+'barplot_merged_protein_distribution.png')
        
        fig=sample.protein.heatmap(attribute='normalized_counts',splitby="label_genotype")
        fig.write_html(str_path_result_path+"heatmap_proteine_splitby_genotype_label.html")
