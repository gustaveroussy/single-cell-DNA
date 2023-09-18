import sys
import missionbio.mosaic as ms
import missionbio.mosaic.io as mio
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import ast
from pathlib import Path
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity
from factor_analyzer.factor_analyzer import FactorAnalyzer
from sklearn.impute import KNNImputer

import argparse

parser = argparse.ArgumentParser(description='Filtering low variants quality & plotting diagnostic plots')

parser.add_argument("--input_h5", help="input h5 path file",required=True,action='store')
parser.add_argument("--ouput_h5", help="output h5 path file",required=True,action='store')
parser.add_argument("--sample_name", help="sample name",required=True,action='store')
parser.add_argument("--config_file", help="sample name",required=True,action='store')

args = parser.parse_args()

lib_to_import=sys.argv[0].replace("/scripts/sc_dna_preprocessing.py", "/")
sys.path.append(lib_to_import+"/common")
from utils import *

sample = ms.load(args.input_h5, raw=False, apply_filter=True, single=True)

config=args.config_file
config=ast.literal_eval(config)

main_output_path=config["output_sample_path"]+"/"+args.sample_name+"/"

Path(main_output_path+"/dna/annotation").mkdir(parents=True, exist_ok=True)

try:
    annotation = sample.dna.get_annotations()

    for col, content in annotation.items():
        sample.dna.add_col_attr(col, content.values)

        annotation.to_csv(main_output_path+"dna/annotation/all_variant_annotation.csv",index=False)
        
    advanced_filtration=True
except NameError:
    print("MB annotation API is not available")
    advanced_filtration=False
try:
    min_dp=config["filtering"]["filtering_variants"]["min_dp"]
except KeyError:
    min_dp=10
try:
    min_gq=config["filtering"]["filtering_variants"]["min_gq"]
except KeyError:
    min_gq=30
try:
    vaf_ref=config["filtering"]["filtering_variants"]["vaf_ref"]
except KeyError:
    vaf_ref=5
try:
    vaf_hom=config["filtering"]["filtering_variants"]["vaf_hom"]
except KeyError:
    vaf_hom=95
try:
    vaf_het=config["filtering"]["filtering_variants"]["vaf_het"]
except KeyError:
    vaf_het=35
try:
    min_prct_cells=config["filtering"]["filtering_variants"]["min_prct_cells"]
except KeyError:
    min_prct_cells=50
try:
    min_mut_prct_cells=config["filtering"]["filtering_variants"]["min_mut_prct_cells"]
except KeyError:
    min_mut_prct_cells=1
try:
    max_vaf_percent=config["filtering"]["max_vaf_percent"]
except KeyError:
    max_vaf_percent=95
try:
    na_filtering_percent=config["filtering"]["filter_na_percent"]
except KeyError:
    na_filtering_percent=25
try:
    bool_predict_missing_value=config["filtering"]["predict_missing_value"]
except KeyError:
    bool_predict_missing_value=False
try:
    filtering_na=config["filtering"]["filter_na"]
except KeyError:
    filtering_na=False

print('\nFiltering low quality variants')
dna_vars = sample.dna.filter_variants(min_dp=min_dp,
min_gq=min_gq,
vaf_ref=vaf_ref,
vaf_hom=vaf_hom,
vaf_het=vaf_het,
min_prct_cells=min_prct_cells,
min_mut_prct_cells=min_mut_prct_cells)

try:
    annotation = sample.dna.get_annotations()

    for col, content in annotation.items():
        sample.dna.add_col_attr(col, content.values)
        
except NameError:
    print("MB annotation API is not available")
    annotation=pd.DataFrame()


# if needed shrunk variants name that are too long and create axis issue
sample=change_too_long_variant_id(sample)

try:
    target_variants=config["filtering"]["whitelist"]
except KeyError:
    target_variants=[]

# Combine whitelisted and filtered variants
# It keep variants that would have been remove or consider as false-positive
final_vars = list(set(list(dna_vars) + target_variants))
sample.dna = sample.dna[sample.dna.barcodes(), final_vars]
# First diagnostic plot to evaluate variant quality

print("\n Saving diagnostic figures")
Path(main_output_path+"dna/diagnostic").mkdir(parents=True, exist_ok=True)

fig_pyplot = sample.dna.stripplot(attribute='AF_MISSING', colorby='GQ')
fig_pyplot.write_html(main_output_path+"dna/diagnostic/stripplot.html")

# Second diagnostic plot: use this heatmap to (de)-select variants with little variance
fig_pyplot = sample.dna.heatmap(attribute='NGT_FILTERED')
fig_pyplot.write_html(main_output_path+"dna/diagnostic/heatmap.NGT_FILTERED.html")

# Third diagnostic plot: use this heatmap to (de) -select variants with low genotype quality (from GATK)

fig_pyplot = sample.dna.heatmap(attribute='GQ')
fig_pyplot.write_html(main_output_path+"dna/diagnostic/heatmap_GQ.html")


print("\nAnnotation & filtering artefactual variants")
try:
    annotation = sample.dna.get_annotations()

    for col, content in annotation.items():
        sample.dna.add_col_attr(col, content.values)
        
except NameError:
    print("\nMB annotation API is not available")
    annotation=pd.DataFrame()
    
genotyped_cells_compute(sample)
mutated_cells_compute(sample)
vaf_by_read_count(sample)
vaf_by_cell_count(sample)

sample.dna.col_attrs['ratio VAF Read Count']=sample.dna.col_attrs['VAF by read count']/(sample.dna.col_attrs['VAF by cell count']-1)
sample.dna.col_attrs['ratio VAF Cell Count']=sample.dna.col_attrs['VAF by cell count']/(sample.dna.col_attrs['VAF by read count']-1)


annotation['Genotyped Cells (%)']=sample.dna.col_attrs['Genotyped Cells']
annotation['Mutated Cells (%)']=sample.dna.col_attrs['Mutated Cells']
annotation['VAF by cell count (%)']=sample.dna.col_attrs['VAF by cell count']
annotation['VAF by read count (%)']=sample.dna.col_attrs['VAF by read count']
annotation['ratio VAF Read Count'] = sample.dna.col_attrs['ratio VAF Read Count']
annotation['ratio VAF Cell Count'] = sample.dna.col_attrs['ratio VAF Cell Count']

annotation = annotation.loc[(annotation['ratio VAF Cell Count'] <= 1.5) & (annotation['ratio VAF Cell Count'] >= -1.5)]

variant_to_keep=list(annotation.index)

final_vars = list(set(list(variant_to_keep) + target_variants))
sample.dna = sample.dna[sample.dna.barcodes(), final_vars]


if config["filtering"]["filter_na"] == True:
    print("Filtering NA")
    df_na=pd.DataFrame(sample.dna.layers["NGT_FILTERED"],columns=sample.dna.ids())
    
    df_na=df_na.replace(3,np.nan)
    
    df_na=pd.DataFrame((df_na.isna().sum()/len(df_na.index))*100,columns=["percent_na"])
    
    variant_to_keep=list(df_na[df_na.percent_na <= int(na_filtering_percent)].index)
    
    final_vars = list(set(list(variant_to_keep) + target_variants))
    sample.dna = sample.dna[sample.dna.barcodes(), variant_to_keep]

mean_percent_vaf=pd.DataFrame(sample.dna.layers["AF_MISSING"],columns=sample.dna.ids()).mean(axis=0,skipna=True).values
df_vaf=pd.DataFrame(mean_percent_vaf,index=sample.dna.ids(),columns=["mean_vaf"])
variant_to_keep=list(df_vaf[df_vaf.mean_vaf <= int(max_vaf_percent)].index)


if bool_predict_missing_value == True:
    print("\nImputing missing VAF & filtering germinal variants")
    X=sample.dna.layers["AF_MISSING"]
    X_nan = np.where(X==-50, np.nan, X)
    imputer = KNNImputer(n_neighbors=5)
    imput_matrix=imputer.fit_transform(X_nan)
    sample.dna.layers["AF_MISSING"]=imput_matrix

final_vars = list(set(list(variant_to_keep) + target_variants))
sample.dna = sample.dna[sample.dna.barcodes(), final_vars]

print("\nSaving Annotation")
annotation.to_csv(main_output_path+"dna/annotation/QC_annotation.csv")


data_for_ann={'Mutated cells':sample.dna.col_attrs['Mutated Cells'],
              'Genotyped cells':sample.dna.col_attrs['Genotyped Cells']}


df=pd.DataFrame(data_for_ann,index=sample.dna.col_attrs['id'])

fig,ax= plt.subplots()
fig.set_size_inches(1*len(sample.dna.ids()), 20)
plt.xticks(rotation=45)
plt.title("Variants per Mutated cells (%)")
sns.barplot(x='index', y='value', hue='percent',
            data=df.reset_index().melt(id_vars='index', var_name='percent'),ax=ax)
fig.savefig(main_output_path+"dna/diagnostic/barplot_mutated.png")
plt.close(fig)

if advanced_filtration == True:

    annotation['DANN'] = np.where(annotation['DANN'] == '', 0.701, annotation['DANN']) 
    # certains variants n'ont pas de score DANN 
    # donc j'assigne une valeur de 0.701 afin qu'il ne soit pas éliminer
    annotation = annotation.loc[(annotation['DANN'] >= 0.7)]
    # on garde les variants dont la valeur de score de DANN est supérieur à 0.7
    #annotation = annotation.loc[(annotation['ratio VAF Cell Count'] <= 1.5) & (annotation['ratio VAF Cell Count'] >= -1.5)]
    #on garde les variants dont la valeur de ratio VAF Cell Count est inférieur à 1.5 ou supérieur à -1.5 
    # ces variants sont considéré comme des potentiels artefacts
    annotation = annotation.loc[(annotation['ClinVar'] != "Benign") & (annotation['ClinVar'] != "Likely Benign")]
    # on garde les variants qui ne sont pas bénin ou probablement bénin
    annotation = annotation.loc[(annotation['Function'] != 'intronic') & (annotation['Coding impact'] != 'synonymous')]
    # on garde les variants dont la mutation n'est pas intronique ou dont l'impact de la mutation n'est pas synonyme
    annotation = annotation.loc[(annotation['Coding impact']) != '']

    annotation.to_csv(main_output_path+"dna/annotation/QC_advanced_annotation.csv")

mio.save(sample=sample,path=args.ouput_h5,raw=False)



chi_square_value,p_value=calculate_bartlett_sphericity(pd.DataFrame(sample.dna.layers["AF_MISSING"]))
chi_square_value, p_value

Path(main_output_path+"dna/elbow").mkdir(parents=True, exist_ok=True)
print('\nElbow plot FA')
if p_value < 0.05:
    print("\nFactor Analysis : The test Bartlett’s test  was statistically significant, indicating that the observed correlation matrix is not an identity matrix.")
    
    fa = FactorAnalyzer(n_factors=len(sample.dna.col_attrs["id"]),rotation="varimax",method="minres")
    fa.fit(pd.DataFrame(sample.dna.layers["AF_MISSING"]))
        
    
    ev, v = fa.get_eigenvalues()
    
    plt.scatter(range(1,pd.DataFrame(sample.dna.layers["AF_MISSING"]).shape[1]+1),ev)
    plt.plot(range(1,pd.DataFrame(sample.dna.layers["AF_MISSING"]).shape[1]+1),ev)
    plt.title('Scree Plot')
    plt.xlabel('Factors')
    plt.ylabel('Eigenvalue')
    plt.grid()
    plt.savefig(main_output_path+"dna/elbow/elbowplot_fa.png")

print('\nElbow plot PCA\n')
variants_number = sample.dna.shape[1]

sample.dna.run_pca(components=variants_number, 
                attribute="AF_MISSING",
                show_plot=True,
                output_label="pca")

plt.savefig(main_output_path+"dna/elbow/elbowplot_pca.png")
plt.close()