import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import OPTICS
import math
import seaborn as sns
import pandas as pd
import os
import yaml
import math
from sklearn.cluster import OPTICS
import itertools

def flatten(l):
    """return in one list a list of list"""
    return [item for sublist in l for item in sublist]

def regex_search_pattern_in_list(pattern,list_to_search):
    """
    allow regular expression re.search() in a list and return the list of results
    """
    regex_obj_list=[re.search(pattern,elem_string) for elem_string in list_to_search]
    regex_obj_list_result=[string_list_element.group() for  string_list_element in regex_obj_list if string_list_element is not None]
    return(regex_obj_list_result)
    
def change_too_long_variant_id(assay):
    """shrunk variants name that are too long and replace the col_attrs["id"] by new names list 
    with shrunken id olders id are in conserve in "orig_id" list in the col_attrs layer"""
    
    old_list = assay.dna.col_attrs["id"]
    tmp_list = []
    for element_txt in assay.dna.col_attrs["id"]:
        if len(element_txt) > 40:
            get_end = re.split("/",element_txt)
            get_chr = re.split(":",element_txt)
            new_element = get_chr[0]+":"+get_chr[1]+":"+get_chr[2][0:18]+"{...}"+"/"+get_end[1]
            tmp_list.append(new_element)
        else:
            tmp_list.append(element_txt)
    assay.dna.col_attrs["id"]=tmp_list
    assay.dna.col_attrs["orig_id"]=old_list
    return(assay)
    
def genotyped_cells_compute(assay):
    NGT_filtered_all_cell_except_missing=[[[assay.dna.layers['NGT_FILTERED'][i][k] 
                    for i in range(len(assay.dna.layers['NGT_FILTERED']))].count(j)
                   for j in range(0,3)]
                  for k in range(len(assay.dna.col_attrs['filtered']))]
    assay.dna.col_attrs['Genotyped Cells']=np.array([sum(i)/assay.dna.shape[0]*100 
                                                     for i in NGT_filtered_all_cell_except_missing])
    
    return(assay)

def mutated_cells_compute(assay):
    NGT_filtered=[[[assay.dna.layers['NGT_FILTERED'][i][k] 
                    for i in range(len(assay.dna.layers['NGT_FILTERED']))].count(j)
                   for j in range(1,3)]
                  for k in range(len(assay.dna.col_attrs['filtered']))]
    sum_NGT_filtered=[sum([[assay.dna.layers['NGT_FILTERED'][i][k] 
                    for i in range(len(assay.dna.layers['NGT_FILTERED']))].count(j)
                   for j in range(0,3)])
                  for k in range(len(assay.dna.col_attrs['filtered']))]
    
    assay.dna.col_attrs['Mutated Cells']=np.array([sum(NGT_filtered[i])/sum_NGT_filtered[i]*100 
                                                   for i in range(len(NGT_filtered))])
    
    return(assay)

def vaf_by_read_count(assay):
    
    af = assay.dna.get_attribute("AF_FILTERED", "row+col")
    #no data is coded by -50
    af = af.replace(-50, np.nan)
    VAF_read= af.mean()
    assay.dna.col_attrs['VAF by read count']=np.array(VAF_read)
    return(assay)

def vaf_by_cell_count(assay):
    #VAF by cell count code, on NGT_FILTERED
    ngt_filtered =  assay.dna.get_attribute("NGT_FILTERED", "row+col")
    #no data is coded by 3
    ngt_filtered = ngt_filtered.replace(3, np.nan)
    VAF_cell = 100 * ngt_filtered.sum() / (2 * ngt_filtered.shape[0])
    assay.dna.col_attrs['VAF by cell count']=np.array(VAF_cell)
    return(assay)

def order_annotation(assay,df,list_to_order,ascending_bool=False):
    """ 
    order variants in function of specific column (usually DANN or coding impact)
    """
    for col, content in df.items():
        assay.dna.add_col_attr(col, content.values)
        ann = df.sort_values(by=list_to_order, ascending=ascending_bool)

    return(ann)
    
def helper_umap(assay,arg_attribute_assay,dimred_attribute_range):
    """
    help make_umap to process
    """
    
    for elem_dimred in dimred_attribute_range:
        # split key reduction dimension label and get the number of components in order to create the umap output label
        if arg_attribute_assay == 'protein':
            dim_string_tmp = re.split(r'_',elem_dimred)
            output_label_umap_str='umap_'+dim_string_tmp[2]+'_ndim_'+str(dim_string_tmp[-1])
        else:
            dim_string_tmp = re.split(r'_',elem_dimred)
            output_label_umap_str='umap_ndim_'+str(dim_string_tmp[-1])
            
        # be careful mosaic UMAP is a bit different from scikit-learn & uwot
        getattr(assay,arg_attribute_assay).run_umap(attribute=elem_dimred,
                                                output_label=output_label_umap_str,
                                                random_state=42,
                                                min_dist=0.0,
                                                n_neighbors=50)

def make_umap(assay,
            arg_attribute_assay):
    """
    make umap and store each different run in a different output_label
    """
    # check sample assay type and look for pca or factor analysis reduction dimension 
    #i n the sample.assay.row_attrs dictionnary
    if arg_attribute_assay in ['dna','cnv']:
        for i in ["^pca_ndim_\d*$","^fa_ndim_\d*$"]:
            dimred_attribute_range = regex_search_pattern_in_list(i,getattr(assay,arg_attribute_assay).row_attrs.keys())
            
            helper_umap(assay,arg_attribute_assay,dimred_attribute_range)
           
    if arg_attribute_assay == 'protein':
        for i in ["^pca_norm_CLR_ndim_\d*$","^pca_norm_DSB_ndim_\d*$","^pca_norm_NSP_ndim_\d*$","^pca_norm_asinh_ndim_\d*$"]:
            dimred_attribute_range = regex_search_pattern_in_list(i, getattr(assay,arg_attribute_assay).row_attrs.keys())
            
            helper_umap(assay,arg_attribute_assay,dimred_attribute_range)
    
    return(assay)


def helper_make_cluster(assay,arg_attribute_assay,arg_method_clustering,clustering_value_range,pattern_element_umap):

    for clustering_value in clustering_value_range:
        
        if arg_method_clustering=='dbscan':

            getattr(assay,arg_attribute_assay).cluster(attribute=pattern_element_umap,method=arg_method_clustering,eps=clustering_value)
            tmp= getattr(assay,arg_attribute_assay).row_attrs['label']
            getattr(assay,arg_attribute_assay).row_attrs[pattern_element_umap+'_dbscan_'+str(round(clustering_value,3))]=tmp
                                
        if arg_method_clustering=='hdbscan':

            getattr(assay,arg_attribute_assay).cluster(attribute=pattern_element_umap,method=arg_method_clustering,min_cluster_size=clustering_value)
            tmp=getattr(assay,arg_attribute_assay).row_attrs['label']
            getattr(assay,arg_attribute_assay).row_attrs[pattern_element_umap+'_hdbscan_'+str(clustering_value)]=tmp

        if arg_method_clustering=="OPTICS":
                                
            model = OPTICS(eps=0.5, min_samples=clustering_value)
            X=pd.DataFrame(getattr(assay,arg_attribute_assay).row_attrs[pattern_element_umap])
            y_label = model.fit_predict(X)
            y_label=np.array([str(i) for i in y_label])
            y_label = np.array([len(np.unique(y_label)) if item == '-1' else item for item in y_label])
            getattr(assay,arg_attribute_assay).row_attrs[pattern_element_umap+"_OPTICS_"+str(clustering_value)]=y_label

def make_clustering(assay,arg_attribute_assay,arg_method_clustering,arg_max_components):
    # check if method of clustering is available
    if arg_method_clustering not in ['dbscan','hdbscan','graph-community','leiden']:
        print('Clustering method not found. Clustering method accepted are : \n - dbscan \n - hdbscan \n - graph-community \n Leiden')

    else:
        # clustering range for dbsan needed np.arange
        if arg_method_clustering=='dbscan':
            
            clustering_value_range=[clustering_value for clustering_value in np.arange(0.1,1.0,0.1)]
            
        else:
            
            clustering_value_range=[clustering_value for clustering_value in range(10,120,10)]

        # if method is not a graph clustering method
        if arg_method_clustering in ['dbscan','hdbscan']:
            # iterate over max number of components
            if arg_attribute_assay in ['protein']:
                for element_in_range in range(2,arg_max_components+1):
                    for pattern_method_norm in ["umap_DSB_ndim_"+str(element_in_range),
                                                        "umap_CLR_ndim_"+str(element_in_range),
                                                        "umap_asinh_ndim_"+str(element_in_range),
                                                        "umap_NSP_ndim_"+str(element_in_range)]:
                        len_pca=len(regex_search_pattern_in_list(pattern_method_norm, getattr(assay,arg_attribute_assay).row_attrs.keys()))
                        
                        if len_pca > 0:
                            pattern_element_umap=pattern_method_norm
                            helper_make_cluster(assay,arg_attribute_assay,arg_method_clustering,clustering_value_range,pattern_element_umap)
                # if assay is not protein  
            if arg_attribute_assay in ["dna","cnv"]:
                for element_in_range in range(2,arg_max_components+1):
                    pattern_element_umap="umap_ndim_"+str(element_in_range)
                    helper_make_cluster(assay,arg_attribute_assay,arg_method_clustering,clustering_value_range,pattern_element_umap)

                    
        elif arg_method_clustering in ["graph-community","leiden"]:
            if arg_method_clustering=='leiden':
                clustering_value_range=[clustering_value for clustering_value in np.arange(0.1,1.2,0.1)]
            if arg_method_clustering=='graph-community':
                clustering_value_range=[clustering_value for clustering_value in range(10,120,10)]

            for i in ["^pca_ndim_\d*$","^fa_ndim_\d*$"]:
                    dimred_list = regex_search_pattern_in_list(i,getattr(assay,arg_attribute_assay).row_attrs.keys())
            if arg_method_clustering=='graph-community':
                for i in dimred_list:
                    for clustering_value in clustering_value_range:
                        getattr(assay,arg_attribute_assay).cluster(attribute=i,method=arg_method_clustering,k=clustering_value)
                        tmp=getattr(assay,arg_attribute_assay).row_attrs['label']
                        getattr(assay,arg_attribute_assay).row_attrs[i+'_graph-community_'+str(clustering_value)]=tmp

            elif arg_method_clustering == "leiden":
                for i in dimred_list:
                    for clustering_value in clustering_value_range:
                        tmp=clustering_leiden(getattr(assay,arg_attribute_assay),clustering_value,k=200)
                        tmp=getattr(assay,arg_attribute_assay).row_attrs['label']
                        getattr(assay,arg_attribute_assay).row_attrs[i+'_leiden_'+str(clustering_value)]=tmp

    return(assay)


def clustering_leiden(assay,attribute,resolution,k,show_plot=False):
    """
    Leiden clustering
    
    modifies the id row attribute of the assay
    """
    
    print("Creating the Shared Nearest Neighbors graph.")
    knn = NearestNeighbors(n_neighbors=k, algorithm="ball_tree").fit(assay.row_attrs[attribute])
    dist, nbrs = knn.kneighbors(assay.row_attrs[attribute])
    # Set the edges
    v1 = 0
    edges = []
    for v1_nbrs in nbrs:
        for v2 in v1_nbrs:
            edges.append((v1, v2))
        v1 += 1

    edges = np.array(edges)
    # Set the weight using Jaccard similarity
    # Optimize this step further, if possible
    print("-" * 50)
    i = 0
    weights = []
    similarities = {}
    nbrs = [set(t) for t in nbrs]
    for e in edges:
        if i % (edges.shape[0] // 50) == 0:
            print("#", flush=True, end="")

        i += 1
        if (e[0], e[1]) not in similarities:
            intersection = len(nbrs[e[0]].intersection(nbrs[e[1]]))
            common = 2 * k - intersection
            similarity = intersection / common
            similarities[(e[1], e[0])] = similarity
        else:
            similarity = similarities[(e[0], e[1])]

        weights.append(similarity)
    # Construct the graphs and identify clusters for various values of k less than the given value
    print("\nIdentifying clusters using Leiden community detection.")
    edges = np.array(edges).reshape(len(nbrs), k, 2)
    weights = np.array(weights).reshape(len(nbrs), k, 1)
    sizes = np.linspace(10, k, 5) if show_plot else [k]
    clusters_found = []
    for size in sizes:
        size = int(size)
        sub_edges = np.array(edges)[:, :size, :]
        sub_edges = sub_edges.reshape(sub_edges.shape[0] * sub_edges.shape[1], sub_edges.shape[2])
        sub_weights = np.array(weights)[:, :size]
        sub_weights = sub_weights.reshape(sub_weights.shape[0] * sub_weights.shape[1])
        
        g = igraph.Graph()
        g.add_vertices(len(assay.barcodes()))
        g.add_edges(sub_edges)
        g.es["weight"] = sub_weights
        
        vc = g.community_leiden(weights="weight",resolution_parameter=resolution)
        vc_labels = np.array(vc.as_cover().membership).flatten()
        
        clusters_found.append((size, len(set(vc_labels))))

        print(f"\nNumber of clusters found: {clusters_found[-1][1]}.")
        print(f"Modularity: {vc.modularity:.3f}")

        if show_plot:
            x = np.array(clusters_found)[:, 0]
            y = np.array(clusters_found)[:, 1]
            plt.figure(figsize=(10, 10))
            plt.scatter(x, y)
            plt.ylabel("Number of clusters found")
            plt.xlabel("Number of nearest neighbors used")
        
        vc_labels=np.array([str(i) for i in vc_labels])
        new_palette={str(i):ms.COLORS[i] for i in range(0,len(np.unique(vc_labels))+1)}
        assay.set_palette(new_palette)
        assay.row_attrs["label"]=vc_labels

def plot_clustering(assay,arg_attribute_assay,arg_method_clustering,arg_max_components,args_directory_result,args_normalization=None):
    for i_pca in range(2,arg_max_components+1):
        if arg_attribute_assay in ['protein']:
            attribute_umap="umap_"+str(args_normalization)+"_ndim_"+str(i_pca)
        if arg_attribute_assay in ['dna','cnv']:
            attribute_umap="umap_ndim_"+str(i_pca)
            
        if arg_method_clustering in ["leiden","dbscan"]:
            pattern_clustering=r"^"+attribute_umap+"_"+arg_method_clustering+"_\d*.\d*$"
        else:
            pattern_clustering=r"^"+attribute_umap+"_"+arg_method_clustering+"_\d*$"
        
        clustering_label_list=regex_search_pattern_in_list(pattern_clustering,getattr(assay,arg_attribute_assay).row_attrs)
        
        nb_row=math.ceil(len(clustering_label_list)/3)
        nb_col=3
        
        tmp_dict= {}
        fig,axes = plt.subplots(nb_row,nb_col,figsize=(6*nb_col,8*nb_row))
        axes = axes.flat
        
        for i_label in range(0,len(clustering_label_list)):
            tmp_df=pd.DataFrame(getattr(assay,arg_attribute_assay).row_attrs[attribute_umap])
            tmp_df['label']=getattr(assay,arg_attribute_assay).row_attrs[clustering_label_list[i_label]]
            tmp_dict[clustering_label_list[i_label]]=tmp_df
            
            
        for ax, (i_name,i_df) in zip(axes, tmp_dict.items()):
            
            sns.scatterplot(ax=ax,data=i_df, x=0, y=1,hue='label',size=0.05,legend=False)
            sns.despine(top=True, right=True)
            ax.set_title(f'{i_name}', fontsize=8)
        
        if args_normalization == None:
            result_path=args_directory_result+arg_method_clustering+"_ndim_"+str(i_pca)+".png"
        else:   
            result_path=args_directory_result+arg_method_clustering+"_ndim_"+str(i_pca)+"_"+args_normalization+".png"
            
        fig.savefig(result_path)
        plt.close(fig)

def make_yaml(config_file):
    
    design_file=pd.read_csv(config_file["design_file"],sep="\t")
    
    
    try:
        panel_path=config_file["panel_path"]
    except KeyError:
        panel_path="/mnt/beegfs/pipelines/single-cell_dna/tapestri_database/v2/panels/Myeloid"
        
    panel_name_dna=panel_path.split("/")[::-1][0]
    
    try:
        reference_genome_path=config_file["reference_genome_path"]
    except KeyError:
        reference_genome_path='/mnt/beegfs/pipelines/single-cell_dna/tapestri_database/v2/hg19/ucsc_hg19.fa'
    try:
        reference_genome=config_file["reference_genome"]
    except KeyError:
        reference_genome='hg19'
    
    pattern_r1 = r'\w*R1\w*\.fastq\.gz'
    pattern_r2 = r'\w*R2\w*\.fastq\.gz'

    for row in range(len(design_file.index)):
        
        if config_file["type_analysis"] == "dna+protein":
            
            try:
                protein_panel=config_file["panel_protein_path"]
                pattern_fasta3=r'\S*adapters_3\S*\.fasta'
                protein_adapter_3=protein_panel+"/"+str([f for f in os.listdir(protein_panel) if re.search(pattern_fasta3,f)][0])
                pattern_fasta5=r'\S*adapters_5\S*\.fasta'
                protein_adapter_5=protein_panel+"/"+str([f for f in os.listdir(protein_panel) if re.search(pattern_fasta5,f)][0])
                pattern_csv_file=r'\S*.csv'
                protein_barcodes=protein_panel+"/"+str([f for f in os.listdir(protein_panel) if re.search(pattern_csv_file,f)][0])
                
            except KeyError:
                protein_adapter_3='/mnt/beegfs/pipelines/single-cell_dna/tapestri_database/v2/panels/protein/ab_adapters_3.fasta'
                protein_adapter_5='/mnt/beegfs/pipelines/single-cell_dna/tapestri_database/v2/panels/protein/ab_adapters_5.fasta'
                protein_barcodes='/mnt/beegfs/pipelines/single-cell_dna/tapestri_database/v2/panels/protein/totalseq-d-heme-oncology.csv'
                
            prot_r1=str([f for f in design_file.iloc[row]["protein_file"].split(",") if re.search(pattern_r1,f)][0])
            prot_r2=str([f for f in design_file.iloc[row]["protein_file"].split(",") if re.search(pattern_r2,f)][0])
        
        dna_r1=str([f for f in design_file.iloc[row]["dna_file"].split(",") if re.search(pattern_r1,f)][0])
        dna_r2=str([f for f in design_file.iloc[row]["dna_file"].split(",") if re.search(pattern_r2,f)][0])
        
        if config_file["type_analysis"] == "dna+protein":
    
            dict_panel={'genome':{'path':reference_genome_path,
                                  'version':reference_genome},
                        'inputs':{'tube1':{'r1':[str(dna_r1)],'r2':[str(dna_r2)]}},
                        'protein_inputs':{'tube1':{'r1':[str(prot_r1)],'r2':[str(prot_r2)]}},
                        'protein_panel':{'adapters_3':protein_adapter_3,
                        'adapters_5':protein_adapter_5,
                        'barcodes':protein_barcodes},
                        'output':{'prefix':design_file.iloc[row]["sample_id"]},
                        'panel':{'name':panel_name_dna,
                                 'path':config_file["panel_path"]}}
        else:
            
            dict_panel={'genome':{'path':reference_genome_path,
                                  'version':reference_genome},
                        'inputs':{'tube1':{'r1':[dna_r1],
                                           'r2':[dna_r2]}},
                        'output':{'prefix':design_file.iloc[row]["sample_id"]},
                        'panel':{'name':panel_name_dna,
                                 'path':config_file["panel_path"]}}

        file_path=str(config_file["input_sample_path"])+"/"+str(design_file.iloc[row]["sample_id"])+"/"+str(design_file.iloc[row]["sample_id"])+"_config_panel.yaml"
        stream=open(file_path,'w')
        yaml.dump(dict_panel,stream,default_flow_style=False)

def cnv_select_barcodes_and_normalization(assay,barcodes_list):

    assay.cnv = assay.cnv[barcodes_list,:]

    assay.cnv.normalize_reads()

    assay.cnv.compute_ploidy(diploid_cells=barcodes_list)

    return(assay)
    
def make_label_for_each_variants(assay,list_of_variants):
    df=pd.DataFrame(assay.dna.layers['AF_MISSING'],index=assay.dna.barcodes(),columns=assay.dna.ids())
    for variants in list_of_variants:
        label_df= pd.DataFrame(df[variants].values,index=assay.dna.barcodes(),columns=[variants])
        value_to_label=df[variants].values
        labelling_tmp=[]
        split_label=variants.split(":")
        for i in range(len(value_to_label)):
            if value_to_label[i] > 20:
                labelling_tmp.append(split_label[0]+"+")
            else:
                labelling_tmp.append(split_label[0]+"-")

        assay.dna.row_attrs[split_label[0]+"_label"]=np.array(labelling_tmp)

def make_cartesian_product_of_variants(assay):
    """
    return a cartesian product from variants of interest 
    """
    list_for_cartesian_product=[list(np.unique(assay.dna.row_attrs[i])) for i in regex_search_pattern_in_list(r"^\w*_label$",assay.dna.row_attrs.keys())]
    
    return [cartproduct for cartproduct in itertools.product(*list_for_cartesian_product)]

def make_label(assay,list_of_cartesian_product_variants):
    """
    automatically return an array of label from variants of interest based on the different genotype possible
    """
    label=[]
    for i in range(len(assay.dna.barcodes())):
        label_row_to_compare=tuple([assay.dna.row_attrs[k][i] for k in regex_search_pattern_in_list(r"^\w*_label$",assay.dna.row_attrs.keys())])
        for j in list_of_cartesian_product_variants:
            if label_row_to_compare == j:
                label.append("_".join(j))
                
    return(np.array(label))