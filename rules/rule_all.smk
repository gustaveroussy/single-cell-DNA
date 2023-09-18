import os

def get_targets(steps):
    target = {}

    try:
        snv_pca_components=config["snv_norm_dimred"]["max_dims"]
    except KeyError:
        snv_pca_components=6
    try:
        cnv_pca_components=config["cnv_norm_dimred"]["pca_max_dims"]
    except KeyError:
        cnv_pca_components=6
    try:
        normalisation_list=config["prot_norm_dimred"]["normalization"]
    except KeyError:
        normalisation_list=["CLR"]

    if "Alignment" in steps:
        if config["type_analysis"] == "dna":

            target["Alignment"]=[
                expand(config["input_sample_path"]+"/{i_sample}/{i_sample}_config_panel.yaml",i_sample=sample_list),
                expand(config["output_sample_path"]+"/{i_sample}/results/{i_sample}.dna.h5",i_sample=sample_list)
            ]
        if config["type_analysis"] == "dna+protein":
                target["Alignment"]=[
                expand(config["input_sample_path"]+"/{i_sample}/{i_sample}_config_panel.yaml",i_sample=sample_list),
                expand(config["output_sample_path"]+"/{i_sample}/{i_sample}.dna+protein.h5",i_sample=sample_list)
            ]

    if "filtering" in steps:

        target["filtering"]=[
        expand(config["output_sample_path"]+"/{i_sample}/h5/dna/filtering/QC_{i_sample}.h5",i_sample=sample_list),
        ]

    if "SNV_CNV" in steps:

        target["SNV_CNV"]=[
        expand(config["output_sample_path"]+"/{i_sample}/h5/dna/dimred_QC_max_n_dim_"+str(snv_pca_components)+".h5",i_sample=sample_list),
        ]
        
    if "PROTEIN" in steps:
    
        if "DSB" in config["prot_norm_dimred"]["normalization"]:
        
            target["PROTEIN"]=[
            expand(config["output_sample_path"]+"/{i_sample}/prot/normalization/reads_counts/{i_sample}_read_counts.tsv",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/prot/normalization/reads_counts/{i_sample}_read_counts_dsb_normalization.csv",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/h5/prot/norm_QC.h5",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/h5/prot/norm_QC_V2.h5",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/h5/prot/dimred_QC_max_n_dim.h5",i_sample=sample_list),
            ]
            
        else:
        
            target["PROTEIN"]=[
            expand(config["output_sample_path"]+"/{i_sample}/h5/prot/dimred_QC_max_n_dim.h5",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/h5/prot/norm_QC_V2.h5",i_sample=sample_list)
            ]
        
    if "ALL" in steps:
        
        target["ALL"]=[
        expand(config["output_sample_path"]+"/{i_sample}/h5/all/all_assays.h5",i_sample=sample_list)
        ]
    
    
    if "phylogeny" in steps:
        
        if len(config["phylogeny"]["methods"]) == 1 and "COMPASS" in config["phylogeny"]["methods"]:
            target["phylogeny"]=[
            expand(config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/whitelist.csv",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/output",i_sample=sample_list)
            ]
            
        if len(config["phylogeny"]["methods"]) == 1 and "infSCITE" in config["phylogeny"]["methods"]:
            target["phylogeny"]=[
            expand(config["output_sample_path"]+"/{i_sample}/dna/infscite/input/",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/dna/infscite/secondary_output/",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/dna/infscite/output/",i_sample=sample_list)
            ]
            
        if len(config["phylogeny"]["methods"]) == 2 and "infSCITE" in config["phylogeny"]["methods"]:
            target["phylogeny"]=[
            expand(config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/whitelist.csv",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/output",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/dna/infscite/input/",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/dna/infscite/secondary_output/",i_sample=sample_list),
            expand(config["output_sample_path"]+"/{i_sample}/dna/infscite/output/",i_sample=sample_list)
            ]
    return(target)