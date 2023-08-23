try:
    max_snv_dim=config["snv_norm_dimred"]["max_dims"]
except KeyError:
    max_snv_dim=6

if config["type_analysis"] == "dna+protein":
    sample=config["output_sample_path"]+"/{i_sample}/h5/prot/dimred_QC_max_n_dim.h5"
if config["type_analysis"] == "dna":
    sample=config["output_sample_path"]+"/{i_sample}/h5/dna/dimred_QC_max_n_dim_"+str(max_snv_dim)+".h5"

rule all_assays:
    input:
        sample=sample
    output:
        new_sample=config["output_sample_path"]+"/{i_sample}/h5/all/all_assays.h5"
    conda:
        CONDA_MOSAIC_ENV
    params:
        workflow_dir=PIPELINE_FOLDER,
        config_file='"'+str(config.copy())+'"',
        sample_name="{i_sample}"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 3072, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 360))
    shell: 
        """
        python3 {params.workflow_dir}/scripts/sc_all_assays.py \
        --input_h5 {input.sample} \
        --output_h5 {output.new_sample} \
        --sample_name {params.sample_name} \
        --config_file {params.config_file} \
        """
