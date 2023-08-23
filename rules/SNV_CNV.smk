try:
    method=config["snv_norm_dimred"]["method_dimred"]
except KeyError:
    method="pca"
try:
    max_components=config["snv_norm_dimred"]["max_dims"]
except KeyError:
    max_components=6
    
rule snv_cnv_dimred:
    input:
        sample=config["output_sample_path"]+"/{i_sample}/h5/dna/filtering/QC_{i_sample}.h5"
    output:
        new_sample=config["output_sample_path"]+"/{i_sample}/h5/dna/dimred_QC_max_n_dim_"+str(max_components)+".h5"
    conda:
        CONDA_MOSAIC_ENV
    params:
        workflow_dir=PIPELINE_FOLDER,
        config_file='"'+str(config.copy())+'"',
        sample_name="{i_sample}"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 2048, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 120, 360))
    shell: 
        """
        python3 {params.workflow_dir}/scripts/sc_dna_snv_cnv.py \
        --input_h5 {input.sample} \
        --output_h5 {output.new_sample} \
        --sample_name {params.sample_name} \
        --config_file {params.config_file} \
        """
