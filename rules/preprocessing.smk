# filtering low quality variants, plotting diagnostic plots and elbow plot for reduction dimension

if config["type_analysis"] == "dna":
    sample_input=config["output_sample_path"]+"/{i_sample}/results/{i_sample}.dna.h5"
if config["type_analysis"] == "dna+protein":
    sample_input=config["output_sample_path"]+"/{i_sample}/{i_sample}.dna+protein.h5"

rule preprocessing_snv:
    input:
        sample=sample_input
    output:
        h5_filtering = config["output_sample_path"]+"/{i_sample}/h5/dna/filtering/QC_{i_sample}.h5",
    conda:
        CONDA_MOSAIC_ENV
    params:
        workflow_dir=PIPELINE_FOLDER,
        sample_name="{i_sample}",
        config_file='"'+str(config.copy())+'"'
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 2048, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 60))
    shell:
        """
        python3 \
        {params.workflow_dir}/scripts/sc_dna_preprocessing.py \
        --input_h5 {input.sample} \
        --ouput_h5 {output.h5_filtering} \
        --sample_name {params.sample_name} \
        --config_file {params.config_file}
        """