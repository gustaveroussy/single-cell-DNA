rule create_panel:
    output:
        panel_name=config["input_sample_path"]+"/{i_sample}/{i_sample}_config_panel.yaml"
    threads:
        1
    conda:
        CONDA_MOSAIC_ENV
    params:
        workflow_dir=PIPELINE_FOLDER,
        config_file='"'+str(config.copy())+'"'
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 2048)),
        time_min = (lambda wildcards, attempt: min(attempt * 10, 60))
    shell:
        """
        python3 {params.workflow_dir}/scripts/create_panel.py {params.config_file}
        """
        
        
if config["type_analysis"] == "dna":    
    rule alignment_dna:
        input:
            configfile_i=config["input_sample_path"]+"/{i_sample}/{i_sample}_config_panel.yaml"
        output:
            file_output=config["output_sample_path"]+"/{i_sample}/results/{i_sample}.dna.h5"
        threads:
            24
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 61440, 122880)),
            time_min = (lambda wildcards, attempt: min(attempt * 4320, 8640))
        params:
            workflow_dir=PIPELINE_FOLDER,
            file_output=config["output_sample_path"]+"/{i_sample}"
        shell:
            """
            bash {params.workflow_dir}/scripts/run_dna_alignment.sh {params.file_output} {input.configfile_i}
            """
        
if config["type_analysis"] == "dna+protein":    
    rule alignment_dna_protein:
        input:
            configfile_i=config["input_sample_path"]+"/{i_sample}/{i_sample}_config_panel.yaml"
        output:
            file_output=config["output_sample_path"]+"/{i_sample}/{i_sample}.dna+protein.h5"
        threads:
            24
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 61440, 122880)),
            time_min = (lambda wildcards, attempt: min(attempt * 4320, 8640))
        params:
            workflow_dir=PIPELINE_FOLDER,
            file_output=config["output_sample_path"]+"/{i_sample}"
        shell:
            """
            bash {params.workflow_dir}/scripts/run_alignment_dna_protein.sh {params.file_output} {input.configfile_i}
            """