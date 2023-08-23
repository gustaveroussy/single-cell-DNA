try:
    normalisation_list=config["prot_norm_dimred"]["normalization"]
except KeyError:
        normalisation_list=["CLR"]

if "DSB" in normalisation_list:
    rule save_in_csv_prot_counts:
        input:
            sample=config["output_sample_path"]+"/{i_sample}/h5/dna/filtering/QC_{i_sample}.h5"
        output:
            read_counts_path_to_save_output=config["output_sample_path"]+"/{i_sample}/prot/normalization/reads_counts/{i_sample}_read_counts.tsv",
        conda:
            CONDA_MOSAIC_ENV
        params:
            workflow_dir=PIPELINE_FOLDER,
            config_file='"'+str(config.copy())+'"'
        threads:
            1
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 2048, 10240)),
            time_min = (lambda wildcards, attempt: min(attempt * 5, 60))
        shell: 
            """
            python3 {params.workflow_dir}/scripts/save_counts_prot.py \
            {input.sample} \
            {output.read_counts_path_to_save_output} \
            {params.workflow_dir} \
            {params.config_file}
            """
    # to add later config["output_sample_path"]+"/{i_sample}/prot/diagnostic/:/WORKDIR/results 
    rule run_R_dsb_script:
        input:
            input_file_path_1=config["output_sample_path"]+"/{i_sample}/prot/normalization/reads_counts/{i_sample}_read_counts.tsv"
        output:
            #hist_path_1=config["output_sample_path"]+"/{i_sample}/prot/diagnostic/hist_{i_sample}_read_counts.png",
            output_file_path_1=config["output_sample_path"]+"/{i_sample}/prot/normalization/reads_counts/{i_sample}_read_counts_dsb_normalization.csv"
        params:
            workflow_dir=PIPELINE_FOLDER,
            sing_bind= config["output_sample_path"]+"/{i_sample}/prot/normalization/reads_counts/:/WORKDIR/normalization/reads_counts -B" +str(PIPELINE_FOLDER)+"/scripts:/WORKDIR/scripts",
            sing_img=SING_IMG_DSB,
            input_file_path="/WORKDIR/normalization/reads_counts/{i_sample}_read_counts.tsv",
            #hist_path="/WORKDIR/results/{i_sample}_read_counts.png",
            output_file_path="/WORKDIR/normalization/reads_counts/{i_sample}_read_counts_dsb_normalization.csv",
            isotype_control=config["prot_norm_dimred"]["isotype_control"]
        threads:
            3
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 2048, 10240)),
            time_min = (lambda wildcards, attempt: min(attempt * 10, 60))
        shell:
            """
            singularity exec --no-home -B {params.sing_bind} \
            {params.sing_img} \
            Rscript /WORKDIR/scripts/dsb_normalization.R \
            --input.file {params.input_file_path} \
            --output.file {params.output_file_path} \
            --isotype.control {params.isotype_control}
            """
            
    rule save_dsb_counts:
        input:
            counts=config["output_sample_path"]+"/{i_sample}/prot/normalization/reads_counts/{i_sample}_read_counts_dsb_normalization.csv",
            sample=config["output_sample_path"]+"/{i_sample}/h5/dna/filtering/QC_{i_sample}.h5"
        output:
            new_sample=config["output_sample_path"]+"/{i_sample}/h5/prot/norm_QC.h5"
        conda:
            CONDA_MOSAIC_ENV
        params:
            workflow_dir=PIPELINE_FOLDER,
            config_file='"'+str(config.copy())+'"'
        threads:
            1
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 2048, 10240)),
            time_min = (lambda wildcards, attempt: min(attempt * 10, 60))
        shell:
            """
            python3 {params.workflow_dir}/scripts/dsb_save.py {input.sample} {input.counts} {output.new_sample} {params.config_file}
            """
    rule prot_make_norm:
        input:
            sample=config["output_sample_path"]+"/{i_sample}/h5/dna/filtering/QC_{i_sample}.h5"
        output:
            new_sample=config["output_sample_path"]+"/{i_sample}/h5/prot/norm_QC_V2.h5",
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
            time_min = (lambda wildcards, attempt: min(attempt * 45, 60))
        shell: 
            """
            python3 {params.workflow_dir}/scripts/prot_normalization.py \
            --input_h5 {input.sample} \
            --output_h5 {output.new_sample} \
            --sample_name {params.sample_name} \
            --config_file {params.config_file}
            """
else:
    rule prot_make_norm:
        input:
            sample=config["output_sample_path"]+"/{i_sample}/h5/dna/filtering/QC_{i_sample}.h5"
        output:
            new_sample=config["output_sample_path"]+"/{i_sample}/h5/prot/norm_QC_V2.h5",
        conda:
            CONDA_MOSAIC_ENV
        params:
            workflow_dir=PIPELINE_FOLDER,
            config_file='"'+str(config.copy())+'"',
            sample_name="{i_sample}"
        threads:
            1
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
            time_min = (lambda wildcards, attempt: min(attempt * 5, 60))
        shell: 
            """
            python3 {params.workflow_dir}/scripts/prot_normalization.py \
            --input_h5 {input.sample} \
            --output_h5 {output.new_sample} \
            --sample_name {params.sample_name} \
            --config_file {params.config_file}
            """
rule prot_norm_dimred:
    input:
        sample=config["output_sample_path"]+"/{i_sample}/h5/prot/norm_QC_V2.h5"
    output:
        new_sample=config["output_sample_path"]+"/{i_sample}/h5/prot/dimred_QC_max_n_dim.h5"
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
        time_min = (lambda wildcards, attempt: min(attempt * 45, 120))
    shell: 
        """
        python3 {params.workflow_dir}/scripts/sc_dna_prot.py \
        --input_h5 {input.sample} \
        --output_h5 {output.new_sample} \
        --sample_name {params.sample_name} \
        --config_file {params.config_file} \
        """