if "COMPASS" in config["phylogeny"]["methods"]:
    try:
        bool_cnv=config["phylogeny"]["COMPASS"]["bool_cnv"]
    except KeyError:
        bool_cnv=0
    
    if config["type_analysis"] == "dna+protein":
        loom_file=config["output_sample_path"]+"/{i_sample}/dna/results/{i_sample}.cells.loom"
    if config["type_analysis"] == "dna":
        loom_file=config["output_sample_path"]+"/{i_sample}/results/{i_sample}.cells.loom"
        
    rule prepare_csv:
        input:
            sample=config["output_sample_path"]+"/{i_sample}/h5/all/all_assays.h5",
            csv=config["output_sample_path"]+"/{i_sample}/dna/annotation/QC_annotation.csv"
        output:
            whitelist_csv=config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/whitelist.csv"
        conda:
            CONDA_MOSAIC_ENV
        params:
            workflow_dir=PIPELINE_FOLDER,
            config_file='"'+str(config.copy())+'"'
        threads:
            1
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
            time_min = (lambda wildcards, attempt: min(attempt * 5, 60))
        shell:
            "python3 {params.workflow_dir}/scripts/sc_compass_variant_to_csv.py {input.sample} {input.csv} {output.whitelist_csv} {params.config_file}"
        
    rule preprocess:
        input:
            sample=loom_file,
            whitelist=config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/whitelist.csv"
        output:
            output_path=directory(config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/preprocessed/")
        conda:
            PIPELINE_COMPASS_PREPROCESS
        params:
            workflow_dir=PIPELINE_FOLDER,
            config_file='"'+str(config.copy())+'"'
        threads:
            1
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
            time_min = (lambda wildcards, attempt: min(attempt * 5, 60))
        shell:
            """
            pyensembl install --release 75 76
            mkdir -p -v {output.output_path}
            python3 {params.workflow_dir}/scripts/sc_compass_preprocess.py \
            -i {input.sample} \
            -o {output.output_path} \
            --whitelist {input.whitelist}
            """
      
    rule run_compass:
        input:
            sample_path=config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/preprocessed/"
        output:
            output_path_dir=directory(config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/output")
        params:
            singularity_env=SING_IMG,
            sing_grp_bind =config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/preprocessed:/mnt/data -B "+config["output_sample_path"]+"/{i_sample}/dna/compass/QC_annotation/output:/mnt/output -B "+str(PIPELINE_FOLDER)+"/scripts:/mnt/scripts",
            chain_length=str(4000),
            nb_core=str(2),
            cnv_bool=str(bool_cnv)
        threads:
            3
        resources:
            mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
            time_min = (lambda wildcards, attempt: min(attempt * 30, 240))
        shell:
            """
            mkdir -p -v {output.output_path_dir} 
            singularity exec --no-home -B {params.sing_grp_bind} \
            {params.singularity_env} \
            bash /mnt/scripts/sc_compass_run.sh \
            /mnt/data/{wildcards.i_sample}.cells \
            /mnt/output/{wildcards.i_sample} \
            {params.nb_core} \
            {params.chain_length} \
            {params.cnv_bool} \
            /mnt/output/{wildcards.i_sample}_tree.png \
            /mnt/output/{wildcards.i_sample}_tree.gv
            """
