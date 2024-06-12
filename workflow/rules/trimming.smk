




rule pre_rename_fastq:
    input:
        r1=lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair="fq1")
    output:
        r1=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/untrimmed/{sample}.fq.gz",)
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/bash.yaml"
        )
    shell:
        "cp {input.r1} {output.r1}"


rule trimming:
    input:
        read1=rules.UMI_tools.output
        #resolve_results_filepath(
         #   config.get("paths").get("results_dir"),"reads/untrimmed/{sample}.fq.gz")
    output:
        read1=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/trimmed/{sample}_trimmed.fq",
            )
        ),
        read1_trimming_report=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/trimmed/{sample}.fq.gz_trimming_report.txt",
        )
    params:
        extra=config.get("params").get("trim_galore").get("arguments"),
        outdir=lambda wildcards, output: os.path.dirname(output.read1),
        qc_dir=resolve_results_filepath(
            config.get("paths").get("results_dir"), "qc/fastqc"
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "logs/trim_galore/{sample}.log"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/trim_galore.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Trimming reads with TRIM_GALORE with {threads} threads for the following files {input.read1}."
    shell:
        "mkdir -p {params.qc_dir}; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input.read1} "
        ">& {log} "


rule rename_trimmed_fastq:
    input:
        read1=rules.trimming.output.read1
    output:
        read1=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/trimmed/{sample}-trimmed.fq")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "logs/bash/{sample}.log"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bash.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Rename fastq files {input.read1}."
    shell:
        "mv {input.read1} {output.read1} "
        ">& {log}"