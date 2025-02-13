rule trimming:
    input:
        read1=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            (
                "reads/umi_extract/{sample}-R1.fq.gz"
                if config.get("params").get("umi").get("included") == "Y"
                else "reads/untrimmed/merged/{sample}-R1.fq.gz"
            ),
        ),
    output:
        read1=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/trimmed/{sample}-R1_trimmed.fq",
            )
        ),
        read1_trimming_report=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt",
        ),
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
        "workflow/envs/trim_galore.yaml"
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
        read1=rules.trimming.output.read1,
    output:
        read1=resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/trimmed/{sample}-trimmed.fq"
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "logs/bash/{sample}.log"
        ),
    conda:
        "workflow/envs/bash.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Rename fastq files {input.read1}."
    shell:
        "mv {input.read1} {output.read1} "
        ">& {log}"
