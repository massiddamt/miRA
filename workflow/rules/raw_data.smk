rule fastq_merge:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair="fq1"),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/untrimmed/merged/{sample}-R1.fq.gz",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bash.yaml")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/merge_fastqs/{sample}_fq_merge.log",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "merge fastq files {input}."
    script:
        "workflow/scripts/merge_units.py"
