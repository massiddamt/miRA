

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
            config.get("paths").get("workdir"), "workflow/envs/bash.yaml"
        )
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
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/merge_units.py"
        )


#
# rule rename_fastq_pre_trim:
#     input:
#         r1=resolve_results_filepath(
#             config.get("paths").get("results_dir"),"reads/umi_extract/{sample}-R1.fq.gz" if config.get("params").get("umi").get("included")=="Y" else "reads/untrimmed/merged/{sample}-R1.fq.gz")
#     output:
#         r1=resolve_results_filepath(
#             config.get("paths").get("results_dir"),"reads/untrimmed/{sample}.fq.gz",)
#     conda:
#         resolve_single_filepath(
#             config.get("paths").get("workdir"),"workflow/envs/bash.yaml"
#         )
#     shell:
#         "mv {input.r1} {output.r1}"
