rule UMI_tools_qiagen:
    input:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/untrimmed/merged/{sample}-R1.fq.gz",
        ),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/umi_extract/{sample}-R1.fq.gz",
        ),
    params:
        bc_pattern=config.get("params").get("umi").get("qiagen_bc_pattern"),
        umi_length=config.get("params").get("umi").get("umi_length"),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/umi_tools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/umi_tools/{sample}_UMI_trimmed.log",
        ),
    shell:
        "umi_tools extract "
        "--stdin={input} "
        "--stdout={output} "
        "--log={log} "
        "--extract-method=regex "
        "--bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT)"
        "{{s<=2}}"
        "(?P<umi_1>.{{12}})(?P<discard_2>.*)'"


rule umi_deduplication_mature_qiagen:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/aligned/{sample}.mirna.bam"
        ),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/aligned/{sample}.mirna.bam.bai",
        ),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/dedup/{sample}.mirna_dedup.bam",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/umi_tools/{sample}_deduplication.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/umi_tools.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "umi_tools dedup "
        "-I {input.bam} "
        "-S {output} "
        "--method=unique "
        "--log={log} "


rule umi_deduplication_hairpin_qiagen:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/aligned/{sample}.mirbase_hairpin.sorted.bam",
        ),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/aligned/{sample}.mirbase_hairpin.sorted.bam.bai",
        ),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/dedup/{sample}.hairpin_dedup.bam",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/umi_tools/{sample}_hairpin_deduplication.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/umi_tools.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "umi_tools dedup "
        "-I {input.bam} "
        "-S {output} "
        "--method=unique "
        "--log={log} "
