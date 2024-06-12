def get_unit_fastqs(wildcards, samples, label='units',read_pair='fq'):
    for unit_set in samples.loc[wildcards.sample,[label]]:
        print(wildcards.sample)
    return [units.loc[x,[read_pair]].dropna()[0] for x in unit_set.split(',')]

rule fastq_merge:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair='fq1')
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/untrimmed/merged/{sample}-R1.fq.gz")
    script:
        "scripts/merge_units.py"


#def get_fastq(wildcards,samples,read_pair='fq'):
#    return samples.loc[wildcards.sample,
#                     [read_pair]].dropna()[0]

rule UMI_tools:
    input:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/untrimmed/merged/{sample}-R1.fq.gz")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/umi_extract/{sample}_umi.fq.gz")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/umi_tools.yaml")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/umi_tools/{sample}_UMI_trimmed.log")
    shell:
        "umi_tools extract "
        "--stdin={input} "
        "--stdout={output} "
        "--log={log} "
        "--extract-method=regex "
        "--bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT)"
        "{{s<=2}}"
        "(?P<umi_1>.{{12}})(?P<discard_2>.+)'"


rule umi_deduplication:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirna.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirna.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/dedup/{sample}.mirna_dedup.bam")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/umi_tools/{sample}_deduplication.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/umi_tools.yaml")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "umi_tools dedup "
        "-I {input.bam} "
        "-S {output} "
        "--method=unique "
        "--log={log} "