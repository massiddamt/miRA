rule count_mirbase_mature:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/dedup/{sample}.mirna_dedup.bam"
            if config.get("params").get("umi").get("included") == "Y"
            else "reads/aligned/{sample}.mirna.bam",
        ),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/dedup/{sample}.mirna_dedup.bam.bai"
            if config.get("params").get("umi").get("included") == "Y"
            else "reads/aligned/{sample}.mirna.bam.bai",
        ),
    output:
        counts=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/counts/{sample}.mirna_mature.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/samtools/idxstats/{sample}_mature_counts.log",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "samtools idxstats {input.bam} | cut -f 1,3 "
        "> {output.counts} "


rule count_mirbase_hairpin:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/dedup/{sample}.hairpin_dedup.bam"
            if config.get("params").get("umi").get("included") == "Y"
            else "reads/aligned/{sample}.mirbase_hairpin.sorted.bam",
        ),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/dedup/{sample}.hairpin_dedup.bam.bai"
            if config.get("params").get("umi").get("included") == "Y"
            else "reads/aligned/{sample}.mirbase_hairpin.sorted.bam.bai",
        ),
    output:
        counts=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/counts/{sample}.mirna_hairpin.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/samtools/idxstats/{sample}_hairpin_counts.log",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "samtools idxstats {input.bam} | cut -f 1,3 "
        "> {output.counts} "


rule count_pirna:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/aligned/{sample}.piRNA.sorted.bam",
        ),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/aligned/{sample}.piRNA.sorted.bam.bai",
        ),
    output:
        counts=resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/counts/{sample}.piRNA.txt"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/samtools/idxstats/{sample}_pirna_counts.log",
        ),
    benchmark:
        "benchmarks/samtools/count/{sample}.pirna.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "samtools idxstats {input.bam} | cut -f 1,3 "
        "> {output.counts} "


rule htseq_genome:
    input:
        bam=rules.samtools_sort_genome.output,
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/aligned/{sample}.genome.sorted.bam.bai",
        ),
    output:
        counts=resolve_results_filepath(
            config.get("paths").get("results_dir"), "htseq/{sample}.counts"
        ),
    params:
        params=config.get("params").get("htseq").get("params"),
        gff=config.get("resources").get("htseq_gff"),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/htseq.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/htseq/{sample}_genome_counts.log",
        ),
    shell:
        "htseq-count "
        "{params.params} "
        "-q {input.bam} "
        "{params.gff} "
        ">{output.counts}"


rule get_counts_all:
    input:
        rules.count_mirbase_mature.output.counts,
        rules.count_mirbase_hairpin.output.counts,
        rules.count_pirna.output.counts,
        rules.htseq_genome.output.counts,
    output:
        touch(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/counts/{sample}.counts.completed",
            )
        ),
