rule samtools_sam_to_bam_mirbase1:
    input:
        sam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature.sam")
    output:
        bam=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature.bam"))
    params:
        genome=config.get("resources").get("mirna_mature_fa"),
        output_fmt="BAM"
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"benchmarks/samtools/sam_to_bam/{sample}.txt")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools view -b "
        "--threads {threads} "
        "-T {params.genome} "
        "-o {output.bam} "
        "-O {params.output_fmt} "
        "{input.sam} "

rule samtools_sam_to_bam_mirbase2:
    input:
        sam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature2.sam")
    output:
        bam=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature2.bam"))
    params:
        genome=config.get("resources").get("mirna_mature_fa"),
        output_fmt="BAM"
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"benchmarks/samtools/sam_to_bam/{sample}.txt")

    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools view -b "
        "--threads {threads} "
        "-T {params.genome} "
        "-o {output.bam} "
        "-O {params.output_fmt} "
        "{input.sam} "


rule samtools_sort_mirbase1:
    input:
        first=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature.bam")
    output:
        first=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature.sorted.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    params:
        tmp_dir=config.get("paths").get("tmp_dir"),
        genome=config.get("resources").get("mirna_mature_fa"),
        output_fmt="BAM"
    benchmark:
        "benchmarks/samtools/sort/{sample}.txt"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools sort "
        "--threads {threads} "
        "-T {params.tmp_dir} "
        "-O {params.output_fmt} "
        "--reference {params.genome} "
        "-o {output.first} "
        "{input.first} "



rule samtools_sort_mirbase2:
    input:
        second=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature2.bam")
    output:
        second=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature2.sorted.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    params:
        tmp_dir=config.get("paths").get("tmp_dir"),
        genome=config.get("resources").get("mirna_mature_fa"),
        output_fmt="BAM"
    benchmark:
        "benchmarks/samtools/sort/{sample}.txt"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools sort "
        "--threads {threads} "
        "-T {params.tmp_dir} "
        "-O {params.output_fmt} "
        "--reference {params.genome} "
        "-o {output.second} "
        "{input.second} "


rule samtools_merge_mirbase:
    input:
        first = resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature.sorted.bam"),
        second = resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_mature2.sorted.bam")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirna.bam")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        "benchmarks/samtools/merge/{sample}.txt"
    params:
        cmd='samtools',
        genome=config.get("resources").get("mirna_mature_fa"),
        output_fmt="BAM"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools merge "
        "--threads {threads} "
        "-O {params.output_fmt} "
        "--reference {params.genome} "
        "{output} "
        "{input.first} {input.second} "


rule index_mirbase:
    input:
        rules.samtools_merge_mirbase.output
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirna.bam.bai")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools index "
        "{input}"



### hairpin
rule samtools_sam_to_bam_hairpin:
    input:
        second=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_hairpin.sam")
    output:
        second=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_hairpin.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        "benchmarks/samtools/sam_to_bam/{sample}.txt"
    params:
        genome=config.get("resources").get("mirna_hairpin_fa"),
        output_fmt="BAM"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools view -b "
        "--threads {threads} "
        "-T {params.genome} "
        "-o {output.second} "
        "-O {params.output_fmt} "
        "{input.second} "


rule samtools_sort_hairpin:
    input:
        first=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_hairpin.bam")
    output:
        first=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_hairpin.sorted.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    params:
        tmp_dir=config.get("paths").get("tmp_dir"),
        genome=config.get("resources").get("mirna_hairpin_fa"),
        output_fmt="BAM"
    benchmark:
        "benchmarks/samtools/sort_hairpin/{sample}.txt"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools sort "
        "--threads {threads} "
        "-T {params.tmp_dir} "
        "-O {params.output_fmt} "
        "--reference {params.genome} "
        "-o {output.first} "
        "{input.first} "


rule index_hairpin:
    input:
        rules.samtools_sort_hairpin.output
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.mirbase_hairpin.sorted.bam.bai")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools index "
        "{input}"

### pirna
rule samtools_sam_to_bam_pirna:
    input:
        second=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.piRNA.sam")
    output:
        second=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.piRNA.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        "benchmarks/samtools/sam_to_bam/{sample}.pirna.txt"
    params:
        genome=config.get("resources").get("pirna_fa"),
        output_fmt="BAM"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools view -b "
        "--threads {threads} "
        "-T {params.genome} "
        "-o {output.second} "
        "-O {params.output_fmt} "
        "{input.second} "


rule samtools_sort_pirna:
    input:
        first=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.piRNA.bam")
    output:
        first=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.piRNA.sorted.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    params:
        tmp_dir=config.get("paths").get("tmp_dir"),
        genome=config.get("resources").get("pirna_fa"),
        output_fmt="BAM"
    benchmark:
        "benchmarks/samtools/sort/{sample}.pirna.txt"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools sort "
        "--threads {threads} "
        "-T {params.tmp_dir} "
        "-O {params.output_fmt} "
        "--reference {params.genome} "
        "-o {output.first} "
        "{input.first} "


rule index_pirna:
    input:
        rules.samtools_sort_pirna.output
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.piRNA.sorted.bam.bai")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools index "
        "{input}"

#### genome


rule samtools_sam_to_bam_genome:
    input:
        second=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.genome.sam")
    output:
        second=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.genome.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        "benchmarks/samtools/sam_to_bam/{sample}.genome.txt"
    params:
        genome=config.get("resources").get("reference"),
        output_fmt="BAM"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools view -b "
        "--threads {threads} "
        "-T {params.genome} "
        "-o {output.second} "
        "-O {params.output_fmt} "
        "{input.second} "


rule samtools_sort_genome:
    input:
        first=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.genome.bam")
    output:
        first=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.genome.sorted.bam"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    params:
        tmp_dir=config.get("paths").get("tmp_dir"),
        genome=config.get("resources").get("reference"),
        output_fmt="BAM"
    benchmark:
        "benchmarks/samtools/sort/{sample}.genome.txt"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools sort "
        "--threads {threads} "
        "-T {params.tmp_dir} "
        "-O {params.output_fmt} "
        "--reference {params.genome} "
        "-o {output.first} "
        "{input.first} "


rule index_genome:
    input:
        rules.samtools_sort_genome.output
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/aligned/{sample}.genome.sorted.bam.bai")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Generating BAI for reads mapped against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools index "
        "{input}"


rule index_deduplicated_bam:
    input:
        rules.umi_deduplication.output
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/dedup/{sample}.mirna_dedup.bam.bai")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "samtools index "
        "{input}"