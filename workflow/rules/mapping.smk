

rule mapping_miRBase_mature:
    input:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/trimmed/{sample}-trimmed.fq"
        ),
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.mirbase_mature.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.mirbase_mature.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "strict_params"
            if config.get("allow_multimapping") == "NO"
            else "strict_multimap_params"
        ),
        basename=config.get("resources").get("mirna_mature"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_mature.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against MIRBASE MATURE DB with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log}"


rule mapping_miRBase_hairpin:
    input:
        rules.mapping_miRBase_mature.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.mirbase_hairpin.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.mirbase_hairpin.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "strict_params"
            if config.get("allow_multimapping") == "NO"
            else "strict_multimap_params"
        ),
        basename=config.get("resources").get("mirna_hairpin"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_hairpin.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against MIRBASE PRECURSOSR DB with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log}"


rule mapping_piRNA:
    input:
        rules.mapping_miRBase_hairpin.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.piRNA.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.piRNA.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "strict_params"
            if config.get("allow_multimapping") == "NO"
            else "strict_multimap_params"
        ),
        basename=config.get("resources").get("pirna"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_pirna.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against PIRBASE MATURE DB with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log} "


rule mapping_tRNA:
    input:
        rules.mapping_piRNA.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.tRNA.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.tRNA.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "strict_params"
            if config.get("allow_multimapping") == "NO"
            else "strict_multimap_params"
        ),
        basename=config.get("resources").get("trna"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_trna.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against tRNA DB with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log}"


rule mapping_rRNA:
    input:
        rules.mapping_tRNA.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.rRNA.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.rRNA.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "strict_params"
            if config.get("allow_multimapping") == "NO"
            else "strict_multimap_params"
        ),
        basename=config.get("resources").get("rrna"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_rrna.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against rRNA DB with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log} "


rule mapping_mRNA:
    input:
        rules.mapping_rRNA.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.mRNA.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.mRNA.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "strict_params"
            if config.get("allow_multimapping") == "NO"
            else "strict_multimap_params"
        ),
        basename=config.get("resources").get("mrna"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_mrna.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against mRNA DB with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log} "


rule mapping_other:
    input:
        rules.mapping_mRNA.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.other_miRNA.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.other_miRNA.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "strict_params"
            if config.get("allow_multimapping") == "NO"
            else "strict_multimap_params"
        ),
        basename=config.get("resources").get("mirna_other"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_other.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against OTHER MIRBASE MATURE DB with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log} "


rule mapping_miRBase_mature2:
    input:
        rules.mapping_other.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.mirbase_mature2.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.mirbase_mature2.fastq",
            )
        ),
    params:
        params=config.get("params")
        .get("bowtie_mapping")
        .get(
            "mismatch_params"
            if config.get("allow_multimapping") == "NO"
            else "mismatch_multimap_params"
        ),
        basename=config.get("resources").get("mirna_mature"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_mature2.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against PIRBASE MATURE DB with 2 mismatches and {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log} "


rule mapping_genome:
    input:
        rules.mapping_miRBase_mature2.output.fastq,
    output:
        sam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/aligned/{sample}.genome.sam",
            )
        ),
        fastq=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "reads/unaligned/{sample}.genome.fastq",
            )
        ),
    params:
        params=config.get("params").get("bowtie_mapping").get("genome_params"),
        basename=config.get("resources").get("genome"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/bowtie/{sample}_mapping_genome.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bowtie.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Mapping reads against GENOME with {threads} threads for the following files {input}."
    shell:
        "bowtie "
        "{params.params} "
        "--threads {threads} "
        "{params.basename} "
        "{input} "
        "--un {output.fastq} "
        "-S {output.sam} "
        ">& {log} "
