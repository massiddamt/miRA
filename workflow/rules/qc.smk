rule fastqc:
    input:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/untrimmed/{sample}.fq.gz")
    output:
        html=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/untrimmed_{sample}.html"),
        zip=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/untrimmed_{sample}_fastqc.zip")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/fastqc/untrimmed/{sample}.log")
    params:
        outdir=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/fastqc.yaml")
    shell:
        "fastqc "
        "{input} "
        "--outdir {params.outdir} "
        "--quiet "
        ">& {log}"


rule fastqc_trimmed:
    input:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-trimmed.fq")
    output:
        html=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample}-trimmed_fastqc.html"),
        zip=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample}-trimmed_fastqc.zip")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/fastqc/trimmed/{sample}.log")
    params:
        outdir=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/fastqc.yaml")
    shell:
        "fastqc "
        "{input} "
        "--outdir {params.outdir} "
        "--quiet "
        ">& {log}"


rule mirtrace:
    input:
        rules.rename_trimmed_fastq.output
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample}/mirtrace-results.json"),
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample}/mirtrace-stats-length.tsv"),
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample}/mirtrace-stats-contamination_basic.tsv"),
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample}/mirtrace-stats-mirna-complexity.tsv")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/mirtrace.yaml")
    params:
        outdir=resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample}"),
        params="--force",
        species=config.get("params").get("mir_trace").get("arguments")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/mirtrace/{sample}.mirtrace_qc.log")
    shell:
        "mirtrace "
        "qc "
        "{params.species} "
        "-o {params.outdir} "
        "{input} "
        "{params.params} "
        ">& {log} "


rule multiqc:
    input:
        expand(resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample.sample}_fastqc.zip"), sample=samples.reset_index().itertuples()),
#        expand("qc/fastqc/{sample.sample}_umi_fastqc.zip", sample=samples.reset_index().itertuples()),
        expand(resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample.sample}-trimmed_fastqc.zip"), sample=samples.reset_index().itertuples()),
        expand(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample.sample}-R1.fq.gz_trimming_report.txt"), sample=samples.reset_index().itertuples()),
        expand(resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample.sample}/mirtrace-results.json"), sample=samples.reset_index().itertuples()),
        expand(resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample.sample}/mirtrace-stats-length.tsv"), sample=samples.reset_index().itertuples()),
        expand(resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample.sample}/mirtrace-stats-contamination_basic.tsv"), sample=samples.reset_index().itertuples()),
        expand(resolve_results_filepath(
            config.get("paths").get("results_dir"),"mir_trace/{sample.sample}/mirtrace-stats-mirna-complexity.tsv"), sample=samples.reset_index().itertuples()),
#        expand("htseq/{sample.sample}.counts", sample=samples.reset_index().itertuples())

    output:
        report(resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/multiqc.html"),caption="../report/multiqc.rst", category="QC")
    params:
        params=config.get("params").get("multiqc").get("arguments"),
        outdir=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc"),
        outname="multiqc.html",
        fastqc=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/"),
        trimming=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/"),
        reheader=config.get("reheader")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/multiqc.yaml")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/multiqc/multiqc.log")
    shell:
        "multiqc "
        "{input} "
        "{params.fastqc} "
        "{params.trimming} "
        "{params.params} "
        "-o {params.outdir} "
        "-n {params.outname} "
        "--sample-names {params.reheader} "
        ">& {log}"