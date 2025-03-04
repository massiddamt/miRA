# Snakemake workflow: miRA - miRNA Analysis
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.27.1-brightgreen.svg)](https://snakemake.bitbucket.io)

This workflow performs miRNA analysis from `fastq` to `counts`.

## Authors

* [Matteo Massidda](https://github.com/massiddamt), Institute for Genetic and Biomedical Research (IRGB) - National Research Council (CNR)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=massiddamt/miRA).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## INSTRUCTIONS
Create a virtual environment with the command:
```commandline
conda create -c conda-forge -c bioconda --name snakemake snakemake=8.27.1 conda=24.7.1 snakedeploy
```
and activate it:
```commandline
conda activate snakemake
```
You can perform the pipeline deploy defining a directory `my_dest_dir` for analysis output and a pipeline tag for a specific version:
```bash
snakedeploy deploy-workflow https://github.com/massiddamt/miRA.git 
                    my_desd_dir 
                    --tag v1.0.2
```
To run the pipeline, go inside the deployed pipeline folder and use the command:
```bash
snakemake --use-conda -p --cores all
```