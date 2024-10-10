# Snakemake workflow: rna-seq-star-deseq2

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737358.svg)](https://doi.org/10.5281/zenodo.4737358)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/rna-seq-star-deseq2/workflows/Tests/badge.svg?branch=master)](https://github.com/snakemake-workflows/rna-seq-star-deseq2/actions?query=branch%3Amaster+workflow%3ATests)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

This workflow performs a differential gene expression analysis with STAR and Deseq2.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).


# Daylily Notes
- I've made some changes to thread assignment and set global defaults much higher.
- Using cookier cutter to create a slurm profile: `cookiecutter https://github.com/Snakemake-Profiles/slurm.git`
- modified the config/units.yaml file to point to the test RNA data
- make a cache dir on the shared FS, set the export SNAKEMAKE_OUTPUT_CACHE=/fsx/resources/environments/containers/ubuntu/cache/
- I make a conda snakemake env, `conda create -n snakmake -f snakemake_env.yaml`
- `snakemake --version` == 8.20.6
```bash
snakemake --use-conda --use-singularity -j 1  --singularity-prefix /fsx/resources/environments/containers/ubuntu/ip-10-0-0-240/ --singularity-args "  -B /tmp:/tmp -B /fsx:/fsx  -B /home/$USER:/home/$USER -B $PWD/:$PWD" --conda-prefix /fsx/resources/environments/containers/ubuntu/ip-10-0-0-240/ --profile ./profiles/slurm_pcluster3 --cache
```
