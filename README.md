# Snakemake workflow: rna-seq-star-deseq2

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.1.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2)

This workflow performs a differential expression analysis with STAR and Deseq2.
It is currently under development. No stable release is available yet.

## Authors

* Johannes Köster (@johanneskoester), https://koesterlab.github.io

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/rna-seq-spew/releases).
If you intend to modify and further develop this workflow, fork this reposity. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, once available, its DOI.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml` and the sheets `samples.tsv` and `units.tsv`.

### Step 3: Execute workflow

All you need to execute this workflow is to install Snakemake via the [Conda package manager](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda). Software needed by this workflow is automatically deployed into isolated environments by Snakemake.

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores. Alternatively, it can be run in cluster of cloud environments (see [the docs](http://snakemake.readthedocs.io/en/stable/executable.html) for details.
