import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")
shell.executable("/bin/bash")

##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.deseq2_ALL.get(**wildcards).output[0]
    return expand(["results/diffexpAll/{i}.diffexp.tsv",
            "results/diffexpAll/{i}.ma-plot.svg"],
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.ma-plot.svg')).i)

def aggregate_input2(wildcards):
    checkpoint_output = checkpoints.deseq2_ALL.get(**wildcards).output[0]
    return expand(["results/diffexpAll2/{i}.diffexp.tsv",
            "results/diffexpAll2/{i}.ma-plot.svg"],
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.ma-plot.svg')).i)

##### target rules #####

rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["advanced"]["contrasts"]),
        aggregate_input,
        aggregate_input2,
        "results/pca.svg",
        "qc/multiqc_report.html"


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
include: "rules/qc.smk"
