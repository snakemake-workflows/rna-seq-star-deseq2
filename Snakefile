import pandas as pd
from snakemake.utils import validate
shell.executable("bash")


configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")


def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])


rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


report: "report/workflow.rst"

include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
