import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####


configfile: "config.yaml"


validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")


##### target rules #####


rule all:
    input:
        expand(
            [
                "results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg",
            ],
            contrast=config["diffexp"]["contrasts"],
        ),
        "results/pca.svg",
        "qc/multiqc_report.html",


##### setup container #####


# this container image defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##### setup report #####


report: "report/workflow.rst"


##### load rules #####


include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
include: "rules/qc.smk"
