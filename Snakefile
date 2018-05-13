import pandas as pd
import yaml
import jsonschema
shell.executable("bash")


config_schema = yaml.load(open("schemas/config.schema.yaml"))
samples_schema = yaml.load(open("schemas/samples.schema.yaml"))
units_schema = yaml.load(open("schemas/units.schema.yaml"))


configfile: "config.yaml"
jsonschema.validate(config, config_schema)

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
jsonschema.validate(samples.to_dict("records"), samples_schema)

units = pd.read_table(config["units"]).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
jsonschema.validate(units.to_dict("records"), units_schema)


def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])


rule all:
    input:
        expand("results/diffexp/{contrast}.diffexp.tsv",
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
