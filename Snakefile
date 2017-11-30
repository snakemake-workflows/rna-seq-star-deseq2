import pandas as pd
shell.executable("bash")


configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index


def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])


rule all:
    input:
        expand("results/diffexp/{contrast}.diffexp.tsv",
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg"


include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
