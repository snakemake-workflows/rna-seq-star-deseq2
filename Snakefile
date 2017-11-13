import pandas as pd


configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col="sample")


rule all:
    input:
        expand("results/diffexp/{contrast}.diffexp.tsv",
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg"


include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
