import pandas as pd


configfile: "config.yaml"
samples = pd.read_table("samples.tsv", index_col=0)


rule all:
    input:
        expand("results/diffexp/{contrast}.diffexp.tsv",
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg"


include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
