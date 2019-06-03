rule count_matrix:
    input:
        expand("star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=units["sample"].tolist()
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

# https://twitter.com/mikelove/status/918770188568363008
def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

def get_deseq2_IA_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_interaction(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule deseq2_init:
    input:
        counts="counts/all.tsv"
    output:
        data = "deseq2/all.rds"
    params:
        samples=config["samples"],
        variables=config["diffexp"]["variables"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

rule pca:
    input:
        "deseq2/all.rds"
    output:
        report("results/pca.svg", "../report/pca.rst")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

def get_interaction(wildcards):
    return config["diffexp"]["interactions"][wildcards.interaction]

rule deseq2:
    input:
        "deseq2/all.rds"
    output:
        table=report("results/diffexp/contrasts/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/contrasts/{contrast}.ma-plot.svg", "../report/ma.rst")
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

rule deseq2_IA:
    input:
        "deseq2/all.rds"
    output:
        table=report("results/diffexp/interactions/{interaction}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/interactions/{interaction}.ma-plot.svg", "../report/ma.rst")
    params:
        interaction=get_interaction
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{interaction}.diffexp.log"
    threads: get_deseq2_IA_threads
    script:
        "../scripts/deseq2IA.R"
