rule count_matrix:
    input:
        expand("results/star/{unit.sample}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "results/counts/all.tsv"
    params:
        units=units,
        ## Is the sequencing data strand-specific? Can be no, yes, reverse. Default is 'no'.
        ## Check https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
        strand=config["params"]["strand-specific"]
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule add_metainfo:
    input:
        "{table}.tsv"
    output:
        "{table}.info.tsv"
    params:
        genetype=config["ref"]["metainfo"]
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/add_metainfo.py"


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        counts="results/counts/all.tsv"
    output:
        "results/deseq2/all.rds"
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "results/deseq2/all.rds"
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


rule deseq2:
    input:
        "results/deseq2/all.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

